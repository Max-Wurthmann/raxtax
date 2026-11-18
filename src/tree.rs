use anyhow::Result;

use indicatif::HumanBytes;
use itertools::Itertools;
use log::{debug, log_enabled, Level};
use logging_timer::time;

use crate::{
    parser::{LabeledSequence, LineageBinPair},
    utils::{seq_to_unique_minenc_canon_kmers_into, KMerEncodingData},
};

#[cfg(feature = "huge_db")]
pub type IndexType = usize;

#[cfg(not(feature = "huge_db"))]
pub type IndexType = u32;

#[cfg(not(feature = "huge_db"))]
fn check_lineage_size(db_size: usize) {
    assert!(
        u32::try_from(db_size).is_ok(),
        "Too many database sequences to run with 32-bit indices!\n
            Re-compile raxtax with '--features huge_db' to enable usize indices."
    );
}

#[cfg(feature = "huge_db")]
fn check_lineage_size(_: usize) {}

/// The taxonomy of the reference database. Reference sequences are indexed in
/// lineage-sorted order, so that every `Node::confidence_range` covers a contiguous
/// range of them.
#[derive(Debug, Eq, PartialEq)]
pub struct Tree {
    pub root: Node,
    pub lineages: Vec<String>,
    pub bins: Vec<String>,

    // maps a reference sequence's position in the reference file to its (lineage-sorted) index
    pub sorted_idx: Vec<IndexType>,

    pub encoding_data: KMerEncodingData,
    pub bin_idx_to_lineage_idxs: Vec<Vec<usize>>,
    pub lineage_idx_to_bin_idx: Vec<Option<usize>>,
    pub num_tips: usize,
}

impl Tree {
    /// Builds the taxonomy from the reference labels, given in reference file order.
    #[time("debug", "Tree::{}")]
    pub fn new(mut labels: Vec<LineageBinPair>, encoding_data: KMerEncodingData) -> Result<Self> {
        check_lineage_size(labels.len());
        let mut root = Node::new(String::from("root"), 0, NodeType::Inner);

        // stable sort of the file positions by label, so the sequences can be streamed in file
        // order later and mapped to their sorted index via `sorted_idx`
        let mut order = (0..labels.len()).collect_vec();
        order.sort_by(|&a, &b| labels[a].cmp(&labels[b]));
        let mut sorted_idx = vec![IndexType::default(); labels.len()];
        for (idx, &file_idx) in order.iter().enumerate() {
            sorted_idx[file_idx] = idx as IndexType;
        }
        let sorted_lineages = order
            .into_iter()
            .map(|file_idx| std::mem::take(&mut labels[file_idx]))
            .collect_vec();

        let mut confidence_idx = 0_usize;
        for (lineage, _) in &sorted_lineages {
            let levels = lineage.split(',').collect_vec();
            let last_level_idx = levels.len() - 1;
            let mut current_node = &mut root;
            for (level, label) in levels.into_iter().enumerate() {
                let node_type = if level == last_level_idx {
                    NodeType::Taxon
                } else {
                    NodeType::Inner
                };
                match &current_node.get_last_child_label() {
                    Some(name) => {
                        if name.as_str() != label {
                            current_node.add_child(Node::new(
                                label.to_string(),
                                confidence_idx,
                                node_type,
                            ));
                        }
                        current_node.confidence_range.1 = confidence_idx + 1;
                    }
                    None => {
                        current_node.add_child(Node::new(
                            label.to_string(),
                            confidence_idx,
                            node_type,
                        ));
                        current_node.confidence_range.1 = confidence_idx + 1;
                    }
                };
                if level == last_level_idx {
                    confidence_idx += 1;
                }
                current_node = current_node.children.last_mut().unwrap();
            }
            current_node.add_child(Node::new(
                current_node.label.clone(),
                confidence_idx - 1,
                NodeType::Sequence,
            ));
            current_node.confidence_range.1 = confidence_idx;
        }
        root.confidence_range.1 = confidence_idx;

        let mut bin_idx_to_lineage_idxs: Vec<Vec<usize>> = Vec::new();
        let bin_id_idx_pairs = sorted_lineages
            .iter()
            .enumerate()
            .flat_map(|(idx, (_, bin_id))| bin_id.clone().map(|b| (b, idx)))
            .sorted()
            .collect_vec();
        let mut lineage_idx_to_bin_idx: Vec<Option<usize>> = vec![None; confidence_idx];
        if let Some(pair) = bin_id_idx_pairs.first() {
            let mut current_bin = &pair.0;
            let mut current_idxs = Vec::new();
            let mut bin_idx = 0_usize;
            for (bin, idx) in bin_id_idx_pairs.iter() {
                if *bin != *current_bin {
                    bin_idx_to_lineage_idxs.push(current_idxs.clone());
                    current_bin = bin;
                    bin_idx += 1;
                    current_idxs.clear();
                }
                current_idxs.push(*idx);
                lineage_idx_to_bin_idx[*idx] = Some(bin_idx);
            }
            bin_idx_to_lineage_idxs.push(current_idxs);
        };
        let (bins, _): (Vec<String>, Vec<usize>) = bin_id_idx_pairs.into_iter().unzip();
        let (lineages, _): (Vec<String>, Vec<Option<String>>) = sorted_lineages.into_iter().unzip();

        Ok(Self {
            root,
            lineages,
            bins: bins.into_iter().unique().collect_vec(),
            sorted_idx,
            encoding_data,
            bin_idx_to_lineage_idxs,
            lineage_idx_to_bin_idx,
            num_tips: confidence_idx,
        })
    }

    pub fn print(&self) {
        self.root.print(0);
    }

    pub fn is_inner_taxon_node(&self, node: &Node) -> bool {
        node.node_type == NodeType::Inner
    }

    pub fn is_taxon_leaf(&self, node: &Node) -> bool {
        node.node_type == NodeType::Taxon
    }
}

/// Maps each minimally-encoded canonical k-mer to the (lineage-sorted) indices of the
/// reference sequences of one reference batch that contain it.
///
/// CSR layout: for k-mer `k`, the sequences (sequenceIDs) containing it are
/// data[offsets[k]..offsets[k + 1]] (see [`KMerMap::kmer_row`])
pub struct KMerMap {
    pub offsets: Vec<usize>,
    pub data: Vec<IndexType>,

    /// Scratch buffer for the unique k-mers of the sequence currently being processed by
    /// `rebuild`. Kept as a field (rather than a local) so its capacity is reused across
    /// `rebuild` calls instead of being regrown from empty.
    kmers: Vec<u32>,
}

impl KMerMap {
    pub fn new(encoding_data: &KMerEncodingData) -> Self {
        Self {
            offsets: vec![0_usize; encoding_data.n_unique_codes as usize + 1],
            data: Vec::new(),
            kmers: Vec::new(),
        }
    }

    /// Rebuilds the map for a batch of reference sequences, where `sequences[0]` is the
    /// `first_file_idx`-th sequence of the reference file. The buffers are reused across
    /// batches, a single offsets-sized vector serves for counting, the prefix sum, and the
    /// fill cursor to avoid any allocation proportional to `n_unique_codes` beyond it.
    #[time("debug", "KMerMap::{}")]
    pub fn rebuild(&mut self, sequences: &[LabeledSequence], first_file_idx: usize, tree: &Tree) {
        let encoding_data = &tree.encoding_data;
        let sorted_idx = &tree.sorted_idx[first_file_idx..first_file_idx + sequences.len()];
        let n = encoding_data.n_unique_codes as usize;
        let offsets = &mut self.offsets;
        offsets.fill(0);

        // pass 1: count unique (k-mer, sequence) occurrences per k-mer, shifted one
        // slot to the right (offsets[k + 1] holds the count for k-mer k)
        // the unique k-mers are recomputed per sequence instead of being kept for all
        // sequences at once, as that would need to coexist in memory with the k-mer map
        // one k-mer buffer is reused both passes across all ref sequences
        let kmers = &mut self.kmers;
        for LabeledSequence { sequence, .. } in sequences {
            seq_to_unique_minenc_canon_kmers_into(sequence, encoding_data, kmers);
            for &k_mer in kmers.iter() {
                offsets[k_mer as usize + 1] += 1;
            }
        }

        // prefix sum in place -> offsets[k] is now the CSR row start for k-mer k
        for i in 0..n {
            offsets[i + 1] += offsets[i];
        }

        self.data.clear();
        self.data.resize(offsets[n], IndexType::default());

        // pass 2: fill data, recomputing the unique k-mers of each sequence,
        // temporarily reusing offsets[0..n] as a per-k-mer write cursor
        // (destructive: offsets[k] ends up equal to its own former row-end,
        // i.e. shifted one slot from where it started)
        // sequences are visited in sorted-index order, so every row ends up ascending, which
        // keeps the scattered increments into the intersection buffers cache-friendly
        let mut order = (0..sequences.len()).collect_vec();
        order.sort_unstable_by_key(|&i| sorted_idx[i]);
        for i in order {
            let idx = sorted_idx[i];
            seq_to_unique_minenc_canon_kmers_into(&sequences[i].sequence, encoding_data, kmers);
            for &k_mer in kmers.iter() {
                let pos = &mut offsets[k_mer as usize];
                self.data[*pos] = idx;
                *pos += 1;
            }
        }

        // undo the shift from pass 2 to restore correct row-start offsets
        for i in (0..n).rev() {
            offsets[i + 1] = offsets[i];
        }
        offsets[0] = 0;

        if log_enabled!(Level::Debug) {
            // log the size of the k-mer map
            let offsets_size = self.offsets.capacity() * size_of::<usize>();
            let data_size = self.data.capacity() * size_of::<IndexType>();

            debug!(
                "size of k-mer map: {} of offsets, {} of data, {} total",
                HumanBytes(offsets_size as u64),
                HumanBytes(data_size as u64),
                HumanBytes((offsets_size + data_size) as u64)
            );
        }
    }

    pub fn kmer_row(&self, kmer: usize) -> &[IndexType] {
        &self.data[self.offsets[kmer]..self.offsets[kmer + 1]]
    }
}

#[derive(PartialEq, Eq, Debug)]
enum NodeType {
    Inner,
    Taxon,
    Sequence,
}

#[derive(Debug, Eq, PartialEq)]
pub struct Node {
    pub label: String,
    pub confidence_range: (usize, usize),
    pub children: Vec<Node>,
    node_type: NodeType,
}

impl Node {
    fn new(label: String, confidence_idx: usize, node_type: NodeType) -> Self {
        Self {
            label,
            confidence_range: (confidence_idx, confidence_idx + 1),
            children: vec![],
            node_type,
        }
    }

    fn add_child(&mut self, child: Node) {
        self.children.push(child);
    }

    fn get_last_child_label(&self) -> Option<&String> {
        match &self.children.last() {
            Some(c) => Some(&c.label),
            None => None,
        }
    }

    fn print(&self, depth: usize) {
        println!(
            "{}{} {:?}",
            "  ".repeat(depth),
            self.label,
            self.confidence_range
        );
        for child in &self.children {
            child.print(depth + 1);
        }
    }
}
