use anyhow::Result;
use std::{
    fs::File,
    io::{BufWriter, Read},
    path::PathBuf,
};

use indicatif::{HumanBytes, ProgressBar, ProgressIterator, ProgressStyle};
use itertools::Itertools;
use log::{debug, log_enabled, Level};
use logging_timer::time;
use serde::{Deserialize, Serialize};

use crate::{
    parser::LineageBinPair,
    utils::{seq_to_unique_minenc_canon_kmers, KMerEncodingData},
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

#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct Tree {
    pub root: Node,
    pub lineages: Vec<String>,
    pub bins: Vec<String>,

    // CSR layout: for k-mer `k`, the sequences (sequenceIDs) containing it are
    // k_mer_map_data[k_mer_map_offsets[k]..k_mer_map_offsets[k + 1]] (see `Tree::kmer_row`)
    pub k_mer_map_offsets: Vec<usize>,
    pub k_mer_map_data: Vec<IndexType>,

    pub encoding_data: KMerEncodingData,
    pub bin_idx_to_lineage_idxs: Vec<Vec<usize>>,
    pub lineage_idx_to_bin_idx: Vec<Option<usize>>,
    pub num_tips: usize,
}

impl Tree {
    #[time("debug", "Tree::{}")]
    pub fn new(
        labels: Vec<(String, Option<String>)>,
        sequences: Vec<Vec<u8>>,
        encoding_data: KMerEncodingData,
    ) -> Result<Self> {
        check_lineage_size(labels.len());
        let mut root = Node::new(String::from("root"), 0, NodeType::Inner);
        let mut lineage_sequence_pairs = labels.into_iter().zip_eq(sequences).collect_vec();

        lineage_sequence_pairs.sort_by(|(l1, _), (l2, _)| l1.cmp(l2));
        let mut confidence_idx = 0_usize;
        for ((lineage, _), _) in &lineage_sequence_pairs {
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
        let (sorted_lineages, sequences): (Vec<LineageBinPair>, Vec<Vec<u8>>) =
            lineage_sequence_pairs.into_iter().unzip();

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

        // Build the k-mer map as a CSR structure (offsets + flat data), reusing a
        // single offsets-sized vector for counting, the prefix sum, and the fill
        // cursor to avoid any allocation proportional to `n_unique_codes` beyond it.
        let n = encoding_data.n_unique_codes as usize;
        let mut k_mer_map_offsets = vec![0_usize; n + 1];

        // pass 1: count unique (k-mer, sequence) occurrences per k-mer, shifted one
        // slot to the right (k_mer_map_offsets[k + 1] holds the count for k-mer k)
        // the unique k-mers are recomputed per sequence instead of being kept for all
        // sequences at once, as that would need to coexist in memory with the k_mer_map
        for sequence in &sequences {
            for k_mer in seq_to_unique_minenc_canon_kmers(sequence, &encoding_data) {
                k_mer_map_offsets[k_mer as usize + 1] += 1;
            }
        }

        // prefix sum in place -> k_mer_map_offsets[k] is now the CSR row start for k-mer k
        for i in 0..n {
            k_mer_map_offsets[i + 1] += k_mer_map_offsets[i];
        }

        let mut k_mer_map_data = vec![IndexType::default(); k_mer_map_offsets[n]];

        let pb = if cfg!(test) {
            ProgressBar::hidden()
        } else {
            ProgressBar::new(sequences.len() as u64).with_style(
                ProgressStyle::with_template(
                    "[{elapsed_precise}] {bar:80.cyan/blue} {pos:>7}/{len:7}[ETA:{eta}] {msg}",
                )
                .unwrap()
                .progress_chars("##-"),
            )
        };

        // pass 2: fill data, recomputing the unique k-mers of each sequence,
        // temporarily reusing k_mer_map_offsets[0..n] as a
        // per-k-mer write cursor (destructive: k_mer_map_offsets[k] ends up equal
        // to its own former row-end, i.e. shifted one slot from where it started)
        for (idx, sequence) in sequences
            .iter()
            .enumerate()
            .progress_with(pb)
            .with_message("Creating k-mer map...")
        {
            for k_mer in seq_to_unique_minenc_canon_kmers(sequence, &encoding_data) {
                let pos = &mut k_mer_map_offsets[k_mer as usize];
                k_mer_map_data[*pos] = idx as IndexType;
                *pos += 1;
            }
        }

        // undo the shift from pass 2 to restore correct row-start offsets
        for i in (0..n).rev() {
            k_mer_map_offsets[i + 1] = k_mer_map_offsets[i];
        }
        k_mer_map_offsets[0] = 0;

        if log_enabled!(Level::Debug) {
            // log the size of the k_mer_map
            let offsets_size = k_mer_map_offsets.capacity() * size_of::<usize>();
            let data_size = k_mer_map_data.capacity() * size_of::<IndexType>();

            debug!(
                "size of k_mer_map: {} of offsets, {} of data, {} total",
                HumanBytes(offsets_size as u64),
                HumanBytes(data_size as u64),
                HumanBytes((offsets_size + data_size) as u64)
            );
        }

        Ok(Self {
            root,
            lineages,
            bins: bins.into_iter().unique().collect_vec(),
            k_mer_map_offsets,
            k_mer_map_data,
            encoding_data,
            bin_idx_to_lineage_idxs,
            lineage_idx_to_bin_idx,
            num_tips: confidence_idx,
        })
    }

    pub fn print(&self) {
        self.root.print(0);
    }

    pub fn kmer_row(&self, kmer: usize) -> &[IndexType] {
        &self.k_mer_map_data[self.k_mer_map_offsets[kmer]..self.k_mer_map_offsets[kmer + 1]]
    }

    #[time("debug")]
    pub fn save_to_file(&self, mut output: BufWriter<File>) -> Result<()> {
        if log_enabled!(Level::Info) {
            eprintln!("[INFO ] Writing database to file...");
        }
        bincode::serialize_into(&mut output, &self)?;
        Ok(())
    }

    pub fn load_from_file(path: &PathBuf) -> Result<Self> {
        if log_enabled!(Level::Info) {
            eprintln!("[INFO ] Trying to read from database file...");
        }
        let mut file = File::open(path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;
        let decoded: Self = bincode::deserialize(&buffer)?;
        Ok(decoded)
    }

    pub fn is_inner_taxon_node(&self, node: &Node) -> bool {
        node.node_type == NodeType::Inner
    }

    pub fn is_taxon_leaf(&self, node: &Node) -> bool {
        node.node_type == NodeType::Taxon
    }
}

#[derive(PartialEq, Eq, Debug, Serialize, Deserialize)]
enum NodeType {
    Inner,
    Taxon,
    Sequence,
}

#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
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
