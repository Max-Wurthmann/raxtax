use anyhow::{bail, Result};
use bincode::Options;
use std::{
    fs::File,
    io::{BufReader, BufWriter},
    path::PathBuf,
    time::Duration,
};

use indicatif::{HumanBytes, ProgressBar, ProgressStyle};
use itertools::Itertools;
use log::{debug, log_enabled, Level};
use logging_timer::{time, timer};
use serde::{Deserialize, Serialize};

use crate::{
    parser::{parse_reference_label, reference_label_regex, LineageBinPair, SequenceReader},
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
    /// Builds the reference database by streaming the reference file twice: once to parse the
    /// labels and count how many sequences contain each k-mer, and once to fill the CSR k-mer map.
    /// The reference sequences are therefore never all in memory at once, which for a database of
    /// whole genomes is the difference between a peak proportional to the uncompressed database and
    /// one proportional to the k-mer map. `open_records` is called once per pass and has to yield
    /// the same records both times.
    #[time("debug", "Tree::{}")]
    pub fn build(
        open_records: impl Fn() -> Result<SequenceReader>,
        encoding_data: KMerEncodingData,
    ) -> Result<Self> {
        // Build the k-mer map as a CSR structure (offsets + flat data), reusing a
        // single offsets-sized vector for counting, the prefix sum, and the fill
        // cursor to avoid any allocation proportional to `n_unique_codes` beyond it.
        let n = encoding_data.n_unique_codes as usize;
        let mut k_mer_map_offsets = vec![0_usize; n + 1];
        let mut labels: Vec<LineageBinPair> = Vec::new();
        // The unique k-mers of a sequence are recomputed in the second pass instead of being kept
        // for all sequences at once, as that would need to coexist in memory with the k_mer_map.
        // One buffer is reused across sequences and across both passes.
        let mut k_mers: Vec<u32> = Vec::new();

        // pass 1: parse the labels and count unique (k-mer, sequence) occurrences per k-mer,
        // shifted one slot to the right (k_mer_map_offsets[k + 1] holds the count for k-mer k)
        {
            let _tmr = timer!(Level::Info; "Read labels and count k-mers");
            let pb = if cfg!(test) {
                ProgressBar::hidden()
            } else {
                ProgressBar::new_spinner().with_message("Parsing Reference...")
            };
            pb.enable_steady_tick(Duration::from_millis(100));
            let label_regex = reference_label_regex()?;
            for record in open_records()? {
                let (label, sequence) = record?;
                labels.push(parse_reference_label(&label_regex, &label)?);
                seq_to_unique_minenc_canon_kmers_into(&sequence, &encoding_data, &mut k_mers);
                for k_mer in &k_mers {
                    k_mer_map_offsets[*k_mer as usize + 1] += 1;
                }
                pb.inc(1);
            }
            pb.finish();
        }
        if labels.is_empty() {
            bail!("File is empty")
        }
        check_lineage_size(labels.len());
        let n_references = labels.len();

        // References are indexed in lineage-sorted order, because every `Node::confidence_range`
        // has to cover a contiguous range of them. Sorting an index vector gives that order
        // without holding on to the sequences: `rank` then maps a sequence's position in the file
        // to its index in the k-mer map. Labels are compared as whole (lineage, bin) pairs, so
        // sequences sharing a lineage are ordered by their bin.
        let mut order: Vec<usize> = (0..n_references).collect();
        order.sort_by(|&a, &b| labels[a].cmp(&labels[b]));
        let mut rank = vec![IndexType::default(); n_references];
        for (sorted_idx, &file_idx) in order.iter().enumerate() {
            rank[file_idx] = sorted_idx as IndexType;
        }
        let sorted_labels: Vec<LineageBinPair> = order
            .into_iter()
            .map(|file_idx| std::mem::take(&mut labels[file_idx]))
            .collect();

        let mut root = Node::new(String::from("root"), 0, NodeType::Inner);
        let mut confidence_idx = 0_usize;
        for (lineage, _) in &sorted_labels {
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
        let (lineages, bin_ids): (Vec<String>, Vec<Option<String>>) =
            sorted_labels.into_iter().unzip();

        let mut bin_idx_to_lineage_idxs: Vec<Vec<usize>> = Vec::new();
        let bin_id_idx_pairs = bin_ids
            .into_iter()
            .enumerate()
            .flat_map(|(idx, bin_id)| bin_id.map(|b| (b, idx)))
            .sorted()
            .collect_vec();
        let mut lineage_idx_to_bin_idx: Vec<Option<usize>> = vec![None; confidence_idx];
        if let Some(pair) = bin_id_idx_pairs.first() {
            let mut current_bin = &pair.0;
            let mut current_idxs = Vec::new();
            let mut bin_idx = 0_usize;
            for (bin, idx) in bin_id_idx_pairs.iter() {
                if *bin != *current_bin {
                    bin_idx_to_lineage_idxs.push(std::mem::take(&mut current_idxs));
                    current_bin = bin;
                    bin_idx += 1;
                }
                current_idxs.push(*idx);
                lineage_idx_to_bin_idx[*idx] = Some(bin_idx);
            }
            bin_idx_to_lineage_idxs.push(current_idxs);
        };
        // the pairs are sorted by bin, so deduplicating adjacent entries leaves the distinct bins
        // in the same order `bin_idx` was handed out above
        let (mut bins, _): (Vec<String>, Vec<usize>) = bin_id_idx_pairs.into_iter().unzip();
        bins.dedup();

        // prefix sum in place -> k_mer_map_offsets[k] is now the CSR row start for k-mer k
        for i in 0..n {
            k_mer_map_offsets[i + 1] += k_mer_map_offsets[i];
        }

        let mut k_mer_map_data = vec![IndexType::default(); k_mer_map_offsets[n]];

        let pb = if cfg!(test) {
            ProgressBar::hidden()
        } else {
            ProgressBar::new(n_references as u64)
                .with_style(
                    ProgressStyle::with_template(
                        "[{elapsed_precise}] {bar:80.cyan/blue} {pos:>7}/{len:7}[ETA:{eta}] {msg}",
                    )
                    .unwrap()
                    .progress_chars("##-"),
                )
                .with_message("Creating k-mer map...")
        };

        // pass 2: stream the file again to fill the data, temporarily reusing
        // k_mer_map_offsets[0..n] as a per-k-mer write cursor (destructive: k_mer_map_offsets[k]
        // ends up equal to its own former row-end, i.e. shifted one slot from where it started)
        let mut n_records = 0_usize;
        for record in open_records()? {
            let (_, sequence) = record?;
            let Some(&idx) = rank.get(n_records) else {
                bail!("Reference file gained sequences while it was being read");
            };
            seq_to_unique_minenc_canon_kmers_into(&sequence, &encoding_data, &mut k_mers);
            for k_mer in &k_mers {
                let pos = &mut k_mer_map_offsets[*k_mer as usize];
                k_mer_map_data[*pos] = idx;
                *pos += 1;
            }
            n_records += 1;
            pb.inc(1);
        }
        pb.finish();
        if n_records != n_references {
            bail!("Reference file changed while it was being read: {n_references} sequences in the first pass, {n_records} in the second");
        }

        // undo the shift from pass 2 to restore correct row-start offsets
        for i in (0..n).rev() {
            k_mer_map_offsets[i + 1] = k_mer_map_offsets[i];
        }
        k_mer_map_offsets[0] = 0;

        // pass 2 writes the sequence indices in file order rather than in index order, so sort
        // every row to restore the ascending order `Tree::kmer_row`'s consumers scan
        for row in k_mer_map_offsets.windows(2) {
            k_mer_map_data[row[0]..row[1]].sort_unstable();
        }

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
            bins,
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
        let file = File::open(path)?;
        // Nothing in a valid database is longer than the file holding it, so bound the read by the
        // file size. Without a limit a truncated or non-bincode file is read as a garbage length
        // prefix and aborts the process on a huge allocation instead of reporting an error.
        // The rest of the configuration is what `bincode::deserialize_from` uses, so the on-disk
        // format is unchanged.
        let limit = file.metadata()?.len();
        Ok(bincode::DefaultOptions::new()
            .with_fixint_encoding()
            .allow_trailing_bytes()
            .with_limit(limit)
            .deserialize_from(BufReader::new(file))?)
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
