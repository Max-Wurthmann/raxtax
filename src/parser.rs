use ahash::HashSet;
use anyhow::{bail, Context, Result};
use indicatif::{ProgressIterator, ProgressStyle};
use log::{log_enabled, warn, Level};
use logging_timer::{time, timer};
use regex::Regex;
use std::{
    io::{BufRead, BufReader, Read},
    path::PathBuf,
};

use crate::{
    tree::Tree,
    utils::{self, KMerEncodingData},
};

pub type LineageBinPair = (String, Option<String>);

fn map_dna_char(ch: char) -> u8 {
    let a: u8 = 0b0001;
    let c: u8 = 0b0010;
    let g: u8 = 0b0100;
    let t: u8 = 0b1000;
    match ch.to_ascii_uppercase() {
        'A' => a,
        'C' => c,
        'G' => g,
        'T' => t,
        'W' => a | t,
        'S' => c | g,
        'M' => a | c,
        'K' => g | t,
        'R' => a | g,
        'Y' => c | t,
        'B' => c | g | t,
        'D' => a | g | t,
        'H' => a | c | t,
        'V' => a | c | g,
        'N' => a | c | g | t,
        _ => panic!("Unexpected character: {ch}"),
    }
}

#[time("info", "Parsing References")]
pub fn parse_reference_fasta_file(
    sequence_path: &PathBuf,
    encoding_data: KMerEncodingData,
) -> Result<(bool, Tree)> {
    if let Ok(tree) = Tree::load_from_file(sequence_path) {
        if tree.encoding_data != encoding_data {
            let notification = format!("k-mer size of loaded tree (k={}) does not match specified k-mer size (k={}). Attempting to reparse reference file using k={} ...", tree.encoding_data.k, encoding_data.k, encoding_data.k);
            warn!("{}", notification);
            if log_enabled!(log::Level::Warn) {
                eprintln!("\x1b[33m[WARN ]\x1b[0m {}", notification);
            }
        } else {
            return Ok((false, tree));
        }
    }
    let mut fasta_str = String::new();
    let _ = utils::get_reader(sequence_path)?.read_to_string(&mut fasta_str);
    Ok((true, parse_reference_fasta_str(&fasta_str, encoding_data)?))
}

fn parse_reference_fasta_str(fasta_str: &str, encoding_data: KMerEncodingData) -> Result<Tree> {
    if fasta_str.is_empty() {
        bail!("File is empty")
    }
    let regex = Regex::new(r"tax=([^;]+);([^;]+)*")?;
    let (labels, sequences) = {
        let _tmr = timer!(Level::Info; "Read file and create k-mer mapping");
        let lines: Vec<String> = fasta_str
            .lines()
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty() && !l.starts_with(';'))
            .collect();
        if !lines[0].starts_with('>') {
            bail!("Not a valid FASTA file")
        }
        let mut labels: Vec<LineageBinPair> = Vec::new();
        let mut sequences: Vec<Vec<u8>> = Vec::new();
        let mut current_sequence = Vec::<u8>::new();
        // let mut bin_id_to_lineages: HashMap<String, Vec<String>> = HashMap::new();

        // create label and sequence vectors
        lines
            .into_iter()
            .progress_with_style(
                ProgressStyle::with_template(
                    "[{elapsed_precise}] {bar:80.cyan/blue} {pos:>7}/{len:7}[ETA:{eta}] {msg}",
                )
                .unwrap()
                .progress_chars("##-"),
            )
            .with_message("Parsing Reference...")
            .map(|line| -> Result<()> {
                if let Some(label) = line.strip_prefix('>') {
                    let caps = regex.captures(label).context(format!(
                        "Unexpected taxonomical annotation detected in label {label}"
                    ))?;
                    let lineage = caps
                        .get(1)
                        .context(format!("No taxonomic string found in label {label}"))?
                        .as_str()
                        .to_owned();
                    let bin = caps.get(2).map(|bin| bin.as_str().to_owned());
                    labels.push((lineage, bin));
                    if !current_sequence.is_empty() {
                        sequences.push(current_sequence.clone());
                        current_sequence = Vec::new();
                    }
                } else {
                    current_sequence.extend(line.chars().map(|c| -> u8 { map_dna_char(c) }));
                }
                Ok(())
            })
            .collect::<Result<Vec<()>>>()?;
        sequences.push(current_sequence);
        if labels.len() != sequences.len() {
            bail!("Number of sequences does not match number of labels")
        }
        (labels, sequences)
    };
    Tree::new(labels, sequences, encoding_data)
}

/// Streams a query FASTA file record by record instead of loading it into
/// memory all at once, handing out fixed-size batches of parsed queries.
pub struct QueryBatchReader<'a> {
    reader: Box<dyn BufRead>,
    pending_header: Option<String>,
    queries_to_skip: &'a HashSet<String>,
}

impl<'a> QueryBatchReader<'a> {
    fn new<R: BufRead + 'static>(
        mut reader: R,
        queries_to_skip: &'a HashSet<String>,
    ) -> Result<Self> {
        let mut line = String::new();
        let pending_header = loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                bail!("File is empty")
            }
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with(';') {
                continue;
            }
            match trimmed.strip_prefix('>') {
                Some(label) => break label.to_string(),
                None => bail!("Not a valid FASTA file"),
            }
        };
        Ok(QueryBatchReader {
            reader: Box::new(reader),
            pending_header: Some(pending_header),
            queries_to_skip,
        })
    }

    /// Reads and returns up to `batch_size` queries. Returns an empty vector
    /// once the file has been fully consumed.
    pub fn next_batch(&mut self, batch_size: usize) -> Result<Vec<(String, Vec<u8>)>> {
        let mut batch = Vec::with_capacity(batch_size);
        let mut line = String::new();
        while batch.len() < batch_size {
            let Some(label) = self.pending_header.take() else {
                break;
            };
            let mut sequence = Vec::new();
            loop {
                line.clear();
                if self.reader.read_line(&mut line)? == 0 {
                    break;
                }
                let trimmed = line.trim();
                if trimmed.is_empty() || trimmed.starts_with(';') {
                    continue;
                }
                if let Some(next_label) = trimmed.strip_prefix('>') {
                    self.pending_header = Some(next_label.to_string());
                    break;
                }
                sequence.extend(trimmed.chars().map(|c| -> u8 { map_dna_char(c) }));
            }
            if !self.queries_to_skip.contains(&label) {
                batch.push((label, sequence));
            }
        }
        Ok(batch)
    }
}

#[time("info", "Opening Query File")]
pub fn open_query_batch_reader<'a>(
    sequence_path: &PathBuf,
    queries_to_skip: &'a HashSet<String>,
) -> Result<QueryBatchReader<'a>> {
    let reader = BufReader::new(utils::get_reader(sequence_path)?);
    QueryBatchReader::new(reader, queries_to_skip)
}

#[cfg(test)]
mod tests {
    use ahash::{HashSet, HashSetExt};

    use itertools::Itertools;

    use crate::{
        tree::Tree,
        utils::{encode, reverse_complement, KMerEncodingData},
    };

    use std::io::Cursor;

    use super::{parse_reference_fasta_str, QueryBatchReader};

    #[test]
    fn test_str_parser() {
        let minenc_8mer = |kmer| {
            encode(
                kmer,
                reverse_complement(kmer, 8),
                &KMerEncodingData::new(8).unwrap(),
            ) as usize
        };
        let fasta_str = r">Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species1;
AAACCCTTTGGGA
>Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species2;
ATACGCTTTGGGA
>Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order4,f:Family5,g:Genus2,s:Species3;
ATCCGCTATGGGA
>Badabing|Badabum;tax=p:Phylum1,c:Class2,o:Order2,f:Family3,g:Genus3,s:Species6;
ATACGCTTTGCGT
>Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species2;
GTGCGCTATGCGA
>Badabing|Badabum;tax=p:Phylum2,c:Class3,o:Order3,f:Family4,g:Genus4,s:Species5;
ATACGCTTTGCGT";

        let tree = parse_reference_fasta_str(fasta_str, KMerEncodingData::new(8).unwrap()).unwrap();
        for (k, v) in tree.k_mer_map.iter().enumerate() {
            if !v.is_empty() {
                println!("{k:b}:\n {v:?}");
            }
        }
        assert_eq!(
            tree.k_mer_map[minenc_8mer(0b0001_0101_1111_1110)]
                .iter()
                .collect_vec(),
            &[&0]
        );
        assert_eq!(
            tree.k_mer_map[minenc_8mer(0b0000_1001_1011_0011)]
                .iter()
                .sorted()
                .collect_vec(),
            &[&1, &4, &5]
        );
        assert_eq!(
            tree.k_mer_map[minenc_8mer(0b0101_0011_0010_0110)]
                .iter()
                .collect_vec(),
            &[&3]
        );
        assert_eq!(tree.num_tips, 6);
        assert_eq!(
            tree.lineages,
            vec![
                String::from("p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species1"),
                "p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species2".into(),
                "p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species2".into(),
                "p:Phylum1,c:Class1,o:Order4,f:Family5,g:Genus2,s:Species3".into(),
                "p:Phylum1,c:Class2,o:Order2,f:Family3,g:Genus3,s:Species6".into(),
                "p:Phylum2,c:Class3,o:Order3,f:Family4,g:Genus4,s:Species5".into(),
            ]
        );
    }

    #[test]
    fn test_query_parser() {
        let skip = HashSet::new();
        let fasta_str = r">label1
AAACCCTTTGGGA";
        let mut reader = QueryBatchReader::new(Cursor::new(fasta_str), &skip).unwrap();
        let (_, sequence) = &reader.next_batch(10).unwrap()[0];
        assert_eq!(sequence, &[1, 1, 1, 2, 2, 2, 8, 8, 8, 4, 4, 4, 1]);

        let fasta_str2 = r">label1
ACGTWSMKRYBDHVN";
        let mut reader2 = QueryBatchReader::new(Cursor::new(fasta_str2), &skip).unwrap();
        let (_, sequence) = &reader2.next_batch(10).unwrap()[0];
        assert_eq!(
            sequence,
            &[1, 2, 4, 8, 9, 6, 3, 12, 5, 10, 14, 13, 11, 7, 15]
        );
    }

    #[test]
    fn test_query_batch_reader_batches_and_skips() {
        let fasta_str = ">q1\nAAAA\n>q2\nCCCC\n>q3\nGGGG\n>q4\nTTTT\n";
        let mut skip = HashSet::new();
        skip.insert("q2".to_string());
        let mut reader = QueryBatchReader::new(Cursor::new(fasta_str), &skip).unwrap();

        let batch1 = reader.next_batch(2).unwrap();
        assert_eq!(
            batch1.iter().map(|(l, _)| l.clone()).collect_vec(),
            vec!["q1".to_string(), "q3".to_string()]
        );

        let batch2 = reader.next_batch(2).unwrap();
        assert_eq!(
            batch2.iter().map(|(l, _)| l.clone()).collect_vec(),
            vec!["q4".to_string()]
        );

        let batch3 = reader.next_batch(2).unwrap();
        assert!(batch3.is_empty());
    }

    #[test]
    fn test_kmers() {
        let minenc_8mer = |kmer| {
            encode(
                kmer,
                reverse_complement(kmer, 8),
                &KMerEncodingData::new(8).unwrap(),
            ) as usize
        };
        let fasta_str = r">Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species1;
AAACCCCGT
>Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus1,s:Species1;
TAACCCCGG
>Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus2,s:Species3;
TTTAAAACC
>Badabing|Badabum;tax=p:Phylum1,c:Class1,o:Order1,f:Family1,g:Genus2,s:Species3;
TTTAAAACA
>Badabing|Badabum;tax=p:Phylum1,c:Class2,o:Order2,f:Family2,g:Genus3,s:Species4;
AAACCCCGG";
        /* reverse complements:
         * 0: AAACCCCGT -> ACGGGGTTT
         * 1: TAACCCCGG -> CCGGGGTTA
         * 2: TTTAAAACC -> GGTTTTAAA
         * 3: TTTAAAACA -> TGTTTTAAA
         * 4: AAACCCCGG -> CCGGGGTTT
         */

        let Tree { k_mer_map, .. } =
            parse_reference_fasta_str(fasta_str, KMerEncodingData::new(8).unwrap()).unwrap();
        for (k, v) in k_mer_map.iter().enumerate() {
            if !v.is_empty() {
                println!("{k:b}:\n {v:?}");
            }
        }

        assert_eq!(
            // AAACCCCG
            k_mer_map[minenc_8mer(0b01_0101_0110)]
                .iter()
                .sorted()
                .collect_vec(),
            &[&0, &4]
        );
        assert_eq!(
            // AACCCCGG
            k_mer_map[minenc_8mer(0b0101_0101_1010)]
                .iter()
                .sorted()
                .collect_vec(),
            &[&1, &4]
        );
        assert_eq!(
            // AACCCCGT
            k_mer_map[minenc_8mer(0b101_0101_1011)]
                .iter()
                .sorted()
                .collect_vec(),
            &[&0]
        );
        assert_eq!(
            // CGGGGTTA
            k_mer_map[minenc_8mer(0b0110_1010_1011_1100)]
                .iter()
                .sorted()
                .collect_vec(),
            &[&1]
        );
        assert_eq!(
            // GGTTTTAA
            k_mer_map[minenc_8mer(0b1010_1111_1111_0000)]
                .iter()
                .sorted()
                .collect_vec(),
            &[&2]
        );
        // GGTTTTAA
        assert_eq!(
            k_mer_map[minenc_8mer(0b1010_1111_1111_0000)]
                .iter()
                .sorted()
                .collect_vec(),
            &[&2]
        );
        // GTTTTAAA
        assert_eq!(
            k_mer_map[minenc_8mer(0b1011_1111_1100_0000)]
                .iter()
                .sorted()
                .collect_vec(),
            &[&2, &3]
        );
    }
}
