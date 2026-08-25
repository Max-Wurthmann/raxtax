use ahash::HashSet;
use anyhow::{anyhow, bail, Context, Result};
use indicatif::{ProgressIterator, ProgressStyle};
use log::{log_enabled, warn, Level};
use logging_timer::{time, timer};
use regex::Regex;
use std::{
    io::{BufRead, Read},
    path::PathBuf,
};

use crate::{
    tree::Tree,
    utils::{self, KMerEncodingData},
};

pub type LineageBinPair = (String, Option<String>);
pub type LabeledSequence = (String, Vec<u8>);

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

/// Streams labeled sequences out of a FASTA file, record by record.
pub struct FastaSequenceReader {
    reader: Box<dyn BufRead>,
    next_header: Option<String>,
    line: String,
}

impl FastaSequenceReader {
    pub fn new(reader: Box<dyn BufRead>) -> Self {
        FastaSequenceReader {
            reader,
            next_header: None,
            line: String::new(),
        }
    }

    /// Advances `self.line` to the next non-empty, non-comment line.
    /// Returns `Ok(Some(()))` when such a line was found, `Ok(None)` on EOF.
    fn read_next_line(&mut self) -> Result<Option<()>> {
        loop {
            self.line.clear();
            if self.reader.read_line(&mut self.line)? == 0 {
                return Ok(None);
            }
            let trimmed = self.line.trim();
            if trimmed.is_empty() || trimmed.starts_with(';') {
                continue;
            }
            return Ok(Some(()));
        }
    }
}

impl Iterator for FastaSequenceReader {
    type Item = Result<LabeledSequence>;

    fn next(&mut self) -> Option<Self::Item> {
        let label = match self.next_header.take() {
            Some(label) => label,
            None => match self.read_next_line() {
                Ok(Some(())) => match self.line.trim().strip_prefix('>') {
                    Some(label) => label.to_string(),
                    None => return Some(Err(anyhow!("Not a valid FASTA file"))),
                },
                Ok(None) => return None,
                Err(e) => return Some(Err(e)),
            },
        };
        let mut sequence = Vec::new();
        loop {
            match self.read_next_line() {
                Ok(Some(())) => {}
                Ok(None) => break,
                Err(e) => return Some(Err(e)),
            }
            let trimmed = self.line.trim();
            if let Some(next_label) = trimmed.strip_prefix('>') {
                self.next_header = Some(next_label.to_string());
                break;
            }
            sequence.extend(trimmed.chars().map(|c| -> u8 { map_dna_char(c) }));
        }
        Some(Ok((label, sequence)))
    }
}

/// Streams labeled sequences out of a FASTQ file, record by record.
/// Each record is expected to span exactly four lines (header, sequence,
/// `+` separator, quality), which covers the near-universal single-line
/// FASTQ convention produced by sequencers and downstream tools.
pub struct FastqSequenceReader {
    reader: Box<dyn BufRead>,
    header: String,
    sequence_line: String,
    plus_line: String,
    quality_line: String,
}

impl FastqSequenceReader {
    pub fn new(reader: Box<dyn BufRead>) -> Self {
        FastqSequenceReader {
            reader,
            header: String::new(),
            sequence_line: String::new(),
            plus_line: String::new(),
            quality_line: String::new(),
        }
    }

    fn read_line(reader: &mut dyn BufRead, line: &mut String) -> Result<Option<()>> {
        line.clear();
        match reader.read_line(line) {
            Err(e) => Err(e.into()),
            Ok(0) => Ok(None),
            Ok(_) => Ok(Some(())),
        }
    }
}

impl Iterator for FastqSequenceReader {
    type Item = Result<LabeledSequence>;

    fn next(&mut self) -> Option<Self::Item> {
        match Self::read_line(&mut *self.reader, &mut self.header) {
            Err(e) => return Some(Err(e)),
            Ok(None) => return None,
            Ok(Some(())) => {}
        }

        for line in [
            &mut self.sequence_line,
            &mut self.plus_line,
            &mut self.quality_line,
        ] {
            match Self::read_line(&mut *self.reader, line) {
                Err(e) => return Some(Err(e)),
                Ok(None) => {
                    return Some(Err(anyhow!(
                        "Unexpected EOF reached while reading FASTQ file"
                    )))
                }
                Ok(Some(())) => {}
            }
        }

        let Some(label) = self.header.trim().strip_prefix('@') else {
            return Some(Err(anyhow!("Not a valid FASTQ file")));
        };

        let Some(_) = self.plus_line.trim().strip_prefix('+') else {
            return Some(Err(anyhow!("Not a valid FASTQ file")));
        };

        let sequence = self
            .sequence_line
            .trim()
            .chars()
            .map(|c| -> u8 { map_dna_char(c) })
            .collect();

        Some(Ok((label.to_string(), sequence)))
    }
}

/// Wraps any iterator of labeled sequences (as produced by
/// [`FastaSequenceReader`] or [`FastqSequenceReader`]) and hands out
/// fixed-size batches instead of loading the whole file into memory.
pub struct BatchedSequenceReader<'a, I: Iterator<Item = Result<LabeledSequence>>> {
    inner: I,
    seq_to_skip: &'a HashSet<String>,
}

impl<'a, I: Iterator<Item = Result<LabeledSequence>>> BatchedSequenceReader<'a, I> {
    pub fn new(inner: I, seq_to_skip: &'a HashSet<String>) -> Self {
        BatchedSequenceReader { inner, seq_to_skip }
    }

    /// Reads and returns up to `batch_size` queries. Returns `Ok(None)` once
    /// the underlying iterator has been fully consumed.
    pub fn next_batch(&mut self, batch_size: usize) -> Result<Option<Vec<LabeledSequence>>> {
        let mut batch = Vec::with_capacity(batch_size);
        while batch.len() < batch_size {
            match self.inner.next() {
                Some(Ok((label, sequence))) => {
                    if !self.seq_to_skip.contains(&label) {
                        batch.push((label, sequence));
                    }
                }
                Some(Err(e)) => return Err(e),
                None => break,
            }
        }
        if batch.is_empty() {
            Ok(None)
        } else {
            Ok(Some(batch))
        }
    }
}

impl<'a> BatchedSequenceReader<'a, Box<dyn Iterator<Item = Result<LabeledSequence>>>> {
    /// Opens `path` (transparently decompressing `.gz`/`.gzip` files) and
    /// picks a FASTA or FASTQ record reader based on the file extension
    /// (`.fastq`/`.fq`, or `.fasta`/`.fa`/anything else defaulting to FASTA).
    pub fn from_file(path: &PathBuf, seq_to_skip: &'a HashSet<String>) -> Result<Self> {
        let reader = utils::get_reader(path)?;
        let inner: Box<dyn Iterator<Item = Result<LabeledSequence>>> = if is_fastq_file(path) {
            Box::new(FastqSequenceReader::new(reader))
        } else {
            Box::new(FastaSequenceReader::new(reader))
        };
        Ok(BatchedSequenceReader::new(inner, seq_to_skip))
    }
}

fn extension_lowercase(path: &std::path::Path) -> String {
    path.extension()
        .and_then(|ext| ext.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase()
}

#[allow(clippy::ptr_arg)]
fn is_fastq_file(path: &PathBuf) -> bool {
    let ext = extension_lowercase(path);
    let ext = if ext == "gz" || ext == "gzip" {
        path.file_stem()
            .map(PathBuf::from)
            .map(|stem| extension_lowercase(&stem))
            .unwrap_or_default()
    } else {
        ext
    };
    matches!(ext.as_str(), "fastq" | "fq")
}

#[cfg(test)]
mod tests {
    use ahash::{HashSet, HashSetExt};

    use itertools::Itertools;

    use crate::{
        tree::Tree,
        utils::{encode, reverse_complement, KMerEncodingData},
    };

    use std::io::{Cursor, Write};

    use super::{
        parse_reference_fasta_str, BatchedSequenceReader, FastaSequenceReader, FastqSequenceReader,
    };

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
        let mut reader = BatchedSequenceReader::new(
            FastaSequenceReader::new(Box::new(Cursor::new(fasta_str))),
            &skip,
        );
        let batch = reader.next_batch(10).unwrap().unwrap();
        let (_, sequence) = &batch[0];
        assert_eq!(sequence, &[1, 1, 1, 2, 2, 2, 8, 8, 8, 4, 4, 4, 1]);

        let fasta_str2 = r">label1
ACGTWSMKRYBDHVN";
        let mut reader2 = BatchedSequenceReader::new(
            FastaSequenceReader::new(Box::new(Cursor::new(fasta_str2))),
            &skip,
        );
        let batch2 = reader2.next_batch(10).unwrap().unwrap();
        let (_, sequence) = &batch2[0];
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
        let mut reader = BatchedSequenceReader::new(
            FastaSequenceReader::new(Box::new(Cursor::new(fasta_str))),
            &skip,
        );

        let batch1 = reader.next_batch(2).unwrap().unwrap();
        assert_eq!(
            batch1.iter().map(|(l, _)| l.clone()).collect_vec(),
            vec!["q1".to_string(), "q3".to_string()]
        );

        let batch2 = reader.next_batch(2).unwrap().unwrap();
        assert_eq!(
            batch2.iter().map(|(l, _)| l.clone()).collect_vec(),
            vec!["q4".to_string()]
        );

        let batch3 = reader.next_batch(2).unwrap();
        assert!(batch3.is_none());
    }

    #[test]
    fn test_fasta_sequence_reader_iterates_records() {
        let fasta_str = ">q1\nAAAA\n>q2\nCCCC\n";
        let mut reader = FastaSequenceReader::new(Box::new(Cursor::new(fasta_str)));
        let (label1, seq1) = reader.next().unwrap().unwrap();
        assert_eq!(label1, "q1");
        assert_eq!(seq1, &[1, 1, 1, 1]);
        let (label2, seq2) = reader.next().unwrap().unwrap();
        assert_eq!(label2, "q2");
        assert_eq!(seq2, &[2, 2, 2, 2]);
        assert!(reader.next().is_none());
    }

    #[test]
    fn test_fastq_sequence_reader_iterates_records() {
        let fastq_str = "@q1\nAAAA\n+\nIIII\n@q2\nCCCC\n+q2\nIIII\n";
        let mut reader = FastqSequenceReader::new(Box::new(Cursor::new(fastq_str)));
        let (label1, seq1) = reader.next().unwrap().unwrap();
        assert_eq!(label1, "q1");
        assert_eq!(seq1, &[1, 1, 1, 1]);
        let (label2, seq2) = reader.next().unwrap().unwrap();
        assert_eq!(label2, "q2");
        assert_eq!(seq2, &[2, 2, 2, 2]);
        assert!(reader.next().is_none());
    }

    #[test]
    fn test_batched_sequence_reader_over_fastq() {
        let fastq_str = "@q1\nAAAA\n+\nIIII\n@q2\nCCCC\n+\nIIII\n";
        let skip = HashSet::new();
        let mut reader = BatchedSequenceReader::new(
            FastqSequenceReader::new(Box::new(Cursor::new(fastq_str))),
            &skip,
        );
        let batch = reader.next_batch(10).unwrap().unwrap();
        assert_eq!(
            batch.iter().map(|(l, _)| l.clone()).collect_vec(),
            vec!["q1".to_string(), "q2".to_string()]
        );
        assert!(reader.next_batch(10).unwrap().is_none());
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

    #[test]
    fn test_batched_sequence_reader_from_file_picks_format_and_decompresses() {
        let dir = std::env::temp_dir();

        let fasta_path = dir.join("raxtax_test_from_file.fasta");
        std::fs::write(&fasta_path, ">q1\nAAAA\n").unwrap();
        let skip = HashSet::new();
        let mut reader = BatchedSequenceReader::from_file(&fasta_path, &skip).unwrap();
        let batch = reader.next_batch(10).unwrap().unwrap();
        assert_eq!(batch, vec![("q1".to_string(), vec![1, 1, 1, 1])]);
        std::fs::remove_file(&fasta_path).unwrap();

        let fastq_gz_path = dir.join("raxtax_test_from_file.fastq.gz");
        let file = std::fs::File::create(&fastq_gz_path).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        encoder.write_all(b"@q1\nCCCC\n+\nIIII\n").unwrap();
        encoder.finish().unwrap();
        let mut reader = BatchedSequenceReader::from_file(&fastq_gz_path, &skip).unwrap();
        let batch = reader.next_batch(10).unwrap().unwrap();
        assert_eq!(batch, vec![("q1".to_string(), vec![2, 2, 2, 2])]);
        std::fs::remove_file(&fastq_gz_path).unwrap();
    }
}
