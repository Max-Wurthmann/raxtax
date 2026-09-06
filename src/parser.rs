use ahash::HashSet;
use anyhow::{anyhow, bail, Context, Result};
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use log::{info, log_enabled, warn, Level};
use logging_timer::{time, timer};
use rayon::prelude::*;
use regex::Regex;
use std::time::Duration;
use std::{
    fs::File,
    io::{BufRead, BufReader, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use crate::{tree::Tree, utils::KMerEncodingData};

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

/// Parses a reference FASTA or FASTQ file into a [`Tree`],
/// next to the generated tree, returns a boolean indicating whether the tree was parsed from the file or loaded from a cached tree file. If the cached tree file is present but has a different k-mer size than specified, it will be ignored and the reference file will be reparsed. The function also logs warnings if the k-mer sizes do not match, and returns an error if there are issues reading the file
#[time("info", "Parsing References")]
pub fn parse_reference_fasta_file(
    sequence_path: &PathBuf,
    encoding_data: KMerEncodingData,
    n_references: usize,
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
    let records = SequenceReader::from_file(sequence_path)?;
    Ok((
        true,
        parse_reference_records(records, encoding_data, n_references)?,
    ))
}

/// Consumes a [`SequenceReader`] of labeled reference sequences, splits each
/// label's `tax=<lineage>;<bin>?` annotation off and builds the resulting
/// `Tree`.
fn parse_reference_records(
    records: SequenceReader,
    encoding_data: KMerEncodingData,
    n_references: usize,
) -> Result<Tree> {
    let regex = Regex::new(r"tax=([^;]+);([^;]+)*")?;
    let (labels, sequences): (Vec<LineageBinPair>, Vec<Vec<u8>>) = {
        let _tmr = timer!(Level::Info; "Read file and create k-mer mapping");
        let pb = ProgressBar::new(n_references as u64)
            .with_style(
                ProgressStyle::with_template(
                    "[{elapsed_precise}] {bar:80.cyan/blue} {pos:>7}/{len:7}[ETA:{eta}] {msg}",
                )
                .unwrap()
                .progress_chars("##-"),
            )
            .with_message("Parsing Reference...");
        pb.enable_steady_tick(Duration::from_millis(100));

        let mut labels: Vec<LineageBinPair> = Vec::new();
        let mut sequences: Vec<Vec<u8>> = Vec::new();
        for record in records {
            let (label, sequence) = record?;
            let caps = regex.captures(&label).context(format!(
                "Unexpected taxonomical annotation detected in label {label}"
            ))?;
            let lineage = caps
                .get(1)
                .context(format!("No taxonomic string found in label {label}"))?
                .as_str()
                .to_owned();
            let bin = caps.get(2).map(|bin| bin.as_str().to_owned());
            labels.push((lineage, bin));
            sequences.push(sequence);
            pb.inc(1);
        }
        pb.finish();
        (labels, sequences)
    };
    if labels.is_empty() {
        bail!("File is empty")
    }
    Tree::new(labels, sequences, encoding_data)
}

/// Streams labeled sequences out of a FASTA file, record by record.
pub struct FastaSequenceReader {
    reader: Box<dyn BufRead + Send>,
    next_header: Option<String>,
    line: String,
}

impl FastaSequenceReader {
    pub fn new(reader: Box<dyn BufRead + Send>) -> Self {
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
    reader: Box<dyn BufRead + Send>,
    header: String,
    sequence_line: String,
    plus_line: String,
    quality_line: String,
}

impl FastqSequenceReader {
    pub fn new(reader: Box<dyn BufRead + Send>) -> Self {
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

/// Either a FASTA or FASTQ record reader. Statically dispatches `Iterator::next`
/// to whichever [`SequenceReader::from_file`] picked, avoiding the heap
/// allocation and dynamic dispatch of a `Box<dyn Iterator>`.
pub enum SequenceReader {
    Fasta(FastaSequenceReader),
    Fastq(FastqSequenceReader),
}

impl SequenceReader {
    /// Opens `path` (transparently decompressing `.gz`/`.gzip` files) and
    /// picks a FASTA or FASTQ record reader based on the file extension
    /// (`.fastq`/`.fq` are identified as FASTQ, everything else defaults to FASTA).
    fn from_file(path: &Path) -> Result<Self> {
        let (format, gzipped) = classify_file(path);
        let reader = get_reader(path, gzipped)?;
        Ok(match format {
            FileFormat::Fastq => SequenceReader::Fastq(FastqSequenceReader::new(reader)),
            FileFormat::Fasta => SequenceReader::Fasta(FastaSequenceReader::new(reader)),
        })
    }
}

impl Iterator for SequenceReader {
    type Item = Result<LabeledSequence>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            SequenceReader::Fasta(reader) => reader.next(),
            SequenceReader::Fastq(reader) => reader.next(),
        }
    }
}

/// Batches labeled sequences from an underlying iterator into fixed-size
/// batches, skipping any sequences whose labels are in `seq_to_skip`.
pub struct BatchedSequenceReader<'a> {
    inner: SequenceReader,
    seq_to_skip: &'a HashSet<String>,
    batch_size: usize,
}

impl<'a> BatchedSequenceReader<'a> {
    /// Opens `path` (transparently decompressing `.gz`/`.gzip` files) and
    /// picks a FASTA or FASTQ record reader based on the file extension
    /// (`.fastq`/`.fq` are identified as FASTQ, everything else defaults to FASTA).
    /// Also batches the records into fixed sizes and ingores any sequences whose labels are in
    /// `seq_to_skip`.
    pub fn from_file(
        path: &Path,
        seq_to_skip: &'a HashSet<String>,
        batch_size: usize,
    ) -> Result<Self> {
        let inner = SequenceReader::from_file(path)?;
        Ok(BatchedSequenceReader::new(inner, seq_to_skip, batch_size))
    }

    pub fn new(inner: SequenceReader, seq_to_skip: &'a HashSet<String>, batch_size: usize) -> Self {
        BatchedSequenceReader {
            inner,
            seq_to_skip,
            batch_size,
        }
    }
}

impl<'a> Iterator for BatchedSequenceReader<'a> {
    type Item = Result<Vec<LabeledSequence>>;

    /// Reads and returns up to `batch_size` queries. Returns `None` once the
    /// underlying iterator has been fully consumed.
    fn next(&mut self) -> Option<Self::Item> {
        let mut batch = Vec::with_capacity(self.batch_size);
        while batch.len() < self.batch_size {
            match self.inner.next() {
                Some(Ok((label, sequence))) => {
                    if !self.seq_to_skip.contains(&label) {
                        batch.push((label, sequence));
                    }
                }
                Some(Err(e)) => return Some(Err(e)),
                None => break,
            }
        }
        if batch.is_empty() {
            None
        } else {
            Some(Ok(batch))
        }
    }
}

#[derive(Clone, Copy)]
enum FileFormat {
    Fasta,
    Fastq,
}

fn extension_lowercase(path: &Path) -> String {
    path.extension()
        .and_then(|ext| ext.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase()
}

/// Classifies `path` by its extension, transparently looking through a
/// trailing `.gz`/`.gzip` to the format extension underneath
/// (`.fastq`/`.fq`, or anything else defaulting to FASTA). Returns the
/// format plus whether the file is gzip-compressed.
fn classify_file(path: &Path) -> (FileFormat, bool) {
    let ext = extension_lowercase(path);
    let (ext, gzipped) = match ext.as_str() {
        "gz" | "gzip" => (
            path.file_stem()
                .map(PathBuf::from)
                .map(|stem| extension_lowercase(&stem))
                .unwrap_or_default(),
            true,
        ),
        _ => (ext, false),
    };
    let format = match ext.as_str() {
        "fastq" | "fq" => FileFormat::Fastq,
        "fasta" | "fa" | "fna" | "faa" => FileFormat::Fasta,
        _ => {
            if log_enabled!(Level::Info) {
                eprintln!("[INFO ] Unrecognized file extension {ext}, attempting to parse as FASTA file...");
                info!("Unrecognized file extension {ext}, attempting to parse as FASTA file...");
            }
            FileFormat::Fasta
        }
    };
    (format, gzipped)
}

fn get_reader(path: &Path, gzipped: bool) -> Result<Box<dyn BufRead + Send>> {
    let file = File::open(path)?;
    if gzipped {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Counts occurrences of `target` in `reader`, scanning in large buffered
/// chunks with SIMD-accelerated counting (mirrors the `linecount` crate's
/// `count_lines`, generalized to an arbitrary byte).
fn count_chars_in_file<R: BufRead>(mut reader: R, target: u8) -> Result<usize> {
    let mut count = 0;
    loop {
        let buf = reader.fill_buf()?;
        if buf.is_empty() {
            break;
        }
        count += bytecount::count(buf, target);
        let len = buf.len();
        reader.consume(len);
    }
    Ok(count)
}

/// Quickly estimates the number of sequences (records) in a FASTA/FASTQ file.
/// Transparently handles `.gz`/`.gzip` compression. Plain (uncompressed) files are scanned in
/// parallel across the rayon thread pool while gzipped files are scanned sequentially
pub fn count_sequences_in_file(path: &Path) -> Result<usize> {
    let (format, gzipped) = classify_file(path);
    let target = match format {
        FileFormat::Fasta => b'>',
        FileFormat::Fastq => b'\n',
    };

    let count = if gzipped {
        count_chars_in_file(get_reader(path, true)?, target)?
    } else {
        let file_len = std::fs::metadata(path)?.len();
        let n_chunks = rayon::current_num_threads()
            .min(file_len.max(1) as usize)
            .max(1);

        (0..n_chunks)
            .into_par_iter()
            .map(|i| -> Result<usize> {
                let start = file_len * i as u64 / n_chunks as u64;
                let end = file_len * (i + 1) as u64 / n_chunks as u64;
                let mut file = File::open(path)?;
                file.seek(SeekFrom::Start(start))?;
                count_chars_in_file(BufReader::new(file.take(end - start)), target)
            })
            .try_reduce(|| 0, |a, b| Ok(a + b))?
    };

    match format {
        // if the file is FASTQ, we counted lines and need to divide by 4 to get the number of sequences
        FileFormat::Fastq => Ok(count / 4),
        FileFormat::Fasta => Ok(count),
    }
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
        count_chars_in_file, count_sequences_in_file, parse_reference_fasta_file,
        parse_reference_records, BatchedSequenceReader, FastaSequenceReader, FastqSequenceReader,
        SequenceReader,
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

        let tree = parse_reference_records(
            SequenceReader::Fasta(FastaSequenceReader::new(Box::new(Cursor::new(fasta_str)))),
            KMerEncodingData::new(8).unwrap(),
            6,
        )
        .unwrap();
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
            SequenceReader::Fasta(FastaSequenceReader::new(Box::new(Cursor::new(fasta_str)))),
            &skip,
            10,
        );
        let batch = reader.next().unwrap().unwrap();
        let (_, sequence) = &batch[0];
        assert_eq!(sequence, &[1, 1, 1, 2, 2, 2, 8, 8, 8, 4, 4, 4, 1]);

        let fasta_str2 = r">label1
ACGTWSMKRYBDHVN";
        let mut reader2 = BatchedSequenceReader::new(
            SequenceReader::Fasta(FastaSequenceReader::new(Box::new(Cursor::new(fasta_str2)))),
            &skip,
            10,
        );
        let batch2 = reader2.next().unwrap().unwrap();
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
            SequenceReader::Fasta(FastaSequenceReader::new(Box::new(Cursor::new(fasta_str)))),
            &skip,
            2,
        );

        let batch1 = reader.next().unwrap().unwrap();
        assert_eq!(
            batch1.iter().map(|(l, _)| l.clone()).collect_vec(),
            vec!["q1".to_string(), "q3".to_string()]
        );

        let batch2 = reader.next().unwrap().unwrap();
        assert_eq!(
            batch2.iter().map(|(l, _)| l.clone()).collect_vec(),
            vec!["q4".to_string()]
        );

        let batch3 = reader.next();
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
            SequenceReader::Fastq(FastqSequenceReader::new(Box::new(Cursor::new(fastq_str)))),
            &skip,
            10,
        );
        let batch = reader.next().unwrap().unwrap();
        assert_eq!(
            batch.iter().map(|(l, _)| l.clone()).collect_vec(),
            vec!["q1".to_string(), "q2".to_string()]
        );
        assert!(reader.next().is_none());
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

        let Tree { k_mer_map, .. } = parse_reference_records(
            SequenceReader::Fasta(FastaSequenceReader::new(Box::new(Cursor::new(fasta_str)))),
            KMerEncodingData::new(8).unwrap(),
            5,
        )
        .unwrap();
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
        let mut reader = BatchedSequenceReader::from_file(&fasta_path, &skip, 10).unwrap();
        let batch = reader.next().unwrap().unwrap();
        assert_eq!(batch, vec![("q1".to_string(), vec![1, 1, 1, 1])]);
        std::fs::remove_file(&fasta_path).unwrap();

        let fastq_gz_path = dir.join("raxtax_test_from_file.fastq.gz");
        let file = std::fs::File::create(&fastq_gz_path).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        encoder.write_all(b"@q1\nCCCC\n+\nIIII\n").unwrap();
        encoder.finish().unwrap();
        let mut reader = BatchedSequenceReader::from_file(&fastq_gz_path, &skip, 10).unwrap();
        let batch = reader.next().unwrap().unwrap();
        assert_eq!(batch, vec![("q1".to_string(), vec![2, 2, 2, 2])]);
        std::fs::remove_file(&fastq_gz_path).unwrap();
    }

    #[test]
    fn test_parse_reference_fasta_file_supports_fastq() {
        let dir = std::env::temp_dir();
        let fastq_path = dir.join("raxtax_test_reference.fastq");
        std::fs::write(
            &fastq_path,
            "@ref1;tax=p:Phylum1,c:Class1;\nAAACCCTTTGGGA\n+\nIIIIIIIIIIIII\n",
        )
        .unwrap();
        let (parsed, tree) =
            parse_reference_fasta_file(&fastq_path, KMerEncodingData::new(8).unwrap(), 1).unwrap();
        assert!(parsed);
        assert_eq!(tree.lineages, vec!["p:Phylum1,c:Class1".to_string()]);
        assert_eq!(tree.num_tips, 1);
        std::fs::remove_file(&fastq_path).unwrap();
    }

    #[test]
    fn test_count_chars_in_file_counts_target_byte() {
        let data = b"some\ntext\nwith\nfour\nlines\n".to_vec();
        let count = count_chars_in_file(Cursor::new(data), b'\n').unwrap();
        assert_eq!(count, 5);
    }

    #[test]
    fn test_count_sequences_in_file_plain_fasta() {
        let dir = std::env::temp_dir();
        let fasta_path = dir.join("raxtax_test_count_sequences.fasta");
        std::fs::write(&fasta_path, ">q1\nAAAA\n>q2\nCCCC\n>q3\nGGGG\n>q4\nTTTT\n").unwrap();
        assert_eq!(count_sequences_in_file(&fasta_path).unwrap(), 4);
        std::fs::remove_file(&fasta_path).unwrap();
    }

    #[test]
    fn test_count_sequences_in_file_plain_fastq() {
        let dir = std::env::temp_dir();
        let fastq_path = dir.join("raxtax_test_count_sequences.fastq");
        std::fs::write(&fastq_path, "@q1\nAAAA\n+\nIIII\n@q2\nCCCC\n+\nIIII\n").unwrap();
        assert_eq!(count_sequences_in_file(&fastq_path).unwrap(), 2);
        std::fs::remove_file(&fastq_path).unwrap();
    }

    #[test]
    fn test_count_sequences_in_file_gzipped_fastq() {
        let dir = std::env::temp_dir();
        let fastq_gz_path = dir.join("raxtax_test_count_sequences.fastq.gz");
        let file = std::fs::File::create(&fastq_gz_path).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        encoder
            .write_all(b"@q1\nAAAA\n+\nIIII\n@q2\nCCCC\n+\nIIII\n@q3\nGGGG\n+\nIIII\n")
            .unwrap();
        encoder.finish().unwrap();
        assert_eq!(count_sequences_in_file(&fastq_gz_path).unwrap(), 3);
        std::fs::remove_file(&fastq_gz_path).unwrap();
    }

    #[test]
    fn test_count_sequences_in_file_tiny_file() {
        let dir = std::env::temp_dir();
        let fasta_path = dir.join("raxtax_test_count_sequences_tiny.fasta");
        std::fs::write(&fasta_path, ">q\nA\n").unwrap();
        assert_eq!(count_sequences_in_file(&fasta_path).unwrap(), 1);
        std::fs::remove_file(&fasta_path).unwrap();

        let empty_path = dir.join("raxtax_test_count_sequences_empty.fasta");
        std::fs::write(&empty_path, "").unwrap();
        assert_eq!(count_sequences_in_file(&empty_path).unwrap(), 0);
        std::fs::remove_file(&empty_path).unwrap();
    }
}
