use ahash::HashSet;
use criterion::{criterion_group, criterion_main, Criterion};
use raxtax::parser::{count_sequences_in_file, parse_reference_taxonomy, BatchedSequenceReader};
use raxtax::tree::{KMerMap, Tree};
use raxtax::utils::KMerEncodingData;
use std::path::{Path, PathBuf};

/// Parses the taxonomy of the reference file and builds the k-mer map over all of its sequences.
fn parse_reference_file(path: &Path, k: u32, n_references: usize) -> (Tree, KMerMap) {
    let tree =
        parse_reference_taxonomy(path, KMerEncodingData::new(k).unwrap(), n_references).unwrap();
    let references = BatchedSequenceReader::from_file(path, &HashSet::default(), tree.num_tips)
        .unwrap()
        .next()
        .unwrap()
        .unwrap();
    let mut kmer_map = KMerMap::new(&tree.encoding_data);
    kmer_map.rebuild(&references, 0, &tree);
    (tree, kmer_map)
}

/// Benchmarks parsing the taxonomy and building the k-mer map of the whole reference against a file named by the
/// `PARSE_REF_BENCH_FILE` env var (falling back to the bundled example
/// reference database), at a k-mer size named by `PARSE_REF_BENCH_K`
/// (falling back to 11), and prints the number of references parsed.
fn bench_parse_reference_file(c: &mut Criterion) {
    let path = std::env::var("PARSE_REF_BENCH_FILE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            println!("The path of the file to parse can be overridden using PARSE_REF_BENCH_FILE=/path/to/file");
            println!("No env var set, falling back to example/diptera_references.fasta");
            PathBuf::from("example/diptera_references.fasta")
        });
    let k: u32 = std::env::var("PARSE_REF_BENCH_K")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or_else(|| {
            println!("The k-mer size can be overridden using PARSE_REF_BENCH_K=<1-16>");
            println!("No env var set, falling back to k=11");
            11
        });

    let n_references = count_sequences_in_file(&path).unwrap();

    c.bench_function("parse_reference_file", |b| {
        b.iter(|| {
            parse_reference_file(
                std::hint::black_box(&path),
                std::hint::black_box(k),
                std::hint::black_box(n_references),
            )
        });
    });

    let (tree, _) = parse_reference_file(&path, k, n_references);
    println!(
        "Parsed {} references from {} at k={k}",
        tree.num_tips,
        path.display()
    );
}

criterion_group!(benches, bench_parse_reference_file);
criterion_main!(benches);
