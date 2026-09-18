use criterion::{criterion_group, criterion_main, Criterion};
use raxtax::parser::{count_sequences_in_file, parse_reference_file};
use raxtax::utils::KMerEncodingData;
use std::path::PathBuf;

/// Benchmarks `parse_reference_file` against a file named by the
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
                std::hint::black_box(KMerEncodingData::new(k).unwrap()),
                std::hint::black_box(n_references),
            )
            .unwrap()
        });
    });

    let (parsed_from_file, _tree) =
        parse_reference_file(&path, KMerEncodingData::new(k).unwrap(), n_references).unwrap();
    println!(
        "Parsed {n_references} references from {} at k={k} (freshly parsed: {parsed_from_file})",
        path.display()
    );
}

criterion_group!(benches, bench_parse_reference_file);
criterion_main!(benches);
