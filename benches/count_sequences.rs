use criterion::{criterion_group, criterion_main, Criterion};
use raxtax::parser::count_sequences_in_file;
use std::path::PathBuf;

/// Benchmarks `count_sequences_in_file` against a file named by the
/// `COUNT_SEQ_BENCH_FILE` env var (falling back to the bundled example
/// reference database), and prints the resulting count.
fn bench_count_sequences(c: &mut Criterion) {
    let path = std::env::var("COUNT_SEQ_BENCH_FILE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            println!("The path of the file to count sequences in can be overridden using COUNT_SEQ_BENCH_FILE=/path/to/file");
            println!("No env var set, falling back to example/diptera_references.fasta");
            PathBuf::from("example/diptera_references.fasta")}
        );

    c.bench_function("count_sequences_in_file", |b| {
        b.iter(|| count_sequences_in_file(std::hint::black_box(&path)).unwrap());
    });

    let count = count_sequences_in_file(&path).unwrap();
    println!("Number of Sequences in ({}) = {count}", path.display());
}

criterion_group!(benches, bench_count_sequences);
criterion_main!(benches);
