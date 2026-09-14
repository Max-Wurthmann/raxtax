use bitvec::prelude::*;
use criterion::{criterion_group, criterion_main, Criterion};
use itertools::Itertools;
use raxtax::parser::SequenceReader;
use raxtax::utils::{
    seq_to_minenc_canon_kmer_iter, seq_to_unique_minenc_canon_kmers, KMerEncodingData,
};
use std::collections::HashSet;
use std::path::PathBuf;

// dedup via hashset
fn seq_to_unique_minenc_canon_kmers_hashset(
    sequence: &[u8],
    encoding_data: &KMerEncodingData,
) -> Vec<u32> {
    let mut hash_set = HashSet::<u32>::new();
    seq_to_minenc_canon_kmer_iter(sequence, encoding_data).for_each(|kmer_code| {
        hash_set.insert(kmer_code);
    });
    hash_set.into_iter().sorted().collect_vec()
}

/// dedup via bitarray, requires possibly large allocation
fn seq_to_unique_minenc_canon_kmers_bitvec(
    sequence: &[u8],
    encoding_data: &KMerEncodingData,
) -> Vec<u32> {
    let mut bitarr = bitvec![0; encoding_data.n_unique_codes as usize];
    seq_to_minenc_canon_kmer_iter(sequence, encoding_data).for_each(|kmer_code| {
        bitarr.set(kmer_code as usize, true);
    });
    bitarr.iter_ones().map(|idx| idx as u32).collect()
}

/// dedup via inplace sort + dedup
fn seq_to_unique_minenc_canon_kmers_inplace(
    sequence: &[u8],
    encoding_data: &KMerEncodingData,
) -> Vec<u32> {
    let mut kmers = seq_to_minenc_canon_kmer_iter(sequence, encoding_data).collect_vec();

    kmers.sort_unstable();
    kmers.dedup();
    kmers.shrink_to_fit();

    kmers
}

/// Benchmarks three dedup strategies for `seq_to_unique_minenc_canon_kmers` (hashset,
/// bitvec, inplace sort+dedup) against the first sequence of a file named by the
/// `KMER_DEDUP_BENCH_FILE` env var (falling back to the bundled example reference
/// database), at a k-mer size named by `KMER_DEDUP_BENCH_K` (falling back to 11), and
/// prints the resulting unique k-mer counts for a manual sanity check.
fn bench_seq_to_unique_minenc_canon_kmers(c: &mut Criterion) {
    let path = std::env::var("KMER_DEDUP_BENCH_FILE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            println!(
                "The path of the file to pull a sequence from can be overridden using KMER_DEDUP_BENCH_FILE=/path/to/file"
            );
            println!("No env var set, falling back to example/diptera_references.fasta");
            PathBuf::from("example/diptera_references.fasta")
        });
    let k: u32 = std::env::var("KMER_DEDUP_BENCH_K")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or_else(|| {
            println!("The k-mer size can be overridden using KMER_DEDUP_BENCH_K=<1-16>");
            println!("No env var set, falling back to k=11");
            11
        });
    let encoding_data = KMerEncodingData::new(k).unwrap();

    let sequence = SequenceReader::from_file(&path)
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .1;

    println!(
        "Using first sequence from {} ({} bases) at k={k}",
        path.display(),
        sequence.len()
    );

    c.bench_function(
        "seq_to_unique_minenc_canon_kmers (inplace sort+dedup, current)",
        |b| {
            b.iter(|| {
                seq_to_unique_minenc_canon_kmers_inplace(
                    std::hint::black_box(&sequence),
                    std::hint::black_box(&encoding_data),
                )
            });
        },
    );

    c.bench_function("seq_to_unique_minenc_canon_kmers (hashset)", |b| {
        b.iter(|| {
            seq_to_unique_minenc_canon_kmers_hashset(
                std::hint::black_box(&sequence),
                std::hint::black_box(&encoding_data),
            )
        });
    });

    c.bench_function("seq_to_unique_minenc_canon_kmers (bitvec)", |b| {
        b.iter(|| {
            seq_to_unique_minenc_canon_kmers_bitvec(
                std::hint::black_box(&sequence),
                std::hint::black_box(&encoding_data),
            )
        });
    });

    let sorted_result = seq_to_unique_minenc_canon_kmers(&sequence, &encoding_data);
    let hashset_result = seq_to_unique_minenc_canon_kmers_hashset(&sequence, &encoding_data);
    let bitvec_result = seq_to_unique_minenc_canon_kmers_bitvec(&sequence, &encoding_data);
    println!(
        "unique k-mers: inplace sort+dedup={}, hashset={}, bitvec={} (all match: {})",
        sorted_result.len(),
        hashset_result.len(),
        bitvec_result.len(),
        sorted_result == hashset_result && sorted_result == bitvec_result
    );
}

criterion_group!(benches, bench_seq_to_unique_minenc_canon_kmers);
criterion_main!(benches);
