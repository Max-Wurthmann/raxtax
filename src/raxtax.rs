use anyhow::{ensure, Result};
use std::{path::Path, time::Duration};

use indicatif::{ProgressBar, ProgressStyle};
use log::Level;
use logging_timer::{time, timer};
use rayon::prelude::*;

use crate::io::{Args, ResultsToPrint};
use crate::lineage;
use crate::parser::{BatchedSequenceReader, LabeledSequence};
use crate::tree::{IndexType, KMerMap, Tree};
use crate::{prob, utils};

#[derive(Debug, Clone, Copy)]
pub struct RaxtaxSettings {
    tsv: bool,
    binning: bool,
}

impl RaxtaxSettings {
    pub fn new(tsv: bool, binning: bool) -> Self {
        Self { tsv, binning }
    }
    pub fn from_args(args: &Args) -> RaxtaxSettings {
        RaxtaxSettings {
            tsv: args.tsv,
            binning: args.binning,
        }
    }
}

/// Classifies the queries batch by batch. For every query batch the reference is streamed
/// in batches of `reference_batch_size` sequences (0 = all at once); for each reference batch
/// a k-mer map is built and the query batch, split into `10 * n_threads` chunks, is
/// intersected against it in parallel. After the full pass over the reference, the results of
/// the query batch are final and are sent to the IO-thread.
/// If the whole reference fits into a single batch, its k-mer map is only built once.
#[time("info")]
pub fn raxtax<I>(
    query_batches: I,
    tree: &Tree,
    reference_path: &Path,
    reference_batch_size: usize,
    sender: &crossbeam::channel::Sender<ResultsToPrint>,
    settings: RaxtaxSettings,
    n_queries: usize,
) -> Result<()>
where
    I: Iterator<Item = Result<Vec<LabeledSequence>>>,
{
    let n_references = tree.num_tips;
    let mut kmer_map = KMerMap::new(&tree.encoding_data);
    let is_single_ref_batch = reference_batch_size == 0 || reference_batch_size >= n_references;

    if is_single_ref_batch {
        // compute the one and only k-mer map
        // if !is_single_ref_batch: one k-mer map per reference batch needs to be computed
        let _tmr = timer!(Level::Info; "Creating k-mer map");
        for_each_reference_batch(reference_path, n_references, tree, &mut kmer_map, |_| {})?;
    }

    let pb = if cfg!(test) {
        ProgressBar::hidden()
    } else {
        ProgressBar::new(n_queries as u64)
            .with_style(
                ProgressStyle::with_template(
                    "[{elapsed_precise}] {bar:80.cyan/blue} {pos:>7}/{len:7}[ETA:{eta}] {msg}",
                )
                .unwrap()
                .progress_chars("##-"),
            )
            .with_message("Running Queries...")
    };
    pb.enable_steady_tick(Duration::from_millis(100));

    let n_chunks = 10 * rayon::current_num_threads();

    // intersection sizes of every query in the batch with every reference sequence,
    // query i owns intersections[i * num_tips..(i + 1) * num_tips]
    let mut intersections: Vec<u32> = Vec::new();

    // unique k-mers of every query in the batch, the inner buffers are reused across batches
    let mut query_kmers: Vec<Vec<u32>> = Vec::new();
    for query_batch in query_batches {
        let query_batch = query_batch?;
        let chunk_size = query_batch.len().div_ceil(n_chunks);

        // process query batch into kmers (in parallel)
        query_kmers.resize_with(query_batch.len(), Vec::new);
        query_kmers
            .par_iter_mut()
            .zip(query_batch.par_iter())
            .for_each(|(kmers, query)| {
                utils::seq_to_unique_minenc_canon_kmers_into(
                    &query.sequence,
                    &tree.encoding_data,
                    kmers,
                )
            });

        // zero the intersection buffer in parallel
        // only for the first and the last (smaller batch) is zeroed lazily by the OS
        if intersections.len() == query_batch.len() * n_references {
            intersections
                .par_chunks_mut(chunk_size * n_references)
                .for_each(|chunk| chunk.fill(0));
        } else {
            intersections = vec![0; query_batch.len() * n_references];
        }

        // compute k-mer intersections of the queries in current query batch
        // with all references
        // requires a pass over the whole reference database
        let tmr = timer!(Level::Debug; "K-mer Intersections");
        if is_single_ref_batch {
            accumulate_intersections(
                &kmer_map,
                &query_kmers,
                &mut intersections,
                n_references,
                chunk_size,
            );
        } else {
            // do pass over reverence database in a batched manner
            for_each_reference_batch(
                reference_path,
                reference_batch_size,
                tree,
                &mut kmer_map,
                |kmer_map| {
                    accumulate_intersections(
                        kmer_map,
                        &query_kmers,
                        &mut intersections,
                        n_references,
                        chunk_size,
                    )
                },
            )?;
        }
        drop(tmr);

        // turn intersection sizes for current query batch into confidence and signal values
        // send to output writer thread to write to result file
        let tmr = timer!(Level::Debug; "Hit-prob and result computation");
        query_batch
            .par_chunks(chunk_size)
            .zip(query_kmers.par_chunks(chunk_size))
            .zip(intersections.par_chunks(chunk_size * n_references))
            .try_for_each(|((queries, k_mers), intersections)| -> Result<()> {
                for ((query, k_mers), intersect_buffer) in queries
                    .iter()
                    .zip(k_mers)
                    .zip(intersections.chunks(n_references))
                {
                    let highest_hit_probs = prob::highest_hit_prob_per_reference(
                        k_mers.len() as u32,
                        k_mers.len() / 2,
                        intersect_buffer,
                    );
                    let (eval_res, bin_res) = lineage::Lineage::new(
                        &query.label,
                        tree,
                        highest_hit_probs,
                        settings.binning,
                    )
                    .evaluate();
                    assert!(!eval_res.is_empty());
                    let primary_results = utils::get_results(&eval_res);
                    let tsv_results = if settings.tsv {
                        Some(utils::get_results_tsv(
                            &eval_res,
                            utils::decompress_sequence(&query.sequence),
                        ))
                    } else {
                        None
                    };
                    let binning_result = if settings.binning {
                        Some(utils::get_results_binning(bin_res))
                    } else {
                        None
                    };
                    sender.send(ResultsToPrint::new(
                        query.label.clone(),
                        primary_results,
                        tsv_results,
                        binning_result,
                    ))?;
                    pb.inc(1);
                }
                Ok(())
            })?;
        drop(tmr)
    }

    pb.finish();
    Ok(())
}

/// Streams the reference file in batches of `batch_size` sequences, rebuilds `kmer_map`
/// for each batch and hands it to `f`.
fn for_each_reference_batch(
    reference_path: &Path,
    batch_size: usize,
    tree: &Tree,
    kmer_map: &mut KMerMap,
    mut f: impl FnMut(&KMerMap),
) -> Result<()> {
    let mut n_references = 0;
    for ref_batch in BatchedSequenceReader::from_file(reference_path, batch_size, None)? {
        let ref_batch = ref_batch?;
        ensure!(
            n_references + ref_batch.len() <= tree.num_tips,
            "{} contains more sequences than when it was first parsed. Was it modified during the run?",
            reference_path.display()
        );
        kmer_map.rebuild(&ref_batch, n_references, tree);
        n_references += ref_batch.len();

        f(kmer_map);
    }
    ensure!(
        n_references == tree.num_tips,
        "{} contains fewer sequences than when it was first parsed. Was it modified during the run?",
        reference_path.display()
    );
    Ok(())
}

/// Adds the intersection sizes of every query with the reference sequences of `kmer_map` to
/// `intersections` (query i owns `intersections[i * num_tips..(i + 1) * num_tips]`).
/// The queries are processed in parallel in chunks of `chunk_size`.
fn accumulate_intersections(
    kmer_map: &KMerMap,
    query_kmers: &[Vec<u32>],
    intersections: &mut [u32],
    num_tips: usize,
    chunk_size: usize,
) {
    query_kmers
        .par_chunks(chunk_size)
        .zip(intersections.par_chunks_mut(chunk_size * num_tips))
        .for_each(|(query_kmers, intersections)| {
            for (k_mers, intersect_buffer) in
                query_kmers.iter().zip(intersections.chunks_mut(num_tips))
            {
                for query_kmer in k_mers {
                    kmer_map.kmer_row(*query_kmer as usize).iter().for_each(
                        |sequence_id: &IndexType| {
                            // cast is necessary, because sequence_id is u32 if not cfg(huge_db)
                            // SAFETY: sequence ids are sorted indices into the tree, so < num_tips
                            #[allow(clippy::unnecessary_cast)]
                            unsafe {
                                *intersect_buffer.get_unchecked_mut(*sequence_id as usize) += 1
                            };
                        },
                    );
                }
            }
        });
}

#[cfg(test)]
mod tests {
    use ahash::HashSet;
    use itertools::Itertools;
    use std::io::Cursor;

    use super::{accumulate_intersections, for_each_reference_batch};
    use crate::parser::{parse_reference_taxonomy, FastaSequenceReader, LabeledSequence};
    use crate::tree::KMerMap;
    use crate::utils::{seq_to_unique_minenc_canon_kmers, KMerEncodingData};

    #[test]
    fn test_batched_reference_intersections_match_single_batch() {
        // lineages deliberately not in sorted order, so file order != sorted order
        let reference_str = ">r0;tax=p:B,c:B2;\nACGTTGCAAGCTTAGC\n\
                             >r1;tax=p:A,c:A1;\nACGTTGCAAGGTTAGA\n\
                             >r2;tax=p:C,c:C1;\nTTTTGGGGCCCCAAAA\n\
                             >r3;tax=p:A,c:A1;\nACGTACGTACGTACGT\n\
                             >r4;tax=p:B,c:B1;\nGCTAAGCTTGCAACGT\n";
        let query_str = ">q0\nACGTTGCAAGCTTAGC\n>q1\nTTTTGGGGCCCCAAAT\n>q2\nACGTACGTACGAACGT\n";
        let path = std::env::temp_dir().join("raxtax_test_batched_reference.fasta");
        std::fs::write(&path, reference_str).unwrap();

        let tree = parse_reference_taxonomy(&path, KMerEncodingData::new(5).unwrap(), 5).unwrap();
        let num_tips = tree.num_tips;
        let kmers =
            |sequence: &[u8]| seq_to_unique_minenc_canon_kmers(sequence, &tree.encoding_data);
        let read = |s: &'static str| -> Vec<LabeledSequence> {
            FastaSequenceReader::new(Box::new(Cursor::new(s)))
                .try_collect()
                .unwrap()
        };
        let query_kmers = read(query_str)
            .iter()
            .map(|query| kmers(&query.sequence))
            .collect_vec();

        // brute force: intersect the k-mer sets directly, placing each reference at its sorted index
        let mut expected = vec![0_u32; query_kmers.len() * num_tips];
        for (file_idx, reference) in read(reference_str).iter().enumerate() {
            let reference_kmers: HashSet<u32> = kmers(&reference.sequence).into_iter().collect();
            for (query_idx, q) in query_kmers.iter().enumerate() {
                expected[query_idx * num_tips + tree.sorted_idx[file_idx] as usize] =
                    q.iter().filter(|k| reference_kmers.contains(k)).count() as u32;
            }
        }
        assert!(expected.iter().any(|&c| c > 0));

        for (batch_size, chunk_size) in [(5, 3), (1, 1), (2, 2), (3, 1), (4, 3)] {
            let mut kmer_map = KMerMap::new(&tree.encoding_data);
            let mut intersections = vec![0_u32; query_kmers.len() * num_tips];
            for_each_reference_batch(&path, batch_size, &tree, &mut kmer_map, |kmer_map| {
                accumulate_intersections(
                    kmer_map,
                    &query_kmers,
                    &mut intersections,
                    num_tips,
                    chunk_size,
                )
            })
            .unwrap();
            assert_eq!(intersections, expected, "reference batch size {batch_size}");
        }
        std::fs::remove_file(&path).unwrap();
    }
}
