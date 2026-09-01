use anyhow::Result;
use std::time::Duration;

use indicatif::{ProgressBar, ProgressStyle};
use log::Level;
use logging_timer::{time, timer};
use rayon::prelude::*;

use crate::io::{Args, ResultsToPrint};
use crate::lineage;
use crate::tree::{IndexType, Tree};
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

#[time("info")]
pub fn raxtax<I>(
    batches: I,
    tree: &Tree,
    sender: &crossbeam::channel::Sender<ResultsToPrint>,
    settings: RaxtaxSettings,
) -> Result<()>
where
    I: Iterator<Item = Result<Vec<(String, Vec<u8>)>>> + Send,
{
    let pb = ProgressBar::new_spinner()
        .with_style(
            ProgressStyle::with_template("[{elapsed_precise}] {pos:>9} queries [{per_sec}] {msg}")
                .unwrap(),
        )
        .with_message("Running Queries...");
    pb.enable_steady_tick(Duration::from_millis(100));

    // One parallel pipeline over the whole query stream: `par_bridge` pulls the
    // next batch from the sequential reader on a free worker (overlapping
    // parsing/decompression with classification) instead of a serial loop of
    // fork-join batches. `try_for_each_init` hands each worker a reusable
    // intersection buffer (intersection size, as u16, with each reference).
    batches.par_bridge().try_for_each_init(
        || vec![0u16; tree.num_tips],
        |intersect_buffer: &mut Vec<u16>, batch| -> Result<()> {
            for (query_label, query_sequence) in &batch? {
                pb.inc(1);
                intersect_buffer.fill(0);
                let tmr = timer!(Level::Debug; "K-mer Intersections");
                let k_mers =
                    utils::seq_to_unique_minenc_canon_kmers(query_sequence, &tree.encoding_data);
                assert!(u16::try_from(k_mers.len()).is_ok());
                let num_trials = k_mers.len() / 2;
                for query_kmer in &k_mers {
                    tree.k_mer_map[*query_kmer as usize].iter().for_each(
                        |sequence_id: &IndexType| {
                            unsafe {
                                *intersect_buffer.get_unchecked_mut(*sequence_id as usize) += 1
                            };
                        },
                    );
                }
                drop(tmr);
                let highest_hit_probs = prob::highest_hit_prob_per_reference(
                    k_mers.len() as u16,
                    num_trials,
                    intersect_buffer,
                );
                let (eval_res, bin_res) =
                    lineage::Lineage::new(query_label, tree, highest_hit_probs, settings.binning)
                        .evaluate();
                assert!(!eval_res.is_empty());
                let primary_results = utils::get_results(&eval_res);
                let tsv_results = if settings.tsv {
                    Some(utils::get_results_tsv(
                        &eval_res,
                        utils::decompress_sequence(query_sequence),
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
                    query_label.clone(),
                    primary_results,
                    tsv_results,
                    binning_result,
                ))?;
            }
            Ok(())
        },
    )?;

    pb.finish();
    Ok(())
}
