use std::process::exit;

use anyhow::{anyhow, Result};
use clap::Parser;
use log::Level;
use logging_timer::timer;
use raxtax::io::FileFingerprint;
use raxtax::io::{self, ResultsToPrint};
use raxtax::parser::{count_sequences_in_file, parse_reference_fasta_file, BatchedSequenceReader};
use raxtax::raxtax::{raxtax, RaxtaxSettings};
use raxtax::utils::{self, KMerEncodingData};
use std::io::Write;

fn main() {
    // Parse args, set up files and other context
    let args = io::Args::parse();
    let (
        io::OutputWriters {
            primary: mut output,
            tsv: mut tsv_output,
            binning: mut binning_output,
            log: mut log_output,
            progress: mut progress_output,
        },
        mut checkpoint,
    ) = args.get_output().unwrap_or_else(|e| {
        if args.verbosity.log_level_filter() >= Level::Error {
            eprintln!("\x1b[31m[ERROR]\x1b[0m {e}");
        }
        exit(exitcode::CANTCREAT);
    });
    if let Err(e) = io::write_build_info(&mut log_output) {
        eprintln!("\x1b[31m[ERROR]\x1b[0m {e}");
    }
    env_logger::Builder::new()
        .target(env_logger::Target::Pipe(log_output))
        .filter_level(args.verbosity.log_level_filter())
        .format_timestamp(None)
        .format_target(false)
        .init();
    if args.pin {
        if let Err(e) = utils::setup_threadpool_pinned(args.threads) {
            utils::report_error(e, "Failed to set up thread pinning! Continuing without");
            if let Err(e) = rayon::ThreadPoolBuilder::new()
                .num_threads(args.threads)
                .build_global()
            {
                utils::report_error(anyhow::Error::from(e), "Failed to set up Multithreading");
                exit(exitcode::OSERR);
            }
        };
    } else if let Err(e) = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
    {
        utils::report_error(anyhow::Error::from(e), "Failed to set up Multithreading");
        exit(exitcode::OSERR);
    };
    let _total_tmr = timer!(Level::Info; "Total Runtime");

    let encoding_data = KMerEncodingData::new(args.kmer_size).unwrap_or_else(|e| {
        utils::report_error(e, "Please provide a valid k-mer size.");
        exit(exitcode::USAGE);
    });

    // Count number of reference sequences
    let n_references = count_sequences_in_file(&args.database_path).unwrap_or_else(|e| {
        utils::report_error(
            e,
            format!(
                "Failed to count sequences in {}",
                args.database_path.display()
            ),
        );
        exit(exitcode::NOINPUT);
    });
    let (store_db, tree) =
        parse_reference_fasta_file(&checkpoint.db_fingerprint.path, encoding_data, n_references)
            .unwrap_or_else(|e| {
                utils::report_error(
                    e,
                    format!(
                        "Failed to parse {}",
                        checkpoint.db_fingerprint.path.display()
                    ),
                );
                exit(exitcode::NOINPUT);
            });
    if store_db && !args.skip_db {
        match args.get_db_output() {
            Ok((db_output, db_path)) => {
                if let Err(e) = tree.save_to_file(db_output) {
                    utils::report_error(e, "Failed to write database");
                    exit(exitcode::IOERR);
                };
                if let Err(e) = FileFingerprint::new(&db_path).and_then(|fp| {
                    checkpoint.db_fingerprint = fp;
                    checkpoint.save()?;
                    Ok(())
                }) {
                    utils::report_error(e, "Failed to write checkpoint! Continuing without...");
                };
            }
            Err(e) => {
                utils::report_error(
                    e,
                    "Could not create database! Rerun with --skip-db to skip this step.",
                );
                exit(exitcode::CANTCREAT);
            }
        }
    } else {
        checkpoint.save().unwrap_or_else(|e| {
            utils::report_error(e, "Failed to write checkpoint! Continuing without...")
        });
    }

    if args.only_db {
        exit(exitcode::OK);
    }
    let settings = RaxtaxSettings::from_args(&args);

    let n_threads = rayon::current_num_threads();

    let query_file = args.query_file.clone().unwrap();
    let n_queries = count_sequences_in_file(&query_file).unwrap_or_else(|e| {
        utils::report_error(
            e,
            format!("Failed to count sequences in {}", query_file.display()),
        );
        exit(exitcode::NOINPUT);
    });
    let n_queries_remaining = n_queries.saturating_sub(checkpoint.processed_queries.len());

    let batch_size = {
        let mut b_size = if args.query_batch_size == 0 {
            n_queries_remaining
        } else {
            args.query_batch_size.min(n_queries_remaining)
        };

        if n_threads > 1 {
            // dividy by n_threads to properly distribute the work across threads
            // divide by 10 to improve load balancing
            b_size /= 10 * n_threads
        }

        b_size.max(100)
    };

    let (sender, receiver) = crossbeam::channel::unbounded::<ResultsToPrint>();
    let writer_handle = std::thread::spawn(move || -> Result<()> {
        for ResultsToPrint {
            query,
            primary,
            tsv,
            binning,
        } in receiver
        {
            if let Some(ref mut tsv_output) = tsv_output {
                writeln!(tsv_output, "{}", tsv.unwrap())?;
            }
            if let Some(ref mut binning_output) = binning_output {
                writeln!(binning_output, "{}\t{}", query, binning.unwrap())?;
            }
            writeln!(output, "{}", primary)?;
            writeln!(progress_output, "{}", query)?;
        }
        Ok(())
    });

    let query_reader =
        BatchedSequenceReader::from_file(&query_file, &checkpoint.processed_queries, batch_size)
            .unwrap_or_else(|e| {
                utils::report_error(e, format!("Failed to open {}", query_file.display()));
                exit(exitcode::NOINPUT);
            });

    if let Err(e) = raxtax(query_reader, &tree, &sender, settings, n_queries_remaining) {
        if e.is::<crossbeam::channel::SendError<ResultsToPrint>>() {
            utils::report_error(
                e,
                "Error while sending results to IO-thread!\n
                        Rerun raxtax to continue from the last checkpoint.\n
                        If the problem persists, please report this issue at: https://github.com/noahares/raxtax/issues",
            );
            exit(exitcode::TEMPFAIL);
        }
        utils::report_error(e, format!("Failed to parse {}", query_file.display()));
        exit(exitcode::NOINPUT);
    }

    drop(sender);
    if writer_handle.join().is_err() {
        utils::report_error(
            anyhow!("IO-thread could not be joined. Check if results are complete!"),
            "",
        );
    };

    if args.clean {
        checkpoint.cleanup().unwrap_or_else(|e| {
            utils::report_error(
                e,
                "Removing checkpoint files failed! Please delete them manually.",
            )
        });
    }

    exit(exitcode::OK);
}
