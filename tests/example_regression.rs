//! Regression test on the example data: compares `raxtax.out` against a committed baseline.
//!
//! Results are not byte-reproducible (output order and tie-breaking between equally good
//! references vary), so the lineage column is ignored: per query, the confidences and signals
//! must match exactly, while lineage-only differences are just reported.
//!
//! Run:   cargo test --release --test example_regression -- --ignored
//! Bless: RAXTAX_BLESS=1 cargo test --release --test example_regression -- --ignored

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::Command;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

const BASELINE: &str = "tests/data/example_k7.out.gz";

/// Elements only in `a` and only in `b`, for sorted inputs (multiset semantics).
fn sorted_diff<'a>(
    a: &'a [String],
    b: &'a [String],
    without_lineage: bool,
) -> (Vec<&'a String>, Vec<&'a String>) {
    assert!(a.is_sorted() && b.is_sorted());

    let (mut i, mut j) = (0, 0);
    let (mut only_a, mut only_b) = (Vec::new(), Vec::new());
    while i < a.len() && j < b.len() {
        let comparison = if without_lineage {
            let a_stripped = remove_lineage(&a[i]);
            let b_stripped = remove_lineage(&b[j]);
            a_stripped.cmp(&b_stripped)
        } else {
            a[i].cmp(&b[j])
        };
        match comparison {
            std::cmp::Ordering::Less => {
                only_a.push(&a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                only_b.push(&b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    only_a.extend(&a[i..]);
    only_b.extend(&b[j..]);
    (only_a, only_b)
}

fn remove_lineage(line: &str) -> String {
    let mut parts = line.split('\t');
    // String consists of query label, linage, confidence vector, signal1, signal2
    assert!(
        parts.clone().count() == 5,
        "line does not have 5 tab-separated fields: {line}"
    );
    let label = parts.next().unwrap();
    let _ = parts.next().unwrap();
    let metrics = parts.collect::<Vec<_>>().join("\t");
    format!("{label}\t{metrics}")
}

#[test]
#[ignore = "slow; run with `cargo test --release --test example_regression -- --ignored`"]
fn example_regression() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("example_regression_out");

    let output = Command::new(env!("CARGO_BIN_EXE_raxtax"))
        .arg("-k")
        .arg("7")
        .arg("-d")
        .arg(root.join("example/diptera_references.fasta"))
        .arg("-i")
        .arg(root.join("example/diptera_queries.fasta"))
        .arg("-o")
        .arg(&out_dir)
        .args(["--skip-db", "--redo"])
        .output()
        .expect("failed to run raxtax");
    assert!(
        output.status.success(),
        "raxtax exited with {}:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );

    let mut actual: Vec<String> = BufReader::new(File::open(out_dir.join("raxtax.out")).unwrap())
        .lines()
        .collect::<Result<_, _>>()
        .unwrap();
    actual.sort();

    // delete the output directory
    if out_dir
        .canonicalize()
        .ok()
        .and_then(|p| std::fs::remove_dir_all(p).ok())
        .is_none()
    {
        eprintln!(
            "Failed to remove temporary output directory {}",
            out_dir.display()
        );
    }

    let baseline_path = root.join(BASELINE);
    if std::env::var_os("RAXTAX_BLESS").is_some() {
        // write the current output as the new baseline
        std::fs::create_dir_all(baseline_path.parent().unwrap()).unwrap();
        let mut enc = GzEncoder::new(File::create(&baseline_path).unwrap(), Compression::best());
        for line in &actual {
            writeln!(enc, "{line}").unwrap();
        }
        enc.finish().unwrap();
        println!(
            "Wrote baseline with {} lines to {}",
            actual.len(),
            baseline_path.display()
        );
        return;
    }

    let mut expected: Vec<String> = BufReader::new(GzDecoder::new(
        File::open(&baseline_path).expect("baseline missing; run with RAXTAX_BLESS=1"),
    ))
    .lines()
    .collect::<Result<_, _>>()
    .unwrap();
    expected.sort();

    let (only_actual_stripped, only_expected_stripped) = sorted_diff(&actual, &expected, true);
    if !only_actual_stripped.is_empty() || !only_expected_stripped.is_empty() {
        println!("Results differ from baseline {}:", baseline_path.display());
        println!(
            "Results only in baseline {}, Results only in actual {}",
            only_expected_stripped.len(),
            only_actual_stripped.len()
        );
        println!("If this is expected, run with RAXTAX_BLESS=1 to update the baseline.");

        println!("Differences (up to 10):\n");
        only_actual_stripped
            .iter()
            .zip(only_expected_stripped.iter())
            .take(10)
            .for_each(|(&a, &e)| {
                println!("actual:   {a}\nexpected: {e}\n");
            });

        panic!("confidences/signals differ from baseline",);
    }

    let (lineage_only, _) = sorted_diff(&actual, &expected, false);

    // regular println! does not print on successful test runs
    let _ = writeln!(
        std::io::stdout(),
        "{} out of {} lines differ in assigned lineage but are identical confidence values, this can safely be ignored as tie-breaking is non-deterministic",
        lineage_only.len(),
        actual.len()
    );
}
