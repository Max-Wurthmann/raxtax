use std::{
    collections::HashSet,
    fs::File,
    io::{BufReader, Read},
    path::PathBuf,
};

use anyhow::{bail, Result};
use flate2::read::GzDecoder;
use itertools::Itertools;
use log::{log_enabled, warn};

use crate::lineage;

pub const F64_OUTPUT_ACCURACY: u32 = 2;

pub fn map_four_to_two_bit_repr(c: u8) -> Option<u16> {
    match c {
        0b0001 => Some(0b00), // A
        0b0010 => Some(0b01), // C
        0b0100 => Some(0b10), // G
        0b1000 => Some(0b11), // T
        _ => None,
    }
}

pub fn seq_to_canonical_kmers(sequence: &[u8]) -> Vec<u16> {
    let mut k_mers = HashSet::new();
    sequence.windows(8).for_each(|vals| {
        let mut kmer = 0_u16;
        let mut rev_compl_kmer = 0_u16;
        for (i, byte) in vals.iter().enumerate() {
            if let Some(repr) = map_four_to_two_bit_repr(*byte) {
                kmer |= repr << (14 - i * 2);
                let complement_repr = repr ^ 0b11; // A <-> T and C <-> G
                rev_compl_kmer |= complement_repr << (i * 2);
            } else {
                // invalid char, skip this k-mer
                return;
            }
        }
        let canonical_kmer = std::cmp::min(kmer, rev_compl_kmer);
        k_mers.insert(canonical_kmer);
    });
    k_mers.into_iter().sorted().collect_vec()
}

pub fn seq_to_kmers(sequence: &[u8]) -> Vec<u16> {
    let mut k_mers = HashSet::new();
    sequence.windows(8).for_each(|vals| {
        if let Some(k_mer) = vals
            .iter()
            .enumerate()
            .map(|(j, v)| map_four_to_two_bit_repr(*v).map(|c| c << (14 - j * 2)))
            .fold_options(0_u16, |acc, c| acc | c)
        {
            k_mers.insert(k_mer);
        }
    });
    k_mers.into_iter().sorted().collect_vec()
}

pub fn get_reader(path: &PathBuf) -> Result<Box<dyn Read>> {
    let file_type = match path.extension() {
        Some(ext) => match ext.to_str() {
            Some(ext_str) => ext_str.to_ascii_lowercase(),
            None => bail!("Extension could not be parsed!"),
        },
        None => "fasta".to_string(),
    };

    let file = File::open(path)?;

    match file_type.as_str() {
        "gz" | "gzip" => {
            let reader = Box::new(GzDecoder::new(file));
            Ok(Box::new(BufReader::new(reader)))
        }
        _ => Ok(Box::new(BufReader::new(file))),
    }
}

pub fn get_results(results: &[lineage::EvaluationResult<'_, '_>]) -> String {
    results
        .iter()
        .map(lineage::EvaluationResult::get_output_string)
        .collect_vec()
        .join("\n")
}

pub fn decompress_sequence(sequence: &[u8]) -> String {
    sequence
        .iter()
        .map(|c| match c {
            0b0001 => 'A',
            0b0010 => 'C',
            0b0100 => 'G',
            0b1000 => 'T',
            _ => '-',
        })
        .join("")
}

pub fn get_results_tsv(results: &[lineage::EvaluationResult<'_, '_>], sequence: String) -> String {
    results
        .iter()
        .map(|er| er.get_tsv_string(&sequence))
        .collect_vec()
        .join("\n")
}

pub fn get_results_binning(result: Option<(String, f64)>) -> String {
    match result {
        Some((bin, conf)) => format!("{}\t{:.5}", bin, conf),
        None => String::from("NO_BIN\t0.0"),
    }
}

pub fn euclidean_distance_l1(a: &[f64], b: &[f64]) -> f64 {
    assert!(a.len() == b.len());
    if a.is_empty() {
        return 0.0;
    };
    let a_sum = a.iter().sum::<f64>();
    let b_sum = b.iter().sum::<f64>();
    assert!(a_sum > 0.0);
    assert!(b_sum > 0.0);
    a.iter()
        .zip(b)
        .map(|(x, y)| (x / a_sum - y / b_sum).powi(2))
        .sum::<f64>()
        .sqrt()
}

pub fn euclidean_norm<I, T>(v: I) -> f64
where
    I: IntoIterator<Item = T>,
    T: std::borrow::Borrow<f64>,
{
    v.into_iter()
        .map(|x| x.borrow() * x.borrow())
        .sum::<f64>()
        .sqrt()
}

pub fn cosine_similarity(vec_a: &[f64], vec_b: &[f64]) -> f64 {
    let norm_a = euclidean_norm(vec_a.iter());
    let norm_b = euclidean_norm(vec_b.iter());
    assert!(norm_a > 0.0);
    assert!(norm_b > 0.0);
    vec_a
        .iter()
        .zip(vec_b.iter())
        .map(|(a, b)| a * b)
        .sum::<f64>()
        / (norm_a * norm_b)
}

pub fn report_error(e: anyhow::Error, message: impl std::fmt::Display) {
    let prefix = "\x1b[31m[ERROR]\x1b[0m";
    log::error!("{}: {}", message, e);
    if log::log_enabled!(log::Level::Error) {
        eprintln!("{prefix} {message}: {e}");
    }
}

pub fn setup_threadpool_pinned(num_threads: usize) -> Result<()> {
    let cpus = get_thread_ids()?;
    if cpus.len() < num_threads {
        warn!("Only at most {} physical cores are available!", cpus.len());
        if log_enabled!(log::Level::Warn) {
            eprintln!(
                "\x1b[33m[WARN ]\x1b[0m Only at most {} physical cores are available!",
                cpus.len()
            );
        }
    };
    let max_num_threads = num_threads.min(cpus.len());
    rayon::ThreadPoolBuilder::new()
        .num_threads(max_num_threads)
        .start_handler(move |index| {
            core_affinity::set_for_current(core_affinity::CoreId { id: cpus[index] });
        })
        .build_global()?;
    Ok(())
}

pub fn get_thread_ids() -> Result<Vec<usize>> {
    if cfg!(target_os = "linux") {
        let mut used_physical = std::collections::HashSet::new();
        let mut all_cpus = Vec::new();
        let mut preferred_cpus = Vec::new();
        let total_num_cores = core_affinity::get_core_ids().unwrap().len();
        for cpu in 0..total_num_cores {
            let core_id_path = format!("/sys/devices/system/cpu/cpu{}/topology/core_id", cpu);
            let socket_id_path = format!(
                "/sys/devices/system/cpu/cpu{}/topology/physical_package_id",
                cpu
            );
            let core_id = std::fs::read_to_string(core_id_path)?
                .trim()
                .parse::<usize>()?;
            let socket_id = std::fs::read_to_string(socket_id_path)?
                .trim()
                .parse::<usize>()?;

            all_cpus.push((core_id, socket_id, cpu));
        }
        all_cpus
            .into_iter()
            .sorted()
            .for_each(|(core, socket, cpu)| {
                if !used_physical.contains(&(core, socket)) {
                    preferred_cpus.push(cpu);
                    used_physical.insert((core, socket));
                }
            });
        Ok(preferred_cpus)
    } else if let Some(available_cores) = core_affinity::get_core_ids() {
        warn!("Thread-pinning used on non-linux system. Avoiding hyper-threading is not implemented for your platform!");
        Ok(available_cores.into_iter().map(|c| c.id).collect_vec())
    } else {
        anyhow::bail!("Failed to get CPU information!")
    }
}

#[cfg(test)]
mod tests {
    use itertools::assert_equal;
    use statrs::assert_almost_eq;

    use crate::utils::{
        cosine_similarity, euclidean_distance_l1, euclidean_norm, seq_to_canonical_kmers,
    };

    use super::{decompress_sequence, map_four_to_two_bit_repr, seq_to_kmers};

    #[test]
    fn test_euclidean_norm() {
        let v = [1.0, 2.0, 3.0, 4.0];
        assert_almost_eq!(euclidean_norm(v.iter()), 30_f64.sqrt(), 1e-7);
        let w = [0.5, 0.5, 0.25, 0.2];
        assert_almost_eq!(euclidean_norm(w.iter()), 0.6025_f64.sqrt(), 1e-7);
    }

    #[test]
    fn test_euclidean_distance() {
        let v = [1.0, 0.0, 0.0];
        let w = [0.0, 1.0, 0.0];
        assert_almost_eq!(euclidean_distance_l1(&v, &w), 2_f64.sqrt(), 1e-7);
        let x = [0.5, 0.1, 0.1];
        let y = [1.0, 1.0, 0.5];
        assert_almost_eq!(euclidean_distance_l1(&x, &y), 0.410_077_145_554_494_9, 1e-7);
    }

    #[test]
    fn test_cosine_similarity() {
        let v = [1.0, 0.0, 0.0];
        let w = [0.0, 1.0, 0.0];
        assert_almost_eq!(cosine_similarity(&v, &w), 0.0, 1e-7);
        let x = [0.5, 0.5];
        let y = [0.5, 0.5];
        assert_almost_eq!(cosine_similarity(&x, &y), 1.0, 1e-7);
    }

    #[test]
    fn test_map() {
        assert_equal(map_four_to_two_bit_repr(1), Some(0));
        assert_equal(map_four_to_two_bit_repr(2), Some(1));
        assert_equal(map_four_to_two_bit_repr(4), Some(2));
        assert_equal(map_four_to_two_bit_repr(8), Some(3));
        assert!(map_four_to_two_bit_repr(10).is_none());
    }

    #[test]
    fn test_seq_to_canonical_kmers() {
        let check_output = |input_seq: &[u8], kmers_expected: Vec<u16>| {
            let output = seq_to_canonical_kmers(input_seq);
            assert!(output.windows(2).all(|w| w[0] <= w[1]));
            assert_equal(output, kmers_expected);
        };

        // no invalid chars
        let seq1 = vec![1, 2, 1, 4, 8, 2, 8, 4, 1, 4, 8, 2, 8, 4, 1, 4];
        let expected1 = vec![
            0b0001_0010_1101_1110,
            0b0001_1101_0010_0001,
            0b0010_0001_1101_0010,
            0b0010_1101_1110_0010,
            0b0100_1000_0111_0100,
            0b0100_1011_0111_1000,
            0b1000_0111_0100_1000,
            0b1000_1011_0111_1000,
        ];
        let expected2 = expected1.clone();
        check_output(&seq1, expected1);

        // seq2 has invalid chars at front and end
        let seq2 = vec![
            12, 13, 1, 2, 1, 4, 8, 2, 8, 4, 1, 4, 8, 2, 8, 4, 1, 4, 17, 3,
        ];
        // expected2 is same as expected1
        check_output(&seq2, expected2);

        let seq3 = vec![1, 1, 2, 2, 4, 4, 8, 8, 11, 17, 1, 1, 2, 2, 4, 4, 8, 8];
        let expected3 = vec![0b0000_0101_1010_1111];
        check_output(&seq3, expected3);
    }

    #[test]
    fn test_seq_to_kmers() {
        let check_output = |input_seq: &[u8], kmers_expected: Vec<u16>| {
            let output = seq_to_kmers(input_seq);
            assert!(output.windows(2).all(|w| w[0] <= w[1]));
            assert_equal(output, kmers_expected);
        };

        // no invalid chars
        let seq1 = vec![1, 2, 1, 4, 8, 2, 8, 4, 1, 4, 8, 2, 8, 4, 1, 4];
        let expected1: Vec<u16> = vec![
            0b0001_0010_1101_1110,
            0b0010_1101_1110_0010,
            0b0100_1011_0111_1000,
            0b0111_1000_1011_0111,
            0b1000_1011_0111_1000,
            0b1011_0111_1000_1011,
            0b1101_1110_0010_1101,
            0b1110_0010_1101_1110,
        ];
        let expected2 = expected1.clone();
        check_output(&seq1, expected1);

        // seq2 has invalid chars at front and end
        let seq2 = vec![
            12, 13, 1, 2, 1, 4, 8, 2, 8, 4, 1, 4, 8, 2, 8, 4, 1, 4, 17, 3,
        ];
        // expected2 is same as expected1
        check_output(&seq2, expected2);

        let seq3 = vec![1, 1, 2, 2, 4, 4, 8, 8, 11, 17, 1, 1, 2, 2, 4, 4, 8, 8];
        let expected3 = vec![0b0000_0101_1010_1111];
        check_output(&seq3, expected3);
    }

    #[test]
    fn test_decompress_sequence() {
        let sequence = vec![1_u8, 2, 1, 4, 8, 2, 8, 4, 1, 4, 8, 2, 8, 4, 1, 4];
        let decompressed = decompress_sequence(&sequence);
        assert_eq!(decompressed, String::from("ACAGTCTGAGTCTGAG"));
    }
}
