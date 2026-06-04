use std::{
    collections::HashSet,
    fs::File,
    hint::unreachable_unchecked,
    io::{BufReader, Read},
    path::PathBuf,
};

use anyhow::{bail, Result};
use bitvec::prelude::*;
use flate2::read::GzDecoder;
use itertools::Itertools;
use log::{log_enabled, warn};

use crate::lineage;

pub const F64_OUTPUT_ACCURACY: u32 = 2;

pub fn map_four_to_two_bit_repr(c: u8) -> Option<u8> {
    match c {
        0b0001 => Some(0b00), // A
        0b0010 => Some(0b01), // C
        0b0100 => Some(0b10), // G
        0b1000 => Some(0b11), // T
        _ => None,
    }
}

pub struct BitMasks {
    k: u8,
    unused_bits_mask: u32,
}

impl BitMasks {
    pub fn new(k: u32) -> Option<Self> {
        if (1..=16).step_by(2).contains(&k) {
            return None;
        }

        Some(Self {
            k: k as u8,
            unused_bits_mask: u32::MAX >> (32 - 2 * k),
        })
    }
}

fn encode(kmer: u32, rev_compl_kmer: u32, bitmasks: &BitMasks) -> u32 {
    let k = bitmasks.k;
    let common_prefix_length = (kmer ^ rev_compl_kmer).leading_zeros();
    let mut code = 0u32;

    if common_prefix_length < k as u32 - 1 {
        // specifieing pair

        // TODO:
        unimplemented!()
    } else if common_prefix_length == k as u32 - 1 {
        // only differs in middle char

        // TODO:
        unimplemented!()
    } else {
        debug_assert!(false, "palindrome can not occur with odd k");
        unsafe { unreachable_unchecked() };
    }

    code
}

/// Extracts canonical k-mers from the given sequence.
/// where k is odd and 1 <= k <= 16.
/// Invalid characters are ignored, any k-mer containing an invalid character is skipped.
/// Canonical k-mers are encoded using a minimal encoding scheme and are thus in the range [0, 1/2 * 4^k).
/// (see Wittler R. General encoding of canonical k-mers. Peer Community Journal. 2023)
/// May yield duplicate k-mers, use `seq_to_unique_canon_kmers` to get unique sorted k-mers.
pub fn seq_to_minenc_canon_kmer_iter<'a, 'b>(
    sequence: &'a [u8],
    bitmasks: &'b BitMasks,
) -> impl Iterator<Item = u32> + use<'a, 'b> {
    let mut seq_iter = sequence.iter();
    let mut kmer = 0_u32;
    let mut rev_compl_kmer = 0_u32;
    let mut filled_bases = 0_u8; // Tracks how many consecutive valid bases we have processed

    std::iter::from_fn(move || {
        for char in seq_iter.by_ref() {
            if let Some(repr) = map_four_to_two_bit_repr(*char) {
                // Add repr to the end of kmer
                kmer = (kmer << 2) | repr as u32;
                // Clear any bits that exceed our k-mer size
                kmer &= bitmasks.unused_bits_mask;

                // A <-> T and C <-> G
                let complement_repr = repr as u32 ^ 0b11;
                // Add complement_repr to the start of rev_compl_kmer
                rev_compl_kmer = (rev_compl_kmer >> 2) | (complement_repr << (bitmasks.k * 2 - 2));

                filled_bases += 1;
                if filled_bases >= bitmasks.k {
                    // Only yield a valid k-mer once our sliding window has 8 consecutive valid bases
                    return Some(encode(kmer, rev_compl_kmer, bitmasks));
                }
            } else {
                // Invalid/ambiguous char encountered
                // We completely flush and reset our window state.
                kmer = 0;
                rev_compl_kmer = 0;
                filled_bases = 0;
            }
        }
        // No more characters left in the sequence
        None
    })
}

/// Extracts canonical 8-mers from the given sequence.
/// Invalid characters are ignored, any k-mer containing an invalid character is skipped.
/// May yield duplicate k-mers, use `seq_to_unique_canon_kmers` to get unique sorted k-mers.
pub fn seq_to_canon_kmer_iter(sequence: &[u8]) -> impl Iterator<Item = u16> + use<'_> {
    let mut seq_iter = sequence.iter();
    let mut kmer = 0_u16;
    let mut rev_compl_kmer = 0_u16;
    let mut filled_bases = 0_u32; // Tracks how many consecutive valid bases we have processed

    std::iter::from_fn(move || {
        while let Some(char) = seq_iter.next() {
            if let Some(repr) = map_four_to_two_bit_repr(*char) {
                // Add repr to the end of kmer
                kmer = (kmer << 2) | repr as u16;

                // A <-> T and C <-> G
                let complement_repr = repr as u16 ^ 0b11;
                // Add complement_repr to the start of rev_compl_kmer
                rev_compl_kmer = (rev_compl_kmer >> 2) | (complement_repr << 14);

                filled_bases += 1;
                if filled_bases >= 8 {
                    // Only yield a valid 8-mer once our sliding window has 8 consecutive valid bases
                    return Some(std::cmp::min(kmer, rev_compl_kmer));
                }
            } else {
                // Invalid/ambiguous char encountered
                // We completely flush and reset our window state.
                kmer = 0;
                rev_compl_kmer = 0;
                filled_bases = 0;
            }
        }
        // No more characters left in the sequence
        None
    })
}

/// Extracts all canonical 8-mers from the given sequence.
/// Invalid characters are ignored, any k-mer containing an invalid character is skipped.
/// The resulting k-mers are sorted and unique.
pub fn seq_to_unique_canon_kmers(sequence: &[u8]) -> Vec<u16> {
    // u16::MAX as usize / 32 + 1 == 2048, which is the number of u32 needed to represent all possible u16 values as bits
    let mut bitarr = BitArray::<[u32; u16::MAX as usize / 32 + 1], Msb0>::ZERO;
    seq_to_canon_kmer_iter(sequence).for_each(|canonical_kmer| {
        bitarr.set(canonical_kmer as usize, true);
    });
    bitarr.iter_ones().map(|idx| idx as u16).collect_vec()
}

pub fn seq_to_kmers(sequence: &[u8]) -> Vec<u16> {
    let mut k_mers = HashSet::new();
    sequence.windows(8).for_each(|vals| {
        if let Some(k_mer) = vals
            .iter()
            .enumerate()
            .map(|(j, v)| map_four_to_two_bit_repr(*v).map(|c| (c as u16) << (14 - j * 2)))
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
        cosine_similarity, euclidean_distance_l1, euclidean_norm, seq_to_unique_canon_kmers,
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
    fn test_seq_to_unique_canon_kmers() {
        let check_output = |input_seq: &[u8], kmers_expected: Vec<u16>| {
            let output = seq_to_unique_canon_kmers(input_seq);
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
