use std::collections::HashSet;

use anyhow::{bail, Result};
use bitvec::prelude::*;
use itertools::Itertools;
use log::{log_enabled, warn};
use serde::{Deserialize, Serialize};

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

// replace patterns for specifieng pair
const REPLACE: [u8; 16] = [
    0x06, 0x05, 0x04, 0x00, 0x08, 0x07, 0x00, 0x04, 0x09, 0x00, 0x07, 0x05, 0x00, 0x09, 0x08, 0x06,
];
// Markers to check if we need to encode the forward or the reverse complement.
const REVERSE: [bool; 16] = [
    false, false, false, false, false, false, false, true, false, false, true, true, false, true,
    true, true,
];

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct KMerEncodingData {
    pub k: u8,
    pub n_unique_codes: u32,
    unused_bits_mask: u32,
    four_to_the_k_half_plus_one: u32,
    twice_four_to_the_k_half: u32,
    remaindermasks: [u32; 19],
}

impl KMerEncodingData {
    pub fn new(k: u32) -> Result<Self> {
        if !(1..=16).contains(&k) {
            bail!("k-mer size k must satisfy 1 <= k <= 16.");
        }

        // size 19 as k + 2 is largest index we need
        // k <= 16 implies k + 2 < 19
        let mut remaindermasks = [0u32; 19];

        remaindermasks[0] = u32::MAX;
        for i in 1..=k {
            let zero_shift = 32 - 2 * k + i;
            let zeromask = u32::MAX >> zero_shift;
            let onemask = u32::MAX << i;
            remaindermasks[i as usize] = zeromask & onemask;
        }

        let n_unique_codes = if k % 2 == 1 {
            // 4u32.pow(k) / 2
            1 << (2 * k - 1)
        } else {
            // (4u64.pow(k) / 2 + 4u64.pow(k / 2) / 2) as u32
            (1 << (2 * k - 1)) + (1 << (k - 1))
        };

        Ok(Self {
            k: k as u8,
            n_unique_codes,
            unused_bits_mask: u32::MAX >> (32 - 2 * k),
            four_to_the_k_half_plus_one: 1 << ((k / 2) * 2 + 2),
            twice_four_to_the_k_half: 1 << ((k / 2) * 2 + 1),
            remaindermasks,
        })
    }
}

/// Encodes the canonical k-mer specified by the given k-mer and its reverse complement
/// Canonical k-mers are mapped to the range [0, 1/2 * 4^k).
/// Uses the encoding scheme described in
/// (Wittler R. General encoding of canonical k-mers. Peer Community Journal. 2023.0)
pub fn encode(kmer: u32, rev_compl_kmer: u32, encoding_data: &KMerEncodingData) -> u32 {
    let k = encoding_data.k;
    // number of common bits rounded down to nearest even number
    let common_prefix_length = ((kmer ^ rev_compl_kmer).trailing_zeros() / 2 * 2) as u8;

    let mut kmer_code = if common_prefix_length > k {
        // is p palindrome
        debug_assert!(k.is_multiple_of(2), "palindrome can not occur with odd k");
        debug_assert!(
            common_prefix_length == 32,
            "palindrome should have common prefix length of 32"
        );
        // use k as common prefix length as actual common prefix length is 32 > k
        encode_prime(kmer, k, encoding_data)
    } else if common_prefix_length < k - 1 {
        // specifieing pair

        // determine which k-mer is lexicographically smaller
        // via a lookup of the specifieng pair
        let mut pattern = 0u8;
        pattern |= (kmer >> (2 * k - common_prefix_length - 4)) as u8 & 0x0C;
        pattern |= (kmer >> common_prefix_length) as u8 & 0x03;

        debug_assert!(pattern < 16);

        let mut kmer_code = if REVERSE[pattern as usize] {
            encode_prime(rev_compl_kmer, common_prefix_length, encoding_data)
        } else {
            encode_prime(kmer, common_prefix_length, encoding_data)
        };

        // insert replace pattern
        kmer_code |= (REPLACE[pattern as usize] as u32) << (2 * k - common_prefix_length - 4);
        kmer_code
    } else {
        // common_prefix_length == k - 1 as
        // - !common_prefix_length < k - 1 and
        // - !common_prefix_length >= k (as this would be a palindrome, see first case)
        debug_assert!(common_prefix_length == k - 1);

        // only differs in middle char
        // look at two bits of middle char, located at positions k and k-1 (0-indexed from the right)
        // Use these bits to encode A/T -> 0 and C/G -> 1:
        let mut kmer_code = encode_prime(kmer, common_prefix_length, encoding_data);

        let bit1 = (kmer & (1 << k)) >> 1;
        let bit2 = kmer & (1 << (k - 1));

        kmer_code |= bit1 ^ bit2;
        kmer_code
    };

    // subtract gaps
    // 2*(k//2-common_prefix_length-1) ones followed by k-2 zeros
    if k >= 4 && common_prefix_length <= k - 4 {
        let gaps = u32::MAX >> (32 - (k / 2 * 2 - common_prefix_length - 2));

        if k % 2 == 1 {
            kmer_code -= gaps << k;
        } else {
            kmer_code -= gaps << (k - 1);
        }
    }

    // subtract gap in code due to specifying middle position
    if k % 2 == 1 && kmer_code >= encoding_data.four_to_the_k_half_plus_one {
        kmer_code -= encoding_data.twice_four_to_the_k_half;
    }

    kmer_code
}

/// encodes a given canonical k-mer which has the given common prefix length with its reverse complement
/// does not subtract gaps
/// is used as a subroutine of `encode`
fn encode_prime(kmer: u32, common_prefix_length: u8, encoding_data: &KMerEncodingData) -> u32 {
    let k = encoding_data.k;

    // This uses a mask of the form 0..01..1 (common_prefix_length trailing ones), to extract
    // the relevant bits on the right, and invert (complement) them.
    let zeromask = u32::MAX.unbounded_shr(32 - common_prefix_length as u32);
    let right = (kmer & zeromask) ^ zeromask;

    // Assert that the values are as expected.
    debug_assert!(common_prefix_length <= k);
    debug_assert!(common_prefix_length.is_multiple_of(2));

    // Use the remainder mask (consisting of ones in the middle) to extract the bits
    // in between the specifying pair, then shift the remainder to the correct position.
    // The mask contains 0 after index k, so that if we have l+2 >= k (no remainder),
    // we just get a zero here, which does nothing to our result.
    let remainder = (kmer & encoding_data.remaindermasks[common_prefix_length as usize + 2]) >> 2;
    right | remainder
}

/// Extracts canonical k-mers from the given sequence.
/// where k is odd and 1 <= k <= 16.
/// Invalid characters are ignored, any k-mer containing an invalid character is skipped.
/// Canonical k-mers are encoded using a minimal encoding scheme and are thus in the range [0, 1/2 * 4^k).
/// (see Wittler R. General encoding of canonical k-mers. Peer Community Journal. 2023)
/// May yield duplicate k-mers, use `seq_to_unique_canon_kmers` to get unique sorted k-mers.
pub fn seq_to_minenc_canon_kmer_iter<'a, 'b>(
    sequence: &'a [u8],
    encoding_data: &'b KMerEncodingData,
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
                kmer &= encoding_data.unused_bits_mask;

                // A <-> T and C <-> G
                let complement_repr = repr as u32 ^ 0b11;
                // Add complement_repr to the start of rev_compl_kmer
                rev_compl_kmer =
                    (rev_compl_kmer >> 2) | (complement_repr << (encoding_data.k * 2 - 2));

                filled_bases += 1;
                if filled_bases >= encoding_data.k {
                    // Only yield a valid k-mer once our sliding window has 8 consecutive valid bases
                    return Some(encode(kmer, rev_compl_kmer, encoding_data));
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

/// Extracts all canonical k-mers from the given sequence, where k is odd and 1 <= k <= 16.
/// Invalid characters are ignored, any k-mer containing an invalid character is skipped.
/// Canonical k-mers are encoded using a minimal encoding scheme and are thus in the range [0, 1/2 * 4^k).
/// (see Wittler R. General encoding of canonical k-mers. Peer Community Journal. 2023)
/// The resulting k-mers are sorted and unique.
pub fn seq_to_unique_minenc_canon_kmers(
    sequence: &[u8],
    encoding_data: &KMerEncodingData,
) -> Vec<u32> {
    // let mut bitvec = bitvec![0; encoding_data.n_unique_codes as usize];
    let mut hash_set = HashSet::<u32>::new();
    seq_to_minenc_canon_kmer_iter(sequence, encoding_data).for_each(|kmer_code| {
        hash_set.insert(kmer_code);
    });
    hash_set.into_iter().sorted().collect_vec()
}

/// Extracts canonical 8-mers from the given sequence.
/// Invalid characters are ignored, any k-mer containing an invalid character is skipped.
/// May yield duplicate k-mers, use `seq_to_unique_canon_kmers` to get unique sorted k-mers.
pub fn seq_to_canon_8mer_iter(sequence: &[u8]) -> impl Iterator<Item = u16> + use<'_> {
    let mut seq_iter = sequence.iter();
    let mut kmer = 0_u16;
    let mut rev_compl_kmer = 0_u16;
    let mut filled_bases = 0_u32; // Tracks how many consecutive valid bases we have processed

    std::iter::from_fn(move || {
        for char in seq_iter.by_ref() {
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
pub fn seq_to_unique_canon_8mers(sequence: &[u8]) -> Vec<u16> {
    // u16::MAX as usize / 32 + 1 == 2048, which is the number of u32 needed to represent all possible u16 values as bits
    let mut bitarr = BitArray::<[u32; u16::MAX as usize / 32 + 1], Msb0>::ZERO;
    seq_to_canon_8mer_iter(sequence).for_each(|canonical_kmer| {
        bitarr.set(canonical_kmer as usize, true);
    });
    bitarr.iter_ones().map(|idx| idx as u16).collect_vec()
}

pub fn seq_to_8mers(sequence: &[u8]) -> Vec<u16> {
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

pub fn get_results(results: &[lineage::EvaluationResult<'_, '_>]) -> String {
    results
        .iter()
        .map(lineage::EvaluationResult::get_output_string)
        .collect_vec()
        .join("\n")
}

/// Compute the reverse complement of a given k-mer using 32-bit parallel bit-hacks.
/// Expects a 2-bit encoded k-mer where (A=00, C=01, G=10, T=11).
/// adapted from https://github.com/gi-bielefeld/MinEncCanKmer/blob/main/canonical.c
pub fn reverse_complement(kmer: u32, k: u8) -> u32 {
    let mut value = kmer;

    // Reverse bit order
    // Swap adjacent 2-bit pairs
    value = ((value & 0xCCCCCCCC) >> 2) | ((value & 0x33333333) << 2);
    // Swap 4-bit nibbles
    value = ((value & 0xF0F0F0F0) >> 4) | ((value & 0x0F0F0F0F) << 4);
    // Swap 8-bit bytes
    value = ((value & 0xFF00FF00) >> 8) | ((value & 0x00FF00FF) << 8);
    // Swap 16-bit halves
    value = ((value & 0xFFFF0000) >> 16) | ((value & 0x0000FFFF) << 16);

    // Complement the bases (A <-> T, C <-> G) and shift down.
    let bitwidth = 32_u8;
    // Shift right to discard the unneeded high-order padding bits,
    // ensuring the resulting k-mer is correctly aligned to the right.
    value = (!value) >> (bitwidth - 2 * k);

    value
}

/// Computes the canonical form of a given k-mer by comparing it to its reverse complement and returning the lexicographically smaller one.
pub fn canonicalize(kmer: u32, k: u8) -> u32 {
    let rev_compl = reverse_complement(kmer, k);
    std::cmp::min(kmer, rev_compl)
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
        cosine_similarity, decompress_sequence, encode, euclidean_distance_l1, euclidean_norm,
        map_four_to_two_bit_repr, reverse_complement, seq_to_8mers, seq_to_unique_canon_8mers,
        KMerEncodingData,
    };

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
            let output = seq_to_unique_canon_8mers(input_seq);
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
            let output = seq_to_8mers(input_seq);
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

    #[test]
    fn test_reverse_complement() {
        let mut k = 5;
        assert_eq!(reverse_complement(0b11_1111_1111, k), 0);
        assert_eq!(reverse_complement(0b10_1101_1110, k), 0b01_0010_0001);
        assert_eq!(reverse_complement(0b01_1101_0011, k), 0b00_1110_0010);
        assert_eq!(reverse_complement(0b01_1110_0010, k), 0b01_1101_0010);
        assert_eq!(reverse_complement(0, k), 0b11_1111_1111);

        k = 3;
        assert_eq!(reverse_complement(0b11_11_11, k), 0b00_00_00);
        assert_eq!(reverse_complement(0b00_10_01, k), 0b10_01_11);
        assert_eq!(reverse_complement(0b01_00_11, k), 0b00_11_10);
        assert_eq!(reverse_complement(0, k), 0b11_11_11);

        k = 16;
        assert_eq!(
            reverse_complement(0b11111111_11111111_11111111_11111111, k),
            0
        );
        assert_eq!(
            reverse_complement(0, k),
            0b11111111_11111111_11111111_11111111
        );
        // Alternating C and G: CGCGCGCGCGCGCGCG = 0b01100110_01100110_01100110_01100110
        // Is palindrome
        let alt_cg = 0x66666666;
        assert_eq!(reverse_complement(alt_cg, k), alt_cg);
        // AAAA_CCCC_GGGG_TTTT is also palindrom
        let mixed_kmer = 0b00000000_01010101_10101010_11111111;
        assert_eq!(reverse_complement(mixed_kmer, k), mixed_kmer);

        k = 1;
        assert_eq!(reverse_complement(0b00, k), 0b11); // A -> T
        assert_eq!(reverse_complement(0b01, k), 0b10); // C -> G
        assert_eq!(reverse_complement(0b10, k), 0b01); // G -> C
        assert_eq!(reverse_complement(0b11, k), 0b00); // T -> A
    }

    #[test]
    fn test_encode_is_bijective() {
        // has also been run for all k <= 16,
        // but takes too long to run every time (k=16 ~ 10 minutes)
        for k in 1..9 {
            let encoding_data = KMerEncodingData::new(k).unwrap();

            let mut counts = vec![0u8; encoding_data.n_unique_codes as usize];

            let max_kmer = (4u64.pow(k) - 1) as u32;

            for kmer in 0..=max_kmer {
                let rev_compl_kmer = reverse_complement(kmer, k as u8);
                let code = encode(kmer, rev_compl_kmer, &encoding_data);

                // Verify that the forward and reverse complement yield identical codes
                assert_eq!(
                    code,
                    encode(rev_compl_kmer, kmer, &encoding_data),
                    "Asymmetry found: k-mer {:b} and its RC do not produce the same code for k={}",
                    kmer,
                    k
                );

                // Verify that the code stays within the strict canonical bounds
                assert!(
                    code < encoding_data.n_unique_codes,
                    "Code {} out of bounds for k={}",
                    code,
                    k
                );

                counts[code as usize] += 1;
            }

            for (code, &count) in counts.iter().enumerate() {
                let expect = if k % 2 == 0 && code < 4usize.pow(k / 2) {
                    1 // palindrome
                } else {
                    2 // kmer and its reverse complement
                };
                assert_eq!(
                    count, expect,
                    "Slot {} was hit {} times instead of {} for k={}",
                    code, count, expect, k
                );
            }
        }
    }
}
