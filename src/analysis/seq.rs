use std::cmp;
use std::collections::HashMap;
use std::collections::HashSet;

use crate::ds::sequence::Sequence;
use crate::error::{BioError, Result};
use crate::ds::builders::suffix_tree::SuffixTreeBuilder;

/// Gets a HashMap containing occurrences of nucleobases in a bio-sequence
///
/// # Arguments
/// * `seq` - sequence of nucleobases
///
/// # Example
/// ```
/// use biotech::analysis::seq::count_nucleotides;
/// use biotech::ds::sequence::Sequence;
/// use std::collections::HashMap;
///
/// let dna = Sequence::from("AGCTTTTCATTCT");
/// let counts = HashMap::<u8, usize>::from([(b'A', 2), (b'C', 3), (b'T', 7), (b'G', 1)]);
/// assert_eq!(count_nucleotides(&dna), counts);
/// ```
pub fn count_nucleotides(seq: &Sequence) -> HashMap<u8, usize> {
    let mut count = HashMap::<u8, usize>::new();
    for c in seq {
        *count.entry(*c).or_insert(0) += 1;
    }
    count
}

/// Calculates percentage of guanine (G) or cytosine (C) nucleotides in a sequence
///
/// Due to base pairing relations, guanine and cytosine will occur with the same frequency.
/// Members of different species might be distinguished by variations in GC-content.
///
/// # Arguments
/// * `seq` - sequence of nucleobases
///
/// # Example
/// ```
/// use biotech::analysis::seq::gc_content;
/// use biotech::ds::sequence::Sequence;
///
/// let dna = Sequence::from("CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTT");
/// assert_eq!(gc_content(&dna), 0.5263157894736842)
/// ```
pub fn gc_content(seq: &Sequence) -> f64 {
    let gc_count = seq
        .into_iter()
        .filter(|x| **x == b'G' || **x == b'C')
        .count();
    gc_count as f64 / seq.len() as f64
}

/// Calculates the number of point mutations, by which two sequences differ.
///
/// The point mutation replaces one base with its complement at a given location.
/// It is the most common form of mutation that occurs in an evolutionary path of two strands.
///
/// # Arguments
/// * `seq1`, `seq2` - sequences to for which the hamming distance is calculated
///
/// # Example
/// ```
/// use biotech::analysis::seq::hamming_distance;
/// use biotech::ds::sequence::Sequence;
///
/// let a = Sequence::from("GAGCCTACTAACGGGAT");
/// let b = Sequence::from("CATCGTAATGACGGCCT");
/// let h_dist = hamming_distance(&a, &b);
/// assert_eq!(h_dist.unwrap(), 7);
/// ```
pub fn hamming_distance(seq1: &Sequence, seq2: &Sequence) -> Result<usize> {
    match seq1.len() == seq2.len() {
        true => Ok(seq1
            .into_iter()
            .zip(seq2.into_iter())
            .filter(|(a, b)| a != b)
            .count()),
        false => Err(BioError::InvalidInputSize),
    }
}

/// Obtain edit-distance between two genetic sequences
///   
/// The Hamming distance provides a way to model point mutations transforming one genetic string into another.
/// In practice, point mutations include insertions and deletions in addition to replacements only.
/// This can produce genetic sequence that vary in length, and cannot be compared using hamming distance.
/// In such scenarios, a measure of the minimum number of replacements / insertions / deletions is provided by edit distance.
///
/// # Arguments
/// * `seq1`, `seq2` - sequences to between which the edit distance is calculated
///
/// # Example
/// ```
/// use biotech::analysis::seq::edit_distance;
/// use biotech::ds::sequence::Sequence;
///
/// let a = Sequence::from("ACTGGATTC");
/// let b = Sequence::from("ACGT");
/// let ed = edit_distance(&a, &b);
/// assert_eq!(ed.unwrap(), 5);
/// ```
pub fn edit_distance(seq1: &Sequence, seq2: &Sequence) -> Result<usize> {
    // Data containers
    let mut memo = vec![vec![0_u128; seq2.len() + 1]; seq1.len() + 1];

    // initialize table
    for i in 0..(seq1.len() + 1) {
        memo[i][0] = i as u128;
    }
    for j in 0..(seq2.len() + 1) {
        memo[0][j] = j as u128;
    }

    // Calculate edit-distance dp table
    for i in 1..seq1.len() + 1 {
        for j in 1..seq2.len() + 1 {
            let minimum = cmp::min(
                memo[i - 1][j - 1] + ((seq1[i - 1] != seq2[j - 1]) as u128),
                cmp::min(memo[i][j - 1] + 1, memo[i - 1][j] + 1),
            );
            // Update edit distance in memoization table
            memo[i][j] = minimum;
        }
    }

    Ok(memo[seq1.len()][seq2.len()] as usize)
}

/// Calculate the ratio of transitions to transversions between two genetic strings
///   
/// Point mutations between genetic strings occur in two flavors: transitions or transversions.
/// Transitions occur for replacements between thymine and cytosine (A <-> C),
/// or adenine and guanine (uracil (U) for RNA) (A <-> G).
/// The number of transitions to transversions is usually a bit higher for coding regions.
/// This makes the transition-transversion ratio a suitable tool to identify those in a long strand.
///
/// # Arguments
/// * `seq1`, `seq2` - sequences to between which transition-transversion ratio is calculated
///
/// # Example
/// ```
/// use biotech::analysis::seq::transition_transversion_ratio;
/// use biotech::ds::sequence::Sequence;
///
/// let a = Sequence::from("TTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT");
/// let b = Sequence::from("TCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT");
/// assert_eq!(transition_transversion_ratio(&a, &b).unwrap(), 1.875);
/// ```
pub fn transition_transversion_ratio(seq1: &Sequence, seq2: &Sequence) -> Result<f64> {
    if seq1.len() != seq2.len() {
        return Err(BioError::InvalidInputSize);
    }

    let missmatch_iter = seq1
        .into_iter()
        .zip(seq2.into_iter())
        .filter(|(a, b)| a != b);
    let hamming_distance = missmatch_iter.clone().count() as f64;
    let mut transitions: f64 = 0.0;
    for (&a, &b) in missmatch_iter {
        transitions += match (a, b) {
            (b'C', b'T') => 1.0,
            (b'T', b'C') => 1.0,
            (b'A', b'G') => 1.0,
            (b'G', b'A') => 1.0,
            _ => 0.0,
        }
    }
    // transversions == hamming_distance - transitions
    Ok(transitions / (hamming_distance - transitions))
}

/// Returns the linguistic complexity of a genetic sequence
///   
/// The human genome contains a large number of repeats, with the ALU repeat being one of the most frequent.
/// ALU repeats can create copies of themselves and reinsert those into the genome.
/// This makes them occur up to a million times within the human genetic sequence.
///
/// The linguistic complexity provides a useful measure that helps quantify how repetitive a genome is.
/// It is defined as the percentage of observed substrings of `s` to the total number that are theoretically possible.
///
/// # Arguments
/// * `seq` - sequence for which linguistic complexity should be calculated
///
/// # Example
/// ```
/// use biotech::analysis::seq::linguistic_complexity;
/// use biotech::ds::sequence::Sequence;
///
/// let a = Sequence::from("ATTTGGATT");
/// assert_eq!(linguistic_complexity(&a).unwrap(), 0.875);
/// ```
pub fn linguistic_complexity(seq: &Sequence) -> Result<f32> {
    // Define alphabet
    let alphabet = HashSet::<u8>::from([b'A', b'C', b'T', b'G']);
    let alphabet_len = alphabet.len();

    // Add the strings to tree and traverse from root to node.
    // Each root to node path will denote suffixes of a string.
    let mut ukonnen_builder = SuffixTreeBuilder::new(alphabet);
    let graph = ukonnen_builder.build(seq);

    // Count unique substrings
    // All the prefixes of these suffixes are unique substrings.
    // Their number can be obtained by summing-up the length of all edges
    let mut num_substrings = 0;
    for eid in graph.edges() {
        // Get start edge
        let start = graph
            .get_edge(eid)
            .data
            .as_ref()
            .ok_or(BioError::ItemNotFound)?
            .suffix_start;
        // Get stop edge
        let stop = graph
            .get_edge(eid)
            .data
            .as_ref()
            .ok_or(BioError::ItemNotFound)?
            .suffix_stop;
        num_substrings += stop - start + 1;
    }

    // The maximum number of k-length substrings for n-letter string is either:
    // limited by the number of substrings that can be formed from a given alphabet (4^k)
    // or by the number of k-windows that can be shifted within n-length string
    let mut max_complexity = 0;
    for k in 1..seq.len() + 1 {
        // Perform 4^k only when 4^k < seq.len(), however 4^k test can result in overflow
        // Use k < log4(seq.len()) instead. Convert log4(seq.len()) => log10(seq.len())/log10(4)
        // This will increase accuracy as log10 is better:
        // Check -> https://doc.rust-lang.org/std/primitive.f64.html#method.log
        if (k as f32) < ((seq.len() as f32).log10() / (4.0_f32).log10()) {
            max_complexity += u128::pow(alphabet_len as u128, k as u32);
        } else {
            max_complexity += (seq.len() - k + 1) as u128;
        }
    }
    Ok(num_substrings as f32 / max_complexity as f32)
}
