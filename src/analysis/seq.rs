use std::cmp;
use std::collections::HashMap;

use crate::ds::sequence::Sequence;
use crate::error::{BioError, Result};

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
/// let counts = HashMap::<u8, usize>::from([(b'A', 2),(b'C', 3),
///                                          (b'T', 7),(b'G', 1)]);
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
/// assert!(h_dist.is_ok());
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

/// Obtain tuple containing edit-distance and edit-alignment of two genetic sequences.
///   
/// The Hamming distance provides a way to model point mutations transforming one genetic string into another.
/// In practice, point mutations include insertions and deletions in addition to replacements only.
/// This can produce genetic sequence that vary in length, and cannot be compared using hamming distance.
/// In such scenarios, a measure of the minimum number of replacements / insertions / deletions between two sequences, is provided by edit distance.
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
/// assert!(ed.is_ok());
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
