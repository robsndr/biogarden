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
/// * `s1`, `s2` - sequences to for which the hamming distance is calculated
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
pub fn hamming_distance(s1: &Sequence, s2: &Sequence) -> Result<usize> {
    match s1.len() == s2.len() {
        true => Ok(s1
            .into_iter()
            .zip(s2.into_iter())
            .filter(|(a, b)| a != b)
            .count()),
        false => Err(BioError::InvalidInputSize),
    }
}
