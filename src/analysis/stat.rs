use ndarray::prelude::*;
use num::pow::pow;
use std::collections::BTreeMap;
use std::collections::HashMap;

use crate::ds::sequence::Sequence;
use crate::ds::tile::Tile;

/// Calculate the average probability that a predefined sequence is found within a random sequence with specified length and gc-content
///
/// Some bacteria are able to fight viruses using specialized enzymes known as 'restriction enzymes'.
/// Those enzymes are programmed with recognition sequences (ie. short intervals that can bind to viral DNA),
/// cutting the DNA strand somewhere in the middle of the interval.
/// Locations of the genome at which a restriction enzyme cuts and modifies DNA are denoted as 'restriction sites'.
/// The short length of the bacterial recognition sequences guarantees that a 'restriction site' can be found in a genome at random.
///
/// # Arguments
/// * `seq` - sequence to be analyzed
/// * `n` - length of the random genetic sequence
/// * `gc_content` - slice containing gc contents of random genomes
///
/// # Example
/// ```
/// use biotech::ds::sequence::Sequence;
/// use biotech::analysis::stat::expected_restriction_sites;
///
/// let recognition_seq = Sequence::from("AG");
/// let n = 10;
/// let gc = [0.25, 0.5, 0.75];
///
/// assert_eq!(expected_restriction_sites(&recognition_seq, n, &gc), [0.422, 0.563, 0.422])
/// ```
pub fn expected_restriction_sites(seq: &Sequence, n: usize, gc_content: &[f64]) -> Vec<f64> {
    // Count number of times A/T or G/C occur in sequence `seq`
    let at_count = seq.into_iter().filter(|&&n| n == b'A' || n == b'T').count();
    let gc_count = seq.into_iter().filter(|&&n| n == b'G' || n == b'C').count();

    // Calculate probability that `seq` occurs in sequence with `gc_content` and length `n`
    let mut expected_occurences = Vec::<f64>::new();
    for gc in gc_content {
        let prop_of_seq = pow((1.0 - gc) / 2.0, at_count) * pow(gc / 2.0, gc_count);
        let num_occurences = n as f64 - seq.len() as f64 + 1.0;
        let e_o = ((prop_of_seq * num_occurences) * 1000.0).ceil() / 1000.0;
        expected_occurences.push(e_o);
    }

    expected_occurences
}

/// Calculate the `NXX`-statistic assembly quality metric
///
/// The process of genome sequencing produces a large number of genetic strings called *reads*.
/// Overlapping reads are combined together forming *contigs*, which are larger constructs of varying length.
/// Contigs are used for reconstruction of the chromosome, which becomes easier if contigs are long.
///
/// The `NXX`-statistic quantifies what percentage of the assembled genome is made up of long contigs.
/// It is defined as the maximum positive integer `L`, such that the total number of nucleotides
/// of all contigs having length `≥L` is at least `XX%` of the sum of contig lengths.
///
/// # Arguments
/// * `seq` - sequence to be analyzed
/// * `n` - length of the random genetic sequence
/// * `gc_content` - slice containing gc contents of random genomes
///
/// # Example
/// ```
/// use biotech::analysis::stat::n_statistic;
/// use biotech::ds::sequence::Sequence;
/// use biotech::ds::tile::Tile;
///
/// let mut contigs = Tile::new();
/// contigs.push(Sequence::from("GATTACA"));
/// contigs.push(Sequence::from("TACTACTAC"));
/// contigs.push(Sequence::from("ATTGAT"));
/// contigs.push(Sequence::from("GAAGA"));
///
/// // Calculate N75 statistic
/// assert_eq!(n_statistic(&contigs, 75), 6);
/// ```
pub fn n_statistic(tile: &Tile, xx: usize) -> usize {
    // Map length to the number of occurences
    let mut length_map = BTreeMap::<usize, usize>::new();
    let mut total_length = 0_usize;

    for seq in tile {
        *length_map.entry(seq.len()).or_insert(0) += 1;
        total_length += seq.len();
    }

    let mut accumulator = 0_usize;
    let mut nxx = 0_usize;
    for (len, cnt) in length_map.into_iter().rev() {
        nxx = len;
        accumulator += len * cnt;
        if accumulator * 100 / total_length >= xx {
            break;
        }
    }
    nxx
}

pub fn infer_number_rna(protein: &Sequence) -> u128 {
    let codon_combs: HashMap<u8, u128> = HashMap::from([
        (b'F', 2),
        (b'I', 3),
        (b'V', 4),
        (b'L', 6),
        (b'S', 6),
        (b'P', 4),
        (b'M', 1),
        (b'T', 4),
        (b'A', 4),
        (b'Y', 2),
        (b'-', 3),
        (b'H', 2),
        (b'N', 2),
        (b'D', 2),
        (b'Q', 2),
        (b'K', 2),
        (b'E', 2),
        (b'C', 2),
        (b'G', 4),
        (b'R', 6),
        (b'W', 1),
    ]);

    // Initialize with 3 as for number of STOP codons
    let mut rna_combinations: u128 = 3;
    // Compute number of combinations
    for amino in protein {
        rna_combinations = (rna_combinations * codon_combs.get(amino).unwrap()) % 1000000;
    }
    rna_combinations
}

pub fn permutations<T: Clone>(n: usize, a: &mut Vec<T>, result: &mut Vec<Vec<T>>) {
    if n == 1 {
        result.push(a.clone());
    } else {
        for i in 0..n - 1 {
            permutations(n - 1, a, result);

            if n % 2 == 0 {
                a.swap(i, n - 1);
            } else {
                a.swap(0, n - 1);
            }
        }
        permutations(n - 1, a, result);
    }
}

pub fn signed_permutations(l: usize, a: &mut Vec<i32>, result: &mut Vec<Vec<i32>>) {
    if a.len() == l {
        result.push(a.clone());
    } else {
        for i in l..a.len() {
            // Swapping done
            a.swap(i, l);
            // Recursion called
            a[l] = -a[l];
            signed_permutations(l + 1, a, result);

            a[l] = -a[l];
            signed_permutations(l + 1, a, result);
            //backtrack
            a.swap(i, l);
        }
    }
}

pub fn p_distance_matrix(matrix: &Tile) -> Array2<f32> {
    let (rows, columns) = matrix.size();
    let mut distances = Array2::<f32>::zeros((rows, rows));
    let mut p_dist: f32 = 0.0;
    for (i, a_row) in matrix.into_iter().enumerate() {
        for (j, b_row) in matrix.into_iter().enumerate() {
            p_dist = 0.0;
            if i != j {
                p_dist = a_row.into_iter().zip(b_row).filter(|(a, b)| a != b).count() as f32;
            }
            distances[(i, j)] = p_dist / (columns as f32);
            distances[(j, i)] = p_dist / (columns as f32);
        }
    }
    distances
}

pub fn calc_profile(arr: &Array3<u16>) -> Array2<u16> {
    // Squash tensor into 2d array
    // Sum the occurs of letters different letters
    // Transpose at the end to get expected shape
    arr.sum_axis(Axis(0)).reversed_axes()
}

pub fn calc_consensus(arr: &Array2<u16>) -> Array1<u8> {
    // Allocate container for result
    let mut consensus = Array1::<u8>::zeros(arr.len_of(Axis(1)));
    // Calculate maximum index for every dimension in array
    for (i, ax) in arr.axis_iter(Axis(1)).enumerate() {
        let mut max: u16 = 0;
        let mut index: usize = 0;
        for (j, it) in ax.indexed_iter() {
            if *it > max {
                max = *it;
                index = j;
            }
        }
        consensus[i] = index as u8;
    }
    consensus
}

pub fn random_substrings(seq: &Sequence, gc_content: &[f64]) -> Vec<f64> {
    let mut probabilities: Vec<f64> = vec![];
    let mut prob_map: HashMap<u8, f64> =
        HashMap::from([(b'A', 1.0), (b'C', 1.0), (b'T', 1.0), (b'G', 1.0)]);

    for x in gc_content {
        *prob_map.get_mut(&(b'G')).unwrap() = x / 2.0;
        *prob_map.get_mut(&(b'C')).unwrap() = x / 2.0;
        *prob_map.get_mut(&(b'T')).unwrap() = (1.0 - x) / 2.0;
        *prob_map.get_mut(&(b'A')).unwrap() = (1.0 - x) / 2.0;
        // calculate probability of given string
        // for provided GC content
        let mut prop = 1.0;
        for c in seq {
            prop *= prob_map[c]
        }
        // probabilities.push(prop.log10());
        probabilities.push(prop);
    }
    probabilities
}

pub fn random_motifs(n: usize, gc_content: f64, seq: &Sequence) -> f64 {
    let prob = random_substrings(seq, &[gc_content])
        .first()
        .cloned()
        .unwrap();
    1.0 - pow(1.0 - prob, n)
}
