use statrs::function::factorial::{binomial, factorial};
use std::collections::HashMap;
use num::*;
use num::ToPrimitive;
use num::pow::pow;
use num::traits::Pow;
use ndarray::prelude::*;
use std::cmp;
use itertools::iproduct;

use crate::ds::sequence::Sequence;
use crate::ds::tile::Tile;

use super::seq::count_nucleotides;

pub fn mendel_first_law(k: u16, m: u16, n: u16 ) -> f32 {
    
    let k_f32 = k as f32;
    let m_f32 = m as f32;
    let n_f32 = n as f32;
    let total : f32 = k_f32 + m_f32 + n_f32;

    // calculate the probability for a homozygous recessive trait
    // heterozygous organisms mating with each other
    let m_p : f32 = m_f32/total * (m_f32 - 1.0 )/(total - 1.0) * 0.25;
    // homozygous organisms mating with each other
    let n_p : f32 = n_f32/total * (n_f32 - 1.0 )/(total - 1.0);
    // heterozygous with homozygous mating
    let mut mn_p : f32 = m_f32/total * n_f32/(total - 1.0) * 0.5;
    mn_p += m_f32/(total - 1.0) * n_f32/total * 0.5;
    
    // probability of a dominant allele is 1 - prob_homozygous
    1.0 - (m_p + n_p + mn_p)
}

pub fn fibo_mortal(n: usize, m: usize) -> usize {
    
    let mut memory: Vec<usize> = vec![1; n];
    let mut val: usize = 0;

    for i in 2..n {    
        val = memory[i - 1] + memory[i - 2];
        if i == m {
            val = val - 1;
        }
        if i > m {
            val = val - memory[i - m - 1];
        }
        memory[i] = val;
    }
    *memory.last().unwrap()
}

pub fn fibo(n: usize, k: u64) -> u64 {
    // F1 == F2 == 1
    // Initialize the vector with 1's 
    let mut memory: Vec<u64> = vec![1; n];
    // Members of (i-1) generation are too young too breed
    // Members of (i-2) generation can breed k new members
    for i in 2..n {
        memory[i] = memory[i-1] + memory[i-2]*k;
    }
    memory[n -1]
}

pub fn expected_offspring(x: u16, y: u16, z: u16, q: u16, p: u16, r: u16 ) -> f32 {
    // The probability of AA-AA, AA-Aa, AA-aa couples
    // having dominant phenotype offspring is 1
    let xyz_e : f32 = 2.0 * (x as f32 + y as f32 + z as f32);
    // P[dominant offspring] = 0.75
    let q_e : f32 = 2.0 * q as f32 * 0.75;
    // P[dominant offspring] = 0.5
    let p_e : f32 = 2.0 * p as f32 *  0.5;
    xyz_e + q_e + p_e
}

pub fn partial_permutations(n: u64, k: u64) -> u64 {
    let combinations = binomial(n, k);
    let permutations = factorial(k);
    (combinations * permutations) as u64 % 1000000
}

pub fn permutations<T: Clone>(n : usize, a : &mut Vec<T>, result : &mut Vec<Vec<T>>) {
    if n == 1 {
        result.push(a.clone());
    }
    else {
        for i in  0 .. n - 1 {
            permutations(n - 1, a, result);

            if n % 2 == 0 {
                a.swap(i, n - 1);
            }
            else {
                a.swap(0, n - 1);
            }
        }
        permutations(n - 1, a, result);
    }
}

pub fn calc_profile(arr: &Array3<u16>) -> Array2::<u16>  {
    // Squash tensor into 2d array 
    // Sum the occurs of letters different letters
    // Transpose at the end to get expected shape
    arr.sum_axis(Axis(0)).reversed_axes()
}

pub fn calc_consensus(arr: &Array2<u16>) -> Array1<u8> {
    // Allocate container for result
    let mut consensus = Array1::<u8>::zeros((arr.len_of(Axis(1))));
    // Calculate maximum index for every dimension in array
    for (i, ax) in arr.axis_iter(Axis(1)).enumerate() {
        let mut max : u16 = 0;
        let mut index : usize = 0;
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

pub fn p_distance_matrix(matrix: &Tile) -> Array2<f32> {

    let (rows, columns) = matrix.size();
    let mut distances =  Array2::<f32>::zeros((rows, rows));
    let mut p_dist : f32 = 0.0;
    
    for (i, a_row) in matrix.into_iter().enumerate() {
        for (j, b_row) in matrix.into_iter().enumerate() {
            p_dist = 0.0;
            if i != j {
                p_dist = a_row.into_iter().zip(b_row).filter(|(a,b)| a != b).count() as f32;
            }
            distances[(i,j)] = p_dist / (columns as f32);
            distances[(j,i)] = p_dist / (columns as f32);
        }
    }
    distances
}

pub fn random_substrings(seq: &Sequence, gc_content: &[f64]) -> Vec<f64> {

    let mut probabilities : Vec<f64> = vec![];
    let mut prob_map : HashMap<u8, f64> = HashMap::from([ 
        (b'A', 1.0), (b'C', 1.0), (b'T', 1.0), (b'G', 1.0)
    ]);

    for x in gc_content {
        *prob_map.get_mut(&(b'G')).unwrap() = x / 2.0;
        *prob_map.get_mut(&(b'C')).unwrap() = x / 2.0;
        *prob_map.get_mut(&(b'T')).unwrap() = (1.0-x) / 2.0;
        *prob_map.get_mut(&(b'A')).unwrap() = (1.0-x) / 2.0;
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

    let prob = random_substrings(seq, &[gc_content]).first().cloned().unwrap();
    1.0- pow((1.0 - prob), n)
}


pub fn number_subsets(n: u64) -> u64 {

    let base = BigUint::from(2 as u64);
    let modulus = BigUint::from(1000000 as u64);
    let subsets = base.pow(n) % modulus;
    ToPrimitive::to_u64(&subsets).unwrap()
}

pub fn signed_permutations(l : usize, a : &mut Vec<i32>, result : &mut Vec<Vec<i32>>) {

    if a.len() == l {
        result.push(a.clone());
    }
    else
    {
        for i in  l..a.len() {
            // Swapping done
            a.swap(i, l);
 
            // Recursion called
            a[l] = -a[l];
            signed_permutations(l+1, a, result);

            a[l] = -a[l];
            signed_permutations(l+1, a, result);
            
            //backtrack
            a.swap(i, l);
        }
    }
}

// Calculate the expected probability of occurrences of sequence `seq`
// within a_i sequence of length n, and gc_content_i
pub fn expected_restriction_sites(seq: &Sequence, n: usize, gc_content: &[f32]) -> Vec<f32> {
    
    // Count number of times A/T
    // or G/C occur in sequence `seq`
    let mut at_count = seq.into_iter()
                          .filter(|n| **n == b'A' || **n == b'T')
                          .count();
    
    let mut gc_count = seq.into_iter()
                          .filter(|n| **n == b'G' || **n == b'C')
                          .count();

    // Calculate probability that `seq` occurs 
    // in sequence with hypothetic `gc_content` and length `n`
    let mut expected_occurences  = Vec::<f32>::new();
    for gc in gc_content {
        let prop_of_seq = pow((1.0 - gc) / 2.0, at_count) * pow(gc / 2.0, gc_count);
        let num_occurences = n as f32 - seq.len() as f32 + 1.0;
        let e_o = ((prop_of_seq * num_occurences) * 1000.0).ceil() / 1000.0;
        expected_occurences.push(e_o);
    }

    expected_occurences
}

pub fn count_basepair_matchings(rna: &Sequence) -> BigUint {

    let cnt = count_nucleotides(rna);

    let gc_min = cmp::min(cnt[&b'G'], cnt[&b'C']);
    let gc_max = cmp::max(cnt[&b'G'], cnt[&b'C']);

    let au_min = cmp::min(cnt[&b'A'], cnt[&b'U']);
    let au_max = cmp::max(cnt[&b'A'], cnt[&b'U']);

    // Number of possible GC/AU matchings is equivalent to:
    // (gc_cnt)(gc_cnt-1)(gc_cnt-2)...(2)(1) = gc_cnt!
    fn factorial(num: u128) -> u128 {
        match num {
            0  => 1,
            1.. => (1..num+1).product(),
        }
    }

    // Generate the number of possible matchings for `GC` and `AU` respectively
    let x1 = factorial(gc_min as u128) * binomial(gc_max as u64, gc_min as u64) as u128;
    let x2 = factorial(au_min as u128) * binomial(au_max as u64, au_min as u64) as u128;

    // Combine the partial results into total number of matchings
    let x1_long = BigUint::from(x1);
    let x2_long = BigUint::from(x2);
    
    x1_long * x2_long
}

/// Return frequency and value of the most frequently occurring shift between peaks of two mass spectrums.
///   
/// Mass spectrography is used to infer protein from their weight.
/// Comparing spectral representations can give useful insight into the similarity of corresponding particles.
/// A given protein will often occur nested in a more complex structure, which will shift the peaks of the mass spectrum.
/// In order to still be able to quantify the similarity of the spectra, the shift must be calculated.
/// 
/// # Arguments
///
/// * `s1` - evaluated mass spectrum
/// * `s2` - evaluated mass spectrum
/// 
pub fn spectral_mass_shift(s1: Vec<f32>, s2: Vec<f32>) -> (usize, f32) {

    // Calculate minkovsky difference for s1 x s2
    let minkovsky_diff : Vec<i32> = iproduct!(s1.iter(), s2.iter())
            .map(|(v, x)| ((v - x) * 1000.0) as i32)
            .collect();

    // Calculate mode of obtained set equvalent to mass shift
    let mut counts = HashMap::<i32, usize>::new();
    let max = minkovsky_diff.iter().copied().max_by_key(|&n| {
        let count = counts.entry(n).or_insert(0);
        *count += 1;
        *count
    }).unwrap_or(0);

    (*counts.get(&max).unwrap(), max as f32 / 1000.0)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mendel_first_law() {
        let k: u16 = 2;
        let m: u16 = 2;
        let n: u16 = 2;
        assert_eq!(0.7833333, 
                    mendel_first_law(k, m, n));
    }

    #[test]
    fn test_fibo() {
        let n: usize = 5;
        let k: u64 = 3;
        assert_eq!(19, fibo(n, k));
    }

    #[test]
    fn test_fibo_mortal() {
        let n: usize = 6;
        let k: usize = 3;
        assert_eq!(4, fibo_mortal(n, k));
    }

    #[test]
    fn test_expected_offspring() {
        let x: u16 = 18137;
        let y: u16 = 16426;
        let z: u16 = 18904;
        let q: u16 = 18674;
        let p: u16 = 18160;
        let r: u16 = 18728;
        assert_eq!(153105.0, expected_offspring(x, y, z, q, p, r));
    }

    #[test]
    fn test_partial_permutations() {
        let part = partial_permutations(21, 7);
        assert_eq!(51200, part);
    }

    #[test]
    fn test_random_substring() {
        let seq = Sequence::from("GACGGACAAGGGCCCCCGTGTATTGTACTTGGGCCCATGTGCC\
                                    CGACCTCGGTAAGTCCATCAGGAGTGCACGAGGACCACCATTTCAAGAAA");
        let arr =  [0.110, 0.127, 0.183, 0.256];
        assert_eq!(vec![ -80.82637701756781, -77.85362263786175, -70.59699957974087, -64.49615401338707], 
                        random_substrings(&seq, &arr).iter().map(|x| x.log10()).collect::<Vec<_>>());
    }

    #[test]
    fn test_random_motifs() {
        assert_eq!(0.3131377954665139, random_motifs(97232, 0.506715, &Sequence::from("GCGTAAGAC")))
    }
}