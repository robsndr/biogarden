use statrs::function::factorial::{binomial, factorial};
use std::collections::HashMap;
use num::*;
use num::ToPrimitive;
use num::pow::pow;
use num::traits::Pow;
use ndarray::prelude::*;

use crate::ds::sequence::Sequence;
use crate::ds::tile::Tile;

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
        ('A' as u8, 1.0), ('C' as u8, 1.0), 
        ('T' as u8, 1.0), ('G' as u8, 1.0)
    ]);

    for x in gc_content {
        *prob_map.get_mut(&('G' as u8)).unwrap() = x / 2.0;
        *prob_map.get_mut(&('C' as u8)).unwrap() = x / 2.0;
        *prob_map.get_mut(&('T' as u8)).unwrap() = (1.0-x) / 2.0;
        *prob_map.get_mut(&('A' as u8)).unwrap() = (1.0-x) / 2.0;
        // calculate probability of given string
        // for provided GC content
        let mut prop = 1.0;
        for c in seq {
            prop *= prob_map[c]
        }
        probabilities.push(prop.log10());
    }
    return probabilities;
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
        {
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
        assert_eq!(vec![ -80.82637701756781, -77.85362263786175, 
                            -70.59699957974087, -64.49615401338707], random_substrings(&seq, &arr));
    }
}