use std::io;
use ndarray::prelude::*;
use super::DNA;

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

pub fn hamming_distance(s1: &String, s2: &String) -> Result<u32, io::Error> {
    let both = s1.chars().zip(s2.chars());
    let mut hamming : u32 = 0;
    for pair in both {
        if pair.0 != pair.1 {
            hamming += 1;
        }
    }
    Ok(hamming)
}

pub fn knuth_morris_pratt(st: &[u8], pat: &[u8]) -> Vec<usize> {
   
    // Build the partial match table
    let mut partial = vec![0];
    for i in 1..pat.len() {
        let mut j = partial[i - 1];
        while j > 0 && pat[j] != pat[i] {
            j = partial[j - 1];
        }
        partial.push(if pat[j] == pat[i] { j + 1 } else { j });
    }

    // Read 'string' to find 'pattern'
    let mut ret = vec![];
    let mut j = 0;

    for (i, &c) in st.iter().enumerate() {
        while j > 0 && c != pat[j] {
            j = partial[j - 1];
        }
        if c == pat[j] {
            j += 1;
        }
        if j == pat.len() {
            ret.push(i + 1 - j);
            j = partial[j - 1];
        }
    }
    ret
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


pub fn overlap_graph<T: Iterator<Item=DNA>>(tile: T) {

    for dna in tile {

        print!("{}\n", dna);
    }


}

// N{P}[ST]{P}
// pub fn generate_motifs(target: &Vec<u8>, result: &mut Vec<String>, i: usize, temp: &mut Vec<u8>) {

//     if i == target.len() {
//         result.push(std::str::from_utf8(temp).unwrap().to_string());
//         return;
//     }

//     let alphabet = String::from("FLSYCWPHQRITMNKVADEG").into_bytes();
//     let mut alternatives = Vec::<u8>::new();
//     let mut j = i;
//     if target[j] == b'[' {
//         j += 1;
//         while j < target.len() && target[j] != b']' {
//             alternatives.push(target[j]);
//             j += 1;
//         }
//         for elem in alternatives {
//             temp.push(elem);
//             generate_motifs(target, result, j + 1, temp);
//             temp.pop();
//         }
//     }
//     else if target[j] == b'{' {
//         j += 1;
//         for elem in alphabet {
//             if elem != target[j] {
//                 temp.push(elem);
//                 generate_motifs(target, result, j + 2, temp);
//                 temp.pop();
//             }
//         }
//         j += 2;
//     }
//     else{
//         temp.push(target[i]);
//         generate_motifs(target, result, i + 1, temp);
//         temp.pop();
//     }
// }

// fn search_motifs(st: String, pat: String) -> Vec<usize> {

//     let st = st.into_bytes();
//     let pat = pat.into_bytes();
//     let mut res = Vec::<String>::new();
//     let mut temp = Vec::<u8>::new();

//     generate_motifs(&pat, &mut res, 0, &mut temp);

//     let mut pos = Vec::<usize>::new();
//     for elem in &res {
//         pos.append(&mut knuth_morris_pratt(&st, elem.as_bytes()));    
//     }
//     pos
// }