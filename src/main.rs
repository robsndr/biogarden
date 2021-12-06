use std::io;
mod dna;
use std::collections::HashMap;
use ndarray::prelude::*;

// count number of chars in dna sequence
fn count_nuclea(dna : &String) -> [u16; 24] {
    let mut count: [u16; 24] = [0; 24];
    for c in dna.chars() {
        count[(c as usize) - 65] += 1;
    }
    return count;
}

// count number of chars in dna sequence
fn transcribe_rna(dna : &String) -> String {
    dna.replace("T", "U")
}

fn mendel_first_law(k: u16, m: u16, n: u16 ) -> f32 {
    
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


fn expected_offspring(x: u16, y: u16, z: u16, q: u16, p: u16, r: u16 ) -> f32 {
    // The probability of AA-AA, AA-Aa, AA-aa couples
    // having dominant phenotype offspring is 1
    let xyz_e : f32 = 2.0 * (x as f32 + y as f32 + z as f32);
    // P[dominant offspring] = 0.75
    let q_e : f32 = 2.0 * q as f32 * 0.75;
    // P[dominant offspring] = 0.5
    let p_e : f32 = 2.0 * p as f32 *  0.5;
    xyz_e + q_e + p_e
}

fn fibo(n: usize, k: u64) -> u64 {
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

// fn fibo_die(n: usize, m: usize) -> u64 {
//     // F1 == F2 == 1
//     // Initialize the vector with 1's 
//     let mut memory: Vec<u64> = vec![1; n];
//     for i in 2..n {
//         if i < m {
//             memory[i] = memory[i-1] + memory[i-2];
//         }
//         else{
//             memory[i] = memory[i-1] + memory[i-2] - memory[i-m];
//         }
//     }
//     memory[n -1]
// }

fn gc_content(s: &String) -> f64 {
    let mut gc_count : u32 = 0;
    for c in s.chars() {
        if c == 'G' ||  c == 'C' {
            gc_count += 1;
        }
    }
    gc_count as f64 / s.chars().count() as f64
}

fn hamming_distance(s1: &String, s2: &String) -> Result<u32, io::Error> {
    let both = s1.chars().zip(s2.chars());
    let mut hamming : u32 = 0;
    for pair in both {
        if pair.0 != pair.1 {
            hamming += 1;
        }
    }
    Ok(hamming)
}

fn translate_protein(mRNA: &String) -> String {
    // mRNA <-> amino-acid translation table (codon table)
    let codon_table: HashMap<&str, &str> = HashMap::from(
        [ ("UUU", "F"),    ("CUU", "L"),   ("AUU", "I"),   ("GUU", "V"),
          ("UUC", "F"),    ("CUC", "L"),   ("AUC", "I"),   ("GUC", "V"),
          ("UUA", "L"),    ("CUA", "L"),   ("AUA", "I"),   ("GUA", "V"),
          ("UUG", "L"),    ("CUG", "L"),   ("AUG", "M"),   ("GUG", "V"),
          ("UCU", "S"),    ("CCU", "P"),   ("ACU", "T"),   ("GCU", "A"),
          ("UCC", "S"),    ("CCC", "P"),   ("ACC", "T"),   ("GCC", "A"),
          ("UCA", "S"),    ("CCA", "P"),   ("ACA", "T"),   ("GCA", "A"),
          ("UCG", "S"),    ("CCG", "P"),   ("ACG", "T"),   ("GCG", "A"),
          ("UAU", "Y"),    ("CAU", "H"),   ("AAU", "N"),   ("GAU", "D"),
          ("UAC", "Y"),    ("CAC", "H"),   ("AAC", "N"),   ("GAC", "D"),
          ("UAA", "Stop"), ("CAA", "Q"),   ("AAA", "K"),   ("GAA", "E"),
          ("UAG", "Stop"), ("CAG", "Q"),   ("AAG", "K"),   ("GAG", "E"),
          ("UGU", "C"),    ("CGU", "R"),   ("AGU", "S"),   ("GGU", "G"),
          ("UGC", "C"),    ("CGC", "R"),   ("AGC", "S"),   ("GGC", "G"),
          ("UGA", "Stop"), ("CGA", "R"),   ("AGA", "R"),   ("GGA", "G"),
          ("UGG", "W"),    ("CGG", "R"),   ("AGG", "R"),   ("GGG", "G") ]);
    // Container for final result of transcription
    let mut amino_acid: String = String::from("");
    // Run the translation 
    let mut z = mRNA.chars().peekable();
    while z.peek().is_some() {
        let chunk: String = z.by_ref().take(3).collect();
        match codon_table.get(&chunk as &str) {
            Some(value) => {
                if value == &"Stop"{
                    break;
                }
                amino_acid.push_str(value);
            },
            None => println!("Codon not found in codon table.")
        }
    }
    amino_acid
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

// N{P}[ST]{P}
pub fn generate_motifs(target: &Vec<u8>, result: &mut Vec<String>, i: usize, temp: &mut Vec<u8>) {

    if i == target.len() {
        result.push(std::str::from_utf8(temp).unwrap().to_string());
        return;
    }

    let alphabet = String::from("FLSYCWPHQRITMNKVADEG").into_bytes();
    let mut alternatives = Vec::<u8>::new();
    let mut j = i;
    if target[j] == b'[' {
        j += 1;
        while j < target.len() && target[j] != b']' {
            alternatives.push(target[j]);
            j += 1;
        }
        for elem in alternatives {
            temp.push(elem);
            generate_motifs(target, result, j + 1, temp);
            temp.pop();
        }
    }
    else if target[j] == b'{' {
        j += 1;
        for elem in alphabet {
            if elem != target[j] {
                temp.push(elem);
                generate_motifs(target, result, j + 2, temp);
                temp.pop();
            }
        }
        j += 2;
    }
    else{
        temp.push(target[i]);
        generate_motifs(target, result, i + 1, temp);
        temp.pop();
    }
}

fn search_motifs(st: String, pat: String) -> Vec<usize> {

    let st = st.into_bytes();
    let pat = pat.into_bytes();
    let mut res = Vec::<String>::new();
    let mut temp = Vec::<u8>::new();

    generate_motifs(&pat, &mut res, 0, &mut temp);

    let mut pos = Vec::<usize>::new();
    for elem in &res {
        pos.append(&mut knuth_morris_pratt(&st, elem.as_bytes()));    
    }
    pos
}

// Consensus and Profile
//
//
pub fn encode_input(arr: &Array2<u8>) -> Array3::<u16> {
    // Tensor with result encoding
    let mut converted = Array3::<u16>::zeros((arr.len_of(Axis(0)), arr.len_of(Axis(1)), 4));
    // Encode 2D input char arry into output tensor
    // A => [1,0,0,0], B => [0,1,0,0]
    // C => [0,0,1,0], D => [0,0,0,1]
    for ((i, j), value) in arr.indexed_iter() {
            match value {
                b'A' => {
                    converted[(i,j,0)] = 1;
                },
                b'C' => {
                    converted[(i,j,1)] = 1
                },
                b'G' => {
                    converted[(i,j,2)] = 1
                },
                b'T' => {
                    converted[(i,j,3)] = 1
                },
                _ => ()
            }
    }
    converted
}

pub fn calc_profile(arr: &Array3<u16>) -> Array2::<u16>  {
    // Squash tensor into 2d array 
    // Sum the occurs of letters different letters
    // Transpose at the end to get expected shape
    arr.sum_axis(Axis(0)).reversed_axes()
}


pub fn calc_consensus(arr: &Array2<u16>) -> Array1<u16> {
    // Allocate container for result
    let mut consensus = Array1::<u16>::zeros((arr.len_of(Axis(1))));
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
        consensus[i] = index as u16;
    }
    consensus
}

pub fn encode_output(arr: &Array1<u16>) -> Vec::<char> {
    // Tensor with result encoding
    let mut out : Vec<char> = vec![];

    for (i, value) in arr.indexed_iter() {
            match value {
                0 => {
                    out.push('A');
                },
                1 => {
                    out.push('C');
                },
                2 => {
                    out.push('G');
                },
                3 => {
                    out.push('T');
                },
                _ => ()
            }
    }
    out
}


// pub fn permute(arr: &, res: &Vec<Array1<u8>>) -> {
//     void perm(char* s, int n, int i){
//         if i >= n-1 {
//             print(s);
//         }
//         else {
//           perm(s, n, i+1);
//           for (int j = i+1; j<n; j++){
//             swap(s[i], s[j]);
//             perm(s, n, i+1);
//             swap(s[i], s[j]);
//           }
//         }
//       }
      
//       perm("ABC", 3, 0);
// }

use dna::DNA;

fn main() {
    
    let mut my_str1 : String = String::from("AAAACCCGGT");
    let mut my_str2 : String = String::from("N{P}[ST]{P}");

    
    let mut d : DNA = DNA::from(my_str1);
    print!("{}", d);
    d.complement();
    print!("{}", d);
    // print!("{}", mendel_first_law(15, 17, 19));
    // print!("{}", expected_offspring(18137, 16426, 18904, 18674, 18160, 18728));
    // print!("{}", fibo_die(6, 3));
    // print!("{}", gc_content(&my_str));
    // print!("{:#?}", find_repeats(&my_str1, &my_str2));
    // print!("{:#?}", knuth_morris_pratt(my_str1, my_str2));

    // match hamming_distance(&my_str1, &my_str2) {
    //     Ok(value) => {
    //         print!("{}", value);
    //     }
    //     Err(err) => {
    //         println!("Error calculating the hamming distance: {}", err);
    //     }
    // }

    // search_motifs(my_str1, my_str2);
    // println!("{:#?}", search_motifs(my_str1, my_str2));

    // for ((x, y, z), val) in a.indexed_iter() {
    //     println!("{:?} {:?} {:?}", x, y, z);
    // }

    // A == [1,0,0,0]
    // C == [0,1,0,0]
    // G == [0,0,1,0]
    // T == [0,0,0,1]

    // let mut arr = array![[b'A', b'T', b'C', b'C', b'A', b'G', b'C', b'T'],
    //                      [b'G', b'G', b'G', b'C', b'A', b'A', b'C', b'T'],
    //                      [b'A', b'T', b'G', b'G', b'A', b'T', b'C', b'T'],
    //                      [b'A', b'A', b'G', b'C', b'A', b'A', b'C', b'C'],
    //                      [b'T', b'T', b'G', b'G', b'A', b'A', b'C', b'T'],
    //                      [b'A', b'T', b'G', b'C', b'C', b'A', b'T', b'T'],
    //                      [b'A', b'T', b'G', b'G', b'C', b'A', b'C', b'T']];
    // print!("{:#?}", encode_output(&calc_consensus(&calc_profile(&encode_input(&arr)))));



}
