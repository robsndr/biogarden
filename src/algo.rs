use std::io;
use ndarray::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt; // Import `fmt`
use statrs::function::factorial::{binomial, factorial};

use super::Sequence;
use super::Tile;
use super::Graph;
use super::GraphProperties;
use super::Dfs;

// Count number of chars in Sequence sequence
// Return array with numbers representing #occur of given char
// count[0] == count ['A'] and count[23] == count['Z']
pub fn count_nucleotides(seq: &Sequence) -> [u16; 4] {
    let mut count: [u16; 24] = [0; 24];
    for c in seq.into_iter() {
        count[(*c as usize) - 65] += 1;
    }
    return [count[0], count[2], count[6], count[19]];
}

// Transcribe the DNS sequence into RNA
pub fn transcribe_dna(dna: Sequence) -> Sequence {
    let temp = String::from(dna);
    Sequence::from(temp.replace("T", "U"))
}

// Complement the Sequence string by reversing in the first step.
// Swap: 'A' <-> 'T' and 'G' <-> 'C'
pub fn complement_dna(seq: Sequence) -> Sequence {
    let mut t = seq.into_iter().rev().map(|c| c as char).collect::<String>();
    // A <-> T
    t = t.replace("A", "X");
    t = t.replace("T", "A");
    t = t.replace("X", "T");
    // G <-> C
    t = t.replace("G", "X");
    t = t.replace("C", "G");
    t = t.replace("X", "C");
    Sequence::from(t)
}

// Percentage of G/C nucleotides in sequence
// Return percentage value
pub fn gc_content(seq: &Sequence) -> f64 {
    let mut gc_count : u32 = 0;
    for c in seq {
        if *c == b'G' ||  *c == b'C' {
            gc_count += 1;
        }
    }
    gc_count as f64 / seq.len() as f64
}

// Make parametrizable number of sequences to find
pub fn translate_rna(rna: Sequence) -> Vec<Sequence> {

    let mut proteins : Vec<Sequence> = vec![];

    // mRNA <-> amino-acid translation table (codon table)
    let codon_table = HashMap::from([   
        ("UUU", "F"),    ("CUU", "L"),   ("AUU", "I"),   ("GUU", "V"),
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
        ("UGG", "W"),    ("CGG", "R"),   ("AGG", "R"),   ("GGG", "G") 
    ]);
    // Container for final result of transcription
    let mut amino_acid = String::new();
    // Run the translation 
    let s = String::from(rna);
    let mut z = s.chars().peekable();
    // Iterate until end of strand is reached
    while z.peek().is_some() {
        amino_acid.clear();
        // Iterate over strand until start codon found
        while z.peek().is_some() {
            // Take 3 characters from strand, that denote codon
            let chunk: String = z.by_ref().take(3).collect();
            // Check for start codon
            if chunk == "AUG"{
                amino_acid.push_str(codon_table.get(&chunk as &str).unwrap());
                break;
            }
        }
        // Copy current iterator to resume search for start codon at that position
        let mut zi = z.clone(); 
        // Decode until stop codon reached
        while zi.peek().is_some() {
            // Take 3 characters from strand, that denote codon
            let chunk: String = zi.by_ref().take(3).collect();
            match codon_table.get(&chunk as &str) {
                Some(value) => {
                    // If stop codon reached, store current protein strand and proceed 
                    if value == &"Stop"{
                        proteins.push(Sequence::from(amino_acid.clone()));
                        break;
                    }
                    else {
                        amino_acid.push_str(value);
                    }
                },
                None => {
                    print!("value: {}\n", &chunk);
                    println!("Codon not found in codon table.");
                    break;
                }
                
            }
        }
    }
    proteins
}

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

pub fn hamming_distance(s1: &Sequence, s2: &Sequence) -> u32 {
    let s1 = s1.to_string();
    let s2 = s2.to_string();
    s1.chars().zip(s2.chars()).filter(|&(a, b)| a != b).count() as u32
}

pub fn knuth_morris_pratt(seq: &Sequence, pat: &Sequence) -> Vec<usize> {
   
    let seq = seq.to_string().into_bytes();
    let pat = pat.to_string().into_bytes();

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

    for (i, &c) in seq.iter().enumerate() {
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

// Find all reverse-palindromes within seq of n <= length <= m
// Return tuples containing position and length of each palindrome O(n^3)?
pub fn reverse_palindromes(seq: &Sequence, n: usize, m: usize) -> Vec<(usize, usize)>{

    let mut palindromes : Vec<(usize, usize)> = vec![];
    let complements = HashMap::from([(b'A', b'T'), (b'T', b'A'),
                                     (b'G', b'C'), (b'C', b'G')]);

    // iterate over every offset within the initial string
    for i in 0..seq.len() {
        // iterate over possible lengths of palindromic substrings
        for j in n..(m+1) {
            // break if potential substring cannot fit 
            if i + j > seq.len() {
                break;
            } 
            // check if substring with length `j` at offset `i` 
            // is a reverse palindrome
            let mut is_palindrome = true;
            for k in 0..j {
                if seq.chain[i+k] != complements[&seq.chain[i+j-1-k]] {
                    is_palindrome = false;
                    break;
                }
            }
            // append (offset, length) into result set
            if is_palindrome {
                palindromes.push((i+1,j));
            }
        }    
    }
    palindromes
}

pub fn rna_splice(mut pre_rna: Sequence, introns: &Tile) -> Sequence {

    for intr in introns {
        let res = knuth_morris_pratt(&pre_rna, intr);
        for index in res {
            pre_rna.chain.drain(index..(index + intr.len()));
        }
    }
    pre_rna
}

pub fn open_reading_frames(dna: &Sequence) -> Vec<Sequence> {

    let mut reading_frames : Vec<Sequence> = vec![];

    let mut strands : Vec<Sequence> = vec![];
    strands.push(dna.clone());
    strands.push(complement_dna(dna.clone()));

    for strand in &strands {
        for i in 0..3 {
            let mut temp: Sequence = strand.clone();
            temp.chain.drain(0..i);
            temp = transcribe_dna(temp);
            reading_frames.extend(translate_rna(temp));
        }    
    }
    reading_frames
}

pub fn infer_number_rna(protein: &Sequence) -> u128 {

    let codon_combs: HashMap<u8, u128> = HashMap::from([   
        ('F' as u8, 2),   ('I' as u8, 3),   ('V' as u8, 4),   ('L' as u8, 6),   
        ('S' as u8, 6),   ('P' as u8, 4),   ('M' as u8, 1),   ('T' as u8, 4),   
        ('A' as u8, 4),   ('Y' as u8, 2),   ('-' as u8, 3),   ('H' as u8, 2),   
        ('N' as u8, 2),   ('D' as u8, 2),   ('Q' as u8, 2),   ('K' as u8, 2),   
        ('E' as u8, 2),   ('C' as u8, 2),   ('G' as u8, 4),   ('R' as u8, 6),      
        ('W' as u8, 1)
    ]);

    // Initialize with 3 as for number of STOP codons
    let mut rna_combinations : u128 = 3;
    // Compute number of combinations
    for amino in protein {
        rna_combinations = (rna_combinations * codon_combs.get(amino).unwrap()) % 1000000;
    }   
    rna_combinations
}

pub fn weighted_mass(protein: &Sequence) -> f64 {

    let monoisotopic_mass_table : HashMap<u8, f64> = HashMap::from([   
        ('F' as u8, 147.06841),   ('I' as u8, 113.08406),   ('V' as u8, 99.06841),   ('L' as u8, 113.08406),   
        ('S' as u8, 87.03203),    ('P' as u8, 97.05276),    ('M' as u8, 131.04049),  ('T' as u8, 101.04768),   
        ('A' as u8, 71.03711),    ('Y' as u8, 163.06333 ),  ('-' as u8, 0.0),        ('H' as u8, 137.05891),   
        ('N' as u8, 114.04293),   ('D' as u8, 115.02694),   ('Q' as u8, 128.05858),  ('K' as u8, 128.09496),   
        ('E' as u8, 129.04259),   ('C' as u8, 103.00919),   ('G' as u8, 57.02146),   ('R' as u8, 156.10111),      
        ('W' as u8, 186.07931)
    ]);

    let mut mass : f64 = 0.0;
    for amino in protein.into_iter() {
        mass += monoisotopic_mass_table.get(amino).unwrap(); 
    }
    mass
}

pub fn random_substrings(seq: &Sequence, gc_content: &[f64]) -> Vec<f64> {

    let mut probabilities : Vec<f64> = vec![];
    let mut prob_map : HashMap<u8, f64> = HashMap::from([ ('A' as u8, 1.0), ('C' as u8, 1.0), ('T' as u8, 1.0), ('G' as u8, 1.0)]);

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

pub fn overlap_graph(sequences: &Tile, k: usize) -> Graph::<Sequence, u8> {

    // Instantiate empty graph
    // let gp = ;
    let mut g = Graph::<Sequence, u8>::new(GraphProperties{directed: true});
    
    // Add all nodes to  graph
    let mut node_ids : Vec<u64> = vec![];
    for seq in sequences {
        node_ids.push(g.add_node(seq.clone()));
    }

    // Connect overlap graph
    for (i, seq) in sequences.into_iter().enumerate() {
        let last = seq.into_iter().rev().take(k).rev();
        for (j, seq2) in sequences.into_iter().enumerate() {
            if j != i {
                let first = seq2.into_iter().take(k);
                // Check if suffix of `seq` is equal to prefix of `seq2`
                if first.zip(last.clone()).filter(|&(a, b)| a != b).count() == 0 {
                    g.add_edge(&node_ids[i], &node_ids[j], None).unwrap();
                }
            }
        }
    }
    g
}

pub fn partial_permutations(n: u64, k: u64) -> u64 {
    let combinations = binomial(n, k);
    let permutations = factorial(k);
    (combinations * permutations) as u64 % 1000000
}

pub fn transition_transversion_ratio(s1: &Sequence, s2: &Sequence) -> f32 {
    let s1 = s1.to_string();
    let s2 = s2.to_string();
    let missmatch_iter = s1.chars().zip(s2.chars()).filter(|&(a, b)| a != b);
    let hamming = missmatch_iter.clone().count() as f32;
    let mut transitions : f32 = 0.0;
    for (a, b) in missmatch_iter {
        // C <-> T 
        if (a == 'C' && b == 'T') || (a == 'T' && b == 'C'){
            transitions += 1.0;
        }
        // A <-> G
        if (a == 'A' && b == 'G') || (a == 'G' && b == 'A'){
            transitions += 1.0;
        }
    }
    let transversions = hamming  - transitions;
    transitions / transversions
}

pub fn connected_components(g: &Graph<u64, u8>) -> (u32, Vec<Vec<u64>>) {

    let mut dfs = Dfs::new(g);

    let mut processed = HashSet::<u64>::new();
    let mut ctr : u32 = 0;
    let mut cc : Vec<u64> = vec![];
    let mut components : Vec<Vec<u64>> = vec![];
    
    for n in g.nodes() {
        if !processed.contains(n) {
            dfs.init(n);
            cc.clear();
            while let Ok(id) = dfs.process_node() {
                processed.insert(id);
                cc.push(id);
            }
            components.push(cc.clone());
            ctr += 1;
        }
    }
    (ctr, components)
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

// Greedy search for shortest common superstring
pub fn subsequences(a: &Sequence, b: &Sequence, limit: Option<usize>) -> Vec<Vec<usize>> {

    let mut result = vec![];
    let mut temp = Vec::<usize>::new();
    let a_idx: usize = 0;
    let b_idx: usize = 0; 

    pub fn subsequences_recursive(a: &Sequence, a_idx: usize,b: &Sequence, b_idx: usize, 
                                     temp: &mut Vec<usize>, result: &mut Vec<Vec<usize>>,
                                     limit: Option<usize>) 
    {
        if b_idx == b.len() {
            result.push(temp.clone());
            return;
        }
        for i in a_idx..a.len() {
            if limit.is_some() && result.len() == limit.unwrap() {
                return;
            }
            if b[b_idx] == a[i] {
                temp.push(i);
                subsequences_recursive(a, i+1, b, b_idx+1, temp, result, limit);
                temp.pop();
            }
        }   
    }

    subsequences_recursive(a, a_idx, b, b_idx, &mut temp, &mut result, limit);
    return result;
}


// Greedy search for shortest common superstring
// pub fn shortest_common_superstring(sequences: &Tile) -> Sequence {

//     let mut subsequence_set = HashSet::<Sequence>::new();
//     let mut k : usize = 0;
//     for seq in sequences {
//         subsequence_set.insert(seq.clone());
//         // find max length of subsequence used
//         if seq.len() > k {
//             k = seq.len();
//         }
//     }

//     k = 150;

//     while k > 1 {
//         let grph = overlap_graph(sequences,  k);
//         if grph.edge_count() > 0 {
//             let e = grph.edges().next().unwrap();
//             // print!("{:#?}", e);
//         }
//         // print!("ALOHA");
//         k = k-1;
//     }

//     Sequence::new()
// }


// N{P}[ST]{P}
// pub fn generate_motifs(target: &Vec<u8>, result: &mut Vec<String>, i: usize, temp: &mut Vec<u8>) {
//
//     if i == target.len() {
//         result.push(std::str::from_utf8(temp).unwrap().to_string());
//         return;
//     }
//
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
//
//     let st = st.into_bytes();
//     let pat = pat.into_bytes();
//     let mut res = Vec::<String>::new();
//     let mut temp = Vec::<u8>::new();
//
//     generate_motifs(&pat, &mut res, 0, &mut temp);
//
//     let mut pos = Vec::<usize>::new();
//     for elem in &res {
//         pos.append(&mut knuth_morris_pratt(&st, elem.as_bytes()));    
//     }
//     pos
// }


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_nucleotides() {
        let input = Sequence::from("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCT\
                                  CTGTGTGGAATTAAAAAAAGAGTGTCTGATGCAGC");
        assert_eq!([20, 12, 17, 21], count_nucleotides(&input));
    }

    #[test]
    fn test_transcribe_dna() {
        let input = Sequence::from("GATGGAACTTGACTACGTAAATT");
        let result = Sequence::from("GAUGGAACUUGACUACGUAAAUU");
        assert_eq!(result, transcribe_dna(input))
    }

    #[test]
    fn test_complement_dna() {
        let input = Sequence::from("AAAACCCGGT");
        let result = Sequence::from("ACCGGGTTTT");
        assert_eq!(result, complement_dna(input));
    }

    #[test]
    fn test_gc_content() {
        let input = Sequence::from("CCTGCGGAAGATCGGCACTAGAATAGCCAG\
                                    AACCGTTTCTCTGAGGCTTCCGGCCTTCCC");
        let result : f64 = 0.5833333333333334;
        assert_eq!(result, gc_content(&input));
    }

    #[test]
    fn test_translate_rna() {
        let input = Sequence::from("AUGGCCAUGGCGCCCAGAACUGAGA\
                                    UCAAUAGUACCCGUAUUAACGGGUGA");
        let result = Sequence::from("MAMAPRTEINSTRING");
        assert_eq!(result, *translate_rna(input).first().unwrap());
    }

    #[test]
    fn test_hamming_distance() {
        let input1 = Sequence::from("GAGCCTACTAACGGGAT");
        let input2 = Sequence::from("CATCGTAATGACGGCCT");
        let result : u32 = 7;
        assert_eq!(result, hamming_distance(&input1, &input2));

    }

    #[test]
    fn test_mendel_first_law() {
        let k: u16 = 2;
        let m: u16 = 2;
        let n: u16 = 2;
        assert_eq!(0.7833333, 
                    mendel_first_law(k, m, n));
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
    fn test_fibo() {
        let n: usize = 5;
        let k: u64 = 3;
        assert_eq!(19, fibo(n, k));
    }

    #[test]
    fn test_substring_positions() {
        let seq = Sequence::from("GATATATGCATATACTT");
        let pat = Sequence::from("ATAT");
        assert_eq!(vec![1, 3, 9], knuth_morris_pratt(&seq, &pat));
    }

    #[test]
    fn test_random_substring() {
        let seq = Sequence::from("GACGGACAAGGGCCCCCGTGTATTGTACTTGGGCCCATGTGCC\
                                    CGACCTCGGTAAGTCCATCAGGAGTGCACGAGGACCACCATTTCAAGAAA");
        let arr =  [0.110, 0.127, 0.183, 0.256];
        assert_eq!(vec![ -80.82637701756781, -77.85362263786175, 
                            -70.59699957974087, -64.49615401338707], random_substrings(&seq, &arr));
    }

    #[test]
    fn test_partial_permutations() {
        let part = partial_permutations(21, 7);
        assert_eq!(51200, part);
    }

    #[test]
    fn test_transition_transversion_ratio() {
        let a = Sequence::from("GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGT\
                                GATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT");
        let b = Sequence::from("TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGC\
                                GATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT");
        let ratio = transition_transversion_ratio(&a, &b);
        assert_eq!(1.21428571429, ratio);
    }
}