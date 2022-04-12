use std::collections::HashMap;
use std::collections::HashSet;
use std::cmp;
use allwords::{Alphabet};

use super::graph::ukonen::Ukonen;
use super::graph::ukonen::UkonenEdge;
use super::graph::ukonen::UkonenNode;

use crate::ds::sequence::Sequence;
use crate::ds::graph::Graph;
use crate::ds::tile::Tile;

use super::graph::trie::{Trie, TrieNode};

// Count number of chars in Sequence sequence
// Return array with numbers representing #occur of given char
// count[0] == count ['A'] and count[23] == count['Z']
pub fn count_nucleotides(seq: &Sequence) -> HashMap<u8, u16> {
    let mut count = HashMap::<u8, u16>::new();
    for c in seq.into_iter() {
        if !count.contains_key(c) {
            count.insert(*c, 1);
        }
        else {
            *count.get_mut(c).unwrap() += 1; 
        }
    }
    count
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

pub fn hamming_distance(s1: &Sequence, s2: &Sequence) -> u32 {
    let s1 = s1.to_string();
    let s2 = s2.to_string();
    s1.chars().zip(s2.chars()).filter(|&(a, b)| a != b).count() as u32
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
        (b'F', 2),   (b'I', 3),   (b'V', 4),   (b'L', 6),   
        (b'S', 6),   (b'P', 4),   (b'M', 1),   (b'T', 4),   
        (b'A', 4),   (b'Y', 2),   (b'-', 3),   (b'H', 2),   
        (b'N', 2),   (b'D', 2),   (b'Q', 2),   (b'K', 2),   
        (b'E', 2),   (b'C', 2),   (b'G', 4),   (b'R', 6),      
        (b'W', 1)
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
        (b'F', 147.06841),   (b'I', 113.08406),   (b'V', 99.06841),   (b'L', 113.08406),   
        (b'S', 87.03203),    (b'P', 97.05276),    (b'M', 131.04049),  (b'T', 101.04768),   
        (b'A', 71.03711),    (b'Y', 163.06333 ),  (b'-', 0.0),        (b'H', 137.05891),   
        (b'N', 114.04293),   (b'D', 115.02694),   (b'Q', 128.05858),  (b'K', 128.09496),   
        (b'E', 129.04259),   (b'C', 103.00919),   (b'G', 57.02146),   (b'R', 156.10111),      
        (b'W', 186.07931)
    ]);

    let mut mass : f64 = 0.0;
    for amino in protein.into_iter() {
        mass += monoisotopic_mass_table.get(amino).unwrap(); 
    }
    mass
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

pub fn rna_splice(mut pre_rna: Sequence, introns: &Tile) -> Sequence {

    for intr in introns {
        let res = knuth_morris_pratt(&pre_rna, intr);
        for index in res {
            pre_rna.chain.drain(index..(index + intr.len()));
        }
    }
    pre_rna
}

pub fn subsequences(a: &Sequence, b: &Sequence, limit: Option<usize>) -> Vec<Vec<usize>> {

    let mut result = vec![];
    let mut temp = Vec::<usize>::new();
    let a_idx: usize = 0;
    let b_idx: usize = 0; 

    pub fn subsequences_recursive( a: &Sequence, a_idx: usize, b: &Sequence, b_idx: usize, 
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

pub fn longest_increasing_subsequence(seq: &[u64]) -> Vec<u64> {

    let mut lis = vec![0; seq.len()];
    let mut pointers = vec![0; seq.len()];

    let mut max_idx : usize = 0;
    let mut max_len : usize = 0;

    lis[0] = 1;
    max_idx = 0;

    for i in 1..lis.len() {

        lis[i] = 1;
        pointers[i] = i;

        for j in 0..i {

            // TODO: Make more generic
            // Pass comparator as parameter
            if seq[i] > seq[j] && lis[i] < lis[j] + 1 {
            
                lis[i] = lis[j] + 1;
                pointers[i] = j;
            
                if lis[i] > max_len {
                    max_idx = i;
                    max_len = lis[i];
                }
            }
        }
    }

    let mut result : Vec<u64> = vec![];

    while max_len > 0 {
        result.push(seq[max_idx]);
        max_idx = pointers[max_idx];
        max_len -= 1;
    }

    result.reverse();
    result

}

pub fn sort_lexicographically(sequences: &Tile, alphabet: &[u8]) -> Tile {

    let mut trie_builder = Trie::new(alphabet);
    let trie = trie_builder.build(sequences).unwrap();

    fn walk_trie_rec(trie: &Graph<TrieNode, u8>, node_id: u64, sorted: &mut Tile) {
        
        if trie.out_neighbors(node_id).count() == 0 {
            return;
        }

        let node = trie.get_node(&node_id);
 
        for c in node.data.children.iter() {
            if *c != -1 {
                if trie.get_node(&(*c as u64)).data.ending == true {
                    let substrings = &trie.get_node(&(*c as u64)).data.substring;
                    substrings.iter().for_each(|s| sorted.push(Sequence::from(s)));
                }
    
                walk_trie_rec(trie, *c as u64, sorted);
            }
        }
    }

    let mut sorted = Tile::new();
    let root = trie.get_root().unwrap();
    walk_trie_rec(trie, root, &mut sorted);

    sorted
}

pub fn longest_common_subsequence(seq1: &Sequence, seq2: &Sequence) -> Sequence {

    let mut match_table = vec![vec![0_usize; seq2.len() + 1 ]; seq1.len() + 1];
    let mut prev_table = vec![vec![(0_usize, 0_usize); seq2.len() + 1 ]; seq1.len() + 1];

    for i in 1..(seq1.len()+1) {
        for j in 1..(seq2.len()+1) {
            if seq1[i-1] == seq2[j-1] {
                match_table[i][j] = match_table[i-1][j-1] + 1;
                prev_table[i][j] = (i-1, j-1);
            }
            else {
                if match_table[i-1][j] > match_table[i][j-1] {
                    match_table[i][j] = match_table[i-1][j];
                    prev_table[i][j] = (i-1, j);
                }
                else {
                    match_table[i][j] = match_table[i][j-1];
                    prev_table[i][j] = (i, j-1);  
                }
            }
        } 
    }

    let mut lcs = Sequence::new();

    let mut i = seq1.len();
    let mut j = seq2.len();

    while match_table[i][j] != 0 {
        let i_next = prev_table[i][j].0;
        let j_next = prev_table[i][j].1;
        if i_next==i-1 && j_next==j-1 {
            lcs.push(seq1[i_next]);
        }
        i = i_next;
        j = j_next
    }

    lcs.reverse();
    lcs
}

pub fn k_mer_composition(seq: &Sequence, k: usize, alphabet: &[u8]) -> Vec<usize> {

    // Generate all possible k-mers from alphabet
    let mut kmers = Tile::new();
    let a = Alphabet::from_chars_in_str(std::str::from_utf8(alphabet).unwrap()).unwrap();
    let words = a.all_words(Some(k)).filter(|x| x.len() == k);
    for a in words {
        kmers.push(Sequence::from(a));
    }

    // Sort k-mers according to ordering from alphabet
    kmers = sort_lexicographically(&kmers, alphabet);

    // Calculate k-mer composition
    let mut kmer_composition = vec![];
    for kmer in &kmers {
        let pos = knuth_morris_pratt(seq, kmer);
        kmer_composition.push(pos.iter().count());
    }

    kmer_composition
}

pub fn longest_common_supersequence(seq1: &Sequence, seq2: &Sequence) -> Sequence {

    let lcs = longest_common_subsequence(seq1, seq2);
    
    let mut superseq = Vec::<u8>::new();

    let mut s1 = seq1.into_iter();
    let mut s2 = seq2.into_iter();

    for c in lcs {

        loop {
            match s1.next() {
                Some(x) if *x != c => {
                    superseq.push(*x);
                }
                _ => {
                    break;
                }
            }
        }

        loop {
            match s2.next() {
                Some(x) if *x != c => {
                    superseq.push(*x);
                }
                _ => {
                    break;
                }
            }
        }

        superseq.push(c);
    }

    superseq.extend(s1);
    superseq.extend(s2);

    Sequence::from(superseq.as_slice())
}


pub fn protein_from_prefix_spectrum(spec: Vec<f32>) -> Sequence {

    let monoisotopic_mass_table : HashMap<u8, f32> = HashMap::from([   
        (b'F', 147.06841),   (b'I', 113.08406),   (b'V', 99.06841),   (b'L', 113.08406),   
        (b'S', 87.03203),    (b'P', 97.05276),    (b'M', 131.04049),  (b'T', 101.04768),   
        (b'A', 71.03711),    (b'Y', 163.06333 ),  (b'-', 0.0),        (b'H', 137.05891),   
        (b'N', 114.04293),   (b'D', 115.02694),   (b'Q', 128.05858),  (b'K', 128.09496),   
        (b'E', 129.04259),   (b'C', 103.00919),   (b'G', 57.02146),   (b'R', 156.10111),      
        (b'W', 186.07931)
    ]);

    let mut result = Sequence::new();

    for i in 1..spec.len() {
        let mut diff = spec[i] - spec[i-1];
        let x = monoisotopic_mass_table.iter()
                        .find(|(key, value)| (*value - diff).abs() < 0.01).unwrap();
        result.push(*x.0);
    }

    result
}

pub fn edit_distance(seq1: &Sequence, seq2: &Sequence) -> usize {

    let mut memo = vec![vec![0_usize; seq2.len() + 1 ]; seq1.len() + 1];

    // initialize table
    for i in 0..seq1.len() {
        memo[i][0] = i;
    }
    for j in 0..seq2.len() {
        memo[0][j] = j;
    }

    for i in 1..seq1.len()+1 {
        for j in 1..seq2.len()+1 {
            if seq1[i-1] == seq2[j-1] {
                memo[i][j] = memo[i-1][j-1];
            }
            else {
                memo[i][j] = cmp::min(memo[i-1][j-1], cmp::min(memo[i][j-1], memo[i-1][j])) + 1;
            }
        }
    }
    
    memo[seq1.len()][seq2.len()]
}




// Greedy search for shortest common superstring
// pub fn shortest_common_superstring(sequences: &Tile) -> Sequence {

//     let num_seq = sequences.len();

//     for i in 0..num_seq {
//         let prefix = &sequences[i];
//         for j in 0..num_seq {
            
//             if i == j {
//                 continue;
//             }
            
//             let suffix = &sequences[j];
//             let num_chars = std::cmp::min(prefix.len(), suffix.len());
            
//             for k in 0..num_chars {
//                 if prefix[k] == suffix[num_chars-k]
//             }



//         }
//     }
    

//     Sequence::new()
// }


// pub fn longest_common_substring(matrix: &Tile) {

//     // TODO: include information about alphabet inside sequence itself
//     let mut alphabet = HashSet::<u8>::from([
//         'A' as u8 , 'C' as u8, 
//         'T' as u8, 'G' as u8
//     ]);

//     let mut separator : u8 = '!' as u8; 
//     let mut temp : Vec<u8> = vec![];
//     let mut wordmap : Vec<(usize, usize)> = vec![];

//     for (idx, a) in matrix.into_iter().enumerate() {
//         temp.extend(a);
//         temp.push(separator);
//         wordmap.extend(vec![(idx, temp.len() - 1); a.len() + 1]);
//         separator += 1;
//         while alphabet.contains(&separator) {
//             separator += 1;
//         }
//     }


//     let seq = Sequence::from(temp.as_slice());
//     let mut ukkokens = Ukonen::<Sequence>::new(seq);
//     let g = ukkokens.process();
//     g.write_dot("abc.dot");


//     fn generate_reachbility_map(graph: &mut Graph<UkonenNode, UkonenEdge>, node_id: u64, 
//                                 discovered: &mut HashSet<u64>, wordmap: &Vec<(usize, usize)>, 
//                                     reachable_suffixes: &mut HashMap<u64, Vec<u8>>)
//     {
//         discovered.insert(node_id);
//         let out_neighbors : Vec<u64> = graph.out_neighbors(node_id).cloned().collect();

//         for t in out_neighbors{
//             if !discovered.contains(&t) {
//                 generate_reachbility_map(graph, t, discovered, wordmap, reachable_suffixes);
//             }
//         }

//         let mut reach = vec![0; wordmap.last().unwrap().0 + 1];
//         let out_edges : Vec<u64> = graph.out_edges(node_id).cloned().collect();

//         for eid in out_edges{

//             let successor_node_id = graph.get_edge(&eid).end;
//             let suffix_start = graph.get_edge(&eid).data.as_ref().unwrap().suffix_start;
//             let suffix_stop = graph.get_edge(&eid).data.as_ref().unwrap().suffix_stop;

//             if suffix_stop == -1 {
//                 reach[wordmap[suffix_start].0] = 1;
//             }
//             else {

//                 for (i, elem) in reachable_suffixes[&successor_node_id].iter().enumerate() {
//                     if *elem == 1 {
//                         reach[i] = 1;
//                     }
//                 };
//             }
//         }

//         reachable_suffixes.insert(node_id, reach);
//     }

//     let mut visited = HashSet::<u64>::new();
//     let mut reachable_suffixes = HashMap::<u64, Vec<u8>>::new();
//     generate_reachbility_map(g, g.get_root().unwrap(), &mut visited, &wordmap, &mut reachable_suffixes);

//     // // Function to perform DFS traversal on the graph
//     fn resolve_suffix_endings(graph: &mut Graph<UkonenNode, UkonenEdge>, node_id: u64, 
//                                 discovered: &mut HashSet<u64>, wordmap: &Vec<(usize, usize)>)
//     {
//         discovered.insert(node_id);

//         let out_neighbors : Vec<u64> = graph.out_neighbors(node_id).cloned().collect();
//         for t in out_neighbors{
//             if !discovered.contains(&t) {
//                 resolve_suffix_endings(graph, t, discovered, wordmap);
//             }
//         }

//         let out_edges : Vec<u64> = graph.out_edges(node_id).cloned().collect();
//         for eid in out_edges {            
//             let suffix_start = graph.get_edge(&eid).data.as_ref().unwrap().suffix_start;
//             let suffix_stop = graph.get_edge(&eid).data.as_ref().unwrap().suffix_stop;
//             if suffix_stop == -1 {
//                 graph.get_edge_mut(&eid).data.as_mut().unwrap().suffix_stop = wordmap[suffix_start].1 as i64;
//             }
//         }
//     }

//     visited.clear();
//     resolve_suffix_endings(g, g.get_root().unwrap(), &mut visited, &wordmap);

//     let mut lcs : Vec<(usize, i64)> = vec![];
//     let mut cur : Vec<(usize, i64)> = vec![];

//     let mut longest = 0;
//     let mut cur_len = 0;

//     fn dfs_recursive_substrings(graph: &Graph<UkonenNode, UkonenEdge>, node_id: u64, 
//                                     discovered: &mut HashSet<u64>,  reachable_suffixes: & HashMap<u64, Vec<u8>>,
//                                         cur_suffix: &mut Vec<(usize, i64)>, cur_length: usize, 
//                                             lcs: &mut Vec<(usize, i64)>, longest: &mut usize, len: u8)
//     {
//         // mark the current node as discovered
//         discovered.insert(node_id);

//         if cur_length > *longest {
//             *longest = cur_length;
//             *lcs = cur_suffix.clone();
//         } 
        
//         // println!("Current length: {:#?}", reachable_suffixes);

//         // do for every edge (v, u)
//         for e in graph.out_edges(node_id){
//             if !discovered.contains(&graph.get_edge(e).end) && (reachable_suffixes[&graph.get_edge(e).end].iter().sum::<u8>() == len  ){
//                 let start = graph.get_edge(e).data.as_ref().unwrap().suffix_start;
//                 let stop = graph.get_edge(e).data.as_ref().unwrap().suffix_stop;
//                 cur_suffix.push((start, stop));
//                 dfs_recursive_substrings(graph, graph.get_edge(e).end, discovered, reachable_suffixes, cur_suffix, cur_length + stop as usize - start + 1, lcs, longest, len);
//                 cur_suffix.pop();
//             }
//         }


//     }

//     visited.clear();
//     dfs_recursive_substrings(g, g.get_root().unwrap(), &mut visited, &reachable_suffixes, &mut cur, cur_len, &mut lcs, &mut longest, matrix.size().0 as u8);



//     let mut result = vec![];

//     println!("{:?}", lcs);

//     for part in lcs {

//         for i in part.0..(part.1 as usize + 1) {
//             result.push(temp[i] as u8);
//         }
//     }

//     // let x =  Sequence::from(result.as_slice());

//     // print!("{}", x);

// }


// Independent Alleles
// sum (p choose i)*0.25^i*0.75^(p-i), i from n to p for p=2^6, n=17






#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_count_nucleotides() {
    //     let input = Sequence::from("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCT\
    //                               CTGTGTGGAATTAAAAAAAGAGTGTCTGATGCAGC");
    //     assert_eq!([20, 12, 17, 21], count_nucleotides(&input));
    // }

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
    fn test_substring_positions() {
        let seq = Sequence::from("GATATATGCATATACTT");
        let pat = Sequence::from("ATAT");
        assert_eq!(vec![1, 3, 9], knuth_morris_pratt(&seq, &pat));
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

    #[test]
    fn test_shortest_supersequence() {
        let s2 = Sequence::from("TTATGTGATATCCCCGCTTCTCACAATGCTCTTAGTTTACCTC\
                                 GAACTAAGTCTGATCGCAGCGGCCGGTATTCCTTTCTACGCG");
        let s1 = Sequence::from("GCTCACGGATTCGAAAGTCGAGTGTCCCCCAGCTGGATGCATTCTT\
                                 TGGGAGTGGCCAAGGAGGGTTATCAGAAGAACAGATTAATTTG");
        let res = Sequence::from("GCTCTACTGTGATATCCCCGCTTCTCACAATGCTCGTTAGTGT\
                                  TACCTCCCAGAACTGGATGCAGTTCTTTGGGAGTGGCGCAAGC\
                                  GAGCCGGTTATTCAGAAGAACAGATTAATCTTACGCG");
        assert_eq!(res, longest_common_supersequence(&s1, &s2));
    }
}