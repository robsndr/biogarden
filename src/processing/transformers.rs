use std::collections::HashMap;
use std::collections::HashSet;
use allwords::{Alphabet};

use crate::ds::sequence::Sequence;
use crate::ds::graph::Graph;
use crate::ds::tile::Tile;
use crate::analysis;
use crate::processing::patterns::knuth_morris_pratt;
use crate::ds::builders::trie::{Trie, TrieNode};

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
                    println!("value: {}", &chunk);
                    println!("Codon not found in codon table.");
                    break;
                }
                
            }
        }
    }
    proteins
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


pub fn rna_splice(mut pre_rna: Sequence, introns: &Tile) -> Sequence {

    for intr in introns {
        let res = knuth_morris_pratt(&pre_rna, intr);
        for index in res {
            pre_rna.chain.drain(index..(index + intr.len()));
        }
    }
    pre_rna
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

/// Returns tuples denoting error corrections that can be applied to a number of input reads.
///   
/// Genome sequencers use a chemical procedure to obtain reads from provided biomaterial.
/// Due to the error-prone character of this approach, multiple reads of the same region are usually taken.
/// Errors can be corrected by splitting the obtained set into correct and faulty reads first.
/// If a read or its complement appears in the set at least `split_margin` times, it is regarded as correct.
/// Faulty reads are corrected based on their hamming distance to one of the correct reads.
/// 
/// # Arguments
///
/// * `reads` - tile containing the analyzed reads
/// * `split_margin` - number of reads/complements required to treat a read as correct
/// * `hd_margin` - hamming distance used for matching faulty reads to correct ones
/// 
pub fn correct_read_errors(reads: &Tile, split_margin: usize, hd_margin: usize) -> Vec<(Sequence, Sequence)> {
    
    // Count number of repeated reads and complements 
    let mut read_counter  = HashMap::<Sequence, usize>::new();
    for read in reads {
        // Insert read if not present already and increment count
        *read_counter.entry(read.clone()).or_insert(0) += 1;
        // Handle complement case
        let complement = complement_dna(read.clone());
        if read_counter.contains_key(&complement) {
            *read_counter.get_mut(&complement).unwrap() += 1;
            *read_counter.get_mut(read).unwrap() += 1;
        }
    }

    // Split according to split margin
    let mut correct_reads = HashSet::<Sequence>::new();
    let mut faulty_reads = HashSet::<Sequence>::new();
    read_counter.iter().for_each(|(fr, cnt)| {
        if *cnt >= split_margin {  
            correct_reads.insert(fr.clone()) 
        } else { 
            faulty_reads.insert(fr.clone()) 
        };
    });
    
    // Compute corrections satisfying hamming margin, applicable to faulty reads 
    let mut corrections = Vec::<(Sequence, Sequence)>::new();
    for fr in faulty_reads.iter() {
        // Find correct reads/complements satisfying `H(x) <= hamming_distance_margin`
        for cr in correct_reads.iter() {
            // H(x) <= hamming_distance_margin
            if analysis::seq::hamming_distance(fr, cr).unwrap() <= hd_margin {
                corrections.push((fr.clone(), cr.clone()));
                break;
            }
            // H(complement(x)) <= hamming_distance_margin
            let complement = complement_dna(cr.clone());
            if analysis::seq::hamming_distance(fr, &complement).unwrap() == hd_margin {
                corrections.push((fr.clone(), complement));
                break;
            }
        }
    }

    corrections
}



pub fn generate_k_mers(seq: &Sequence, k: usize) -> Vec<Sequence> {

    // TODO: return proper error
    if k > seq.len() {
        panic!("LEN(SEQ1) !< LEN(SEQ2)");
    }

    (0..seq.len()-k+1).into_iter()
                    .map(|i| Sequence::from(&seq[i..k+i]) )
                    .collect()
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
