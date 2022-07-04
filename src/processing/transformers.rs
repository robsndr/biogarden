use allwords::Alphabet;
use std::collections::HashMap;
use std::collections::HashSet;

use crate::ds::builders::trie::{Trie, TrieNode};
use crate::ds::graph::Graph;
use crate::ds::sequence::Sequence;
use crate::ds::tile::Tile;
use crate::error::{BioError, Result};
use crate::processing::constants;
use crate::processing::patterns::find_motif;
use crate::analysis::seq::hamming_distance;

/// Transcribe the DNS sequence into RNA
///
/// During RNA transcription a strand of DNA is converted into RNA, where thymine (T) is replaced by uracil (U).
///
/// # Arguments
/// * `dna` - DNA string to transcribe into RNA
///
/// # Example
/// ```
/// use biotech::processing::transformers::transcribe_dna;
/// use biotech::ds::sequence::Sequence;
///
/// let dna = Sequence::from("GATGGAACTTGACTACGTAAATT");
/// let rna = Sequence::from("GAUGGAACUUGACUACGUAAAUU");
///
/// assert_eq!(transcribe_dna(dna), rna);
/// ```
pub fn transcribe_dna(dna: Sequence) -> Sequence {
    dna.into_iter()
        .map(|x| if x == b'T' { b'U' } else { x })
        .collect()
}

/// Complement DNA by reversing in the first step.
///
/// In the double helix, adenine (A) always bonds with Thymine (T), and cytosine (C) always bonds with guanine (G).
/// To generate the complementary strand of the primary, a strand must be reversed and bases must be swapped: A-T and G-C.
///
/// # Arguments
/// * `dna` - DNA string to transcribe into RNA
///
/// # Example
/// ```
/// use biotech::processing::transformers::complement_dna;
/// use biotech::ds::sequence::Sequence;
///
/// let seq = Sequence::from("AAAACCCGGT");
/// let complement = Sequence::from("ACCGGGTTTT");
///
/// assert_eq!(complement_dna(seq), complement);
/// ```
pub fn complement_dna(dna: Sequence) -> Sequence {
    dna.into_iter()
        .rev()
        .map(|x| match x {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            _ => x,
        })
        .collect()
}

/// Translate an RNA sequence into proteins
///
/// One of the central dogma in molecular biology is that DNA is transcribed into RNA, which in turn is translated into Protein.
/// The regions of DNA/RNA that encode protein are of particular interest, with groups of 3 nucleobases (codons) mapping to amino acids.
/// There are two special types of codons: the `AUG` start codon is used to indicate the start of translation,
/// with `UAA`, `UAG`, `UGA` stop codons used to mark where translation should end.
///
/// # Arguments
/// * `rna` - RNA sequence to translate into proteins
/// * `limit` - maximum number of protein to decode (many protein possible in one RNA)
///
/// # Example
/// ```
/// use biotech::processing::transformers::translate_rna;
/// use biotech::ds::sequence::Sequence;
/// use biotech::ds::tile::Tile;
///
/// let rna = Sequence::from("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA");
/// let mut result = Tile::new();
/// result.push(Sequence::from("MAMAPRTEINSTRING"));
/// result.push(Sequence::from("MAPRTEINSTRING"));
///
/// assert_eq!(translate_rna(rna, Some(2)), result);
/// ```
pub fn translate_rna(rna: Sequence, limit: Option<usize>) -> Tile {
    // TODO: Refactor into iterator/generator
    let mut proteins = Tile::new();

    // Container for final result of transcription
    let mut amino_acid = String::new();
    // Run the translation
    // TODO: Do not use String type
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
            if chunk == "AUG" {
                amino_acid.push_str(constants::CODON_TABLE.get(&chunk as &str).unwrap());
                break;
            }
        }
        // Copy current iterator to resume search for start codon at that position
        let mut zi = z.clone();
        // Decode until stop codon reached
        while zi.peek().is_some() {
            // Take 3 characters from strand, that denote codon
            let chunk: String = zi.by_ref().take(3).collect();
            match constants::CODON_TABLE.get(&chunk as &str) {
                Some(value) => {
                    // If stop codon reached, store current protein strand and proceed
                    if value == &"Stop" {
                        proteins.push(Sequence::from(amino_acid.clone()));
                        break;
                    } else {
                        amino_acid.push_str(value);
                    }
                }
                None => {
                    println!("value: {}", &chunk);
                    println!("Codon not found in codon table.");
                    break;
                }
            }
        }
        // Check if the upper limit of regions to translate has been reached
        if let Some(l) = limit {
            if l <= proteins.len() {
                break;
            };
        }
    }
    proteins
}

/// Decode protein directly from a DNA sequence
///
/// DNA can be transcribed into RNA, which in turn might be translated to protein.
/// As a result, one might consider the direct transcription of DNA into protein.
///
/// # Arguments
/// * `dna` - dna sequence to decode into proteins
///
/// # Example
/// ```
/// use biotech::processing::transformers::open_reading_frames;
/// use biotech::ds::sequence::Sequence;
/// use biotech::ds::tile::Tile;
///
/// let dna = Sequence::from("AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACT\
///                             TGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG");
///
/// let mut result = Tile::new();
/// result.push(Sequence::from("MGMTPRLGLESLLE"));
/// result.push(Sequence::from("MTPRLGLESLLE"));
/// result.push(Sequence::from("M"));
/// result.push(Sequence::from("M"));
/// result.push(Sequence::from("MLLGSFRLIPKETLIQVAGSSPCNLS"));
///
/// let orf = open_reading_frames(&dna);
///
/// assert_eq!(orf, result);
/// ```
pub fn open_reading_frames(dna: &Sequence) -> Tile {
    let mut reading_frames = Tile::new();
    for strand in [dna.clone(), complement_dna(dna.clone())] {
        // Translation of DNA into Protein might begin
        // at any position, so consider all start positions for codons (len==4)
        for i in 0..3 {
            let mut temp = strand.clone();
            temp.chain.drain(0..i);
            temp = transcribe_dna(temp);
            reading_frames.extend(translate_rna(temp, None));
        }
    }
    reading_frames
}

/// Remove introns from a sequence of pre-RNA
///
/// A RNA sequence can be subdivided into a set of introns and exons.
/// Introns are intervals that are ignored during protein translation, such that they don't influence the
/// expression of a gene. Exons on the other hand, are bases that form coding region of a gene.
/// During `splicing`, introns are removed from a sequence, with exons stitched together as a result.
///
/// # Arguments
/// * `pre_rna` - preRNA to splice
/// * `intons` - introns to be spliced-out
///
/// # Example
/// ```
/// use biotech::processing::transformers::splice_introns;
/// use biotech::ds::sequence::Sequence;
/// use biotech::ds::tile::Tile;
///
/// let prna = Sequence::from("ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTC\
///                             GAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG");
///
/// let mut introns = Tile::new();
/// introns.push(Sequence::from("ATCGGTCGAA"));
/// introns.push(Sequence::from("ATCGGTCGAGCGTGT"));
///
/// let spliced = Sequence::from("ATGGTCTACATAGCTGACAAACAGCACGTAGCATCTCG\
///                                 AGAGGCATATGGTCACATGTTCAAAGTTTGCGCCTAG");
/// assert_eq!(splice_introns(prna, &introns), spliced);
/// ```
pub fn splice_introns(mut pre_rna: Sequence, introns: &Tile) -> Sequence {
    for intr in introns {
        let positions = find_motif(&pre_rna, intr);
        for index in positions {
            pre_rna.chain.drain(index..(index + intr.len()));
        }
    }
    pre_rna
}

/// Calculate the k-mer composition of a sequence, using a set of k-mers obtained from a predefined alphabet
///
/// A `k-mer` is a substring of a genetic string with a length of `k`.
/// A genetic string of length n can be seen as composed of n−k+1 overlapping k-mers.
/// The k-mer composition of a genetic string encodes the number of times that each possible k-mer occurs in the string.
/// It is a generalization of the GC-Content, which can be seen as a 1-mer composition.
///
/// # Arguments
/// * `seq` - genetic sequence for which k-mer composition should be calculated
/// * `k` - length of the `k-kmers`
/// * `alphabet` - slice containing set of nucleobases used that can occur in the given genetic sequence
///
/// # Example
/// ```
/// use biotech::processing::transformers::k_mer_composition;
/// use biotech::ds::sequence::Sequence;
///
/// let seq = Sequence::from("CTTCGAAAGTTTGGGCCGAGTCTTACAGTCGGTCTTGAAGCAAAGTAACGAACTCCACGG");
/// let k = 2; // calculate 2-mer composition
/// let alphabet = [b'A', b'C', b'T', b'G']; // use the DNA alphabet
///
/// assert_eq!(
///     k_mer_composition(&seq, k, &alphabet).unwrap(),
///     Vec::from([7, 4, 0, 5, 3, 2, 4, 5, 2, 5, 5, 2, 4, 2, 5, 4])
/// );
/// ```
pub fn k_mer_composition(seq: &Sequence, k: usize, alphabet: &[u8]) -> Result<Vec<usize>> {
    // k-mer cannot be longer than the sequence itself
    if k > seq.len() {
        return Err(BioError::InvalidInputSize);
    }
    // Generate all possible k-mers from alphabet
    let ltrs = std::str::from_utf8(alphabet).map_err(|_| BioError::TypeConversionError)?;
    let alpha = Alphabet::from_chars_in_str(ltrs).map_err(|_| BioError::TypeConversionError)?;
    let words = alpha.all_words(Some(k)).filter(|x| x.len() == k);
    let mut kmers = Tile::new();
    for w in words {
        kmers.push(Sequence::from(w));
    }
    // Sort k-mers according to ordering from alphabet
    kmers = sort_lexicographically(&kmers, alphabet)?;

    // Calculate k-mer composition
    let mut kmer_composition = vec![];
    for kmer in &kmers {
        let pos = find_motif(seq, kmer);
        kmer_composition.push(pos.len());
    }
    Ok(kmer_composition)
}

/// Generate all k-length substrings (k-mers) of a genetic sequence
///
/// # Arguments
/// * `seq` - top-level supersequence
/// * `k` - length of generated k-mers
///
/// # Example
/// ```
/// use biotech::processing::transformers::generate_k_mers;
/// use biotech::ds::sequence::Sequence;
///
/// let seq = Sequence::from("AGTCGGTCTTG");
/// let k = 9;
///
/// assert_eq!(
///     generate_k_mers(&seq, k).unwrap(),
///     Vec::from([
///         Sequence::from("AGTCGGTCT"),
///         Sequence::from("GTCGGTCTT"),
///         Sequence::from("TCGGTCTTG")
///     ])
/// );
/// ```
pub fn generate_k_mers(seq: &Sequence, k: usize) -> Result<Vec<Sequence>> {
    if k > seq.len() {
        return Err(BioError::InvalidInputSize);
    }
    Ok((0..seq.len() - k + 1)
        .into_iter()
        .map(|i| Sequence::from(&seq[i..k + i]))
        .collect())
}

/// Returns tuples with error corrections that can be applied to a set of reads
///   
/// Genome sequencers use a chemical procedure to obtain reads from provided biomaterial.
/// Due to the error-prone character of this approach, multiple reads of the same region are usually taken.
/// Errors can be corrected by splitting the obtained set into correct and faulty reads first.
/// If a read or its complement appears in the set at least `split_margin` times, it is regarded as correct.
/// Faulty reads are corrected based on their hamming distance to one of the correct reads.
///
/// # Arguments
///
/// * `reads` - tile containing reads
/// * `split_margin` - number of reads/complements required to treat a read as correct
/// * `hd_margin` - hamming distance used for matching faulty reads to correct ones
///
/// # Example
/// ```
/// use biotech::processing::transformers::correct_read_errors;
/// use biotech::ds::sequence::Sequence;
/// use biotech::ds::tile::Tile;
///
/// let mut reads = Tile::new();
/// reads.push(Sequence::from("TTCAT"));
/// reads.push(Sequence::from("TGAAA"));
/// reads.push(Sequence::from("GAGGA"));
/// reads.push(Sequence::from("ATCAA"));
/// reads.push(Sequence::from("TTGAT"));
///
/// correct_read_errors(&reads, 2, 1);
/// ```
pub fn correct_read_errors(
    reads: &Tile,
    split_margin: usize,
    hd_margin: usize,
) -> Result<Vec<(Sequence, Sequence)>> {
    // Count number of repeated reads and complements
    let mut read_counter = HashMap::<Sequence, usize>::new();
    for read in reads {
        // Insert read if not present already and increment count
        *read_counter.entry(read.clone()).or_insert(0) += 1;
        // Handle complement case
        let compl = complement_dna(read.clone());
        if read_counter.contains_key(&compl) {
            *read_counter.get_mut(&compl).ok_or(BioError::ItemNotFound)? += 1;
            *read_counter.get_mut(read).ok_or(BioError::ItemNotFound)? += 1;
        }
    }

    // Split according to split margin
    let mut correct_reads = HashSet::<Sequence>::new();
    let mut faulty_reads = HashSet::<Sequence>::new();
    read_counter.iter().for_each(|(fr, &cnt)| {
        if cnt >= split_margin {
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
            if hamming_distance(fr, cr)? <= hd_margin {
                corrections.push((fr.clone(), cr.clone()));
                break;
            }
            // H(complement(x)) <= hamming_distance_margin
            let complement = complement_dna(cr.clone());
            if hamming_distance(fr, &complement)? == hd_margin {
                corrections.push((fr.clone(), complement));
                break;
            }
        }
    }

    Ok(corrections)
}

/// Sort sequences of varying length lexicographically, according to the order of bases define in an alphabet
///   
/// # Arguments
///
/// * `reads` - tile containing reads
/// * `split_margin` - number of reads/complements required to treat a read as correct
/// * `hd_margin` - hamming distance used for matching faulty reads to correct ones
///
/// # Example
/// ```
/// use biotech::processing::transformers::sort_lexicographically;
/// use biotech::ds::sequence::Sequence;
/// use biotech::ds::tile::Tile;
///
/// let mut reads = Tile::new();
/// reads.push(Sequence::from("TAT"));
/// reads.push(Sequence::from("TGAAA"));
/// reads.push(Sequence::from("GAGGA"));
/// reads.push(Sequence::from("ATCAA"));
/// reads.push(Sequence::from("TTGA"));
/// 
/// // Define custom alphabet
/// let alphabet = [b'G', b'T', b'A', b'C'];
/// 
/// // Tile sorted according to custom lexicographic order
/// let mut result = Tile::new();
/// result.push(Sequence::from("GAGGA"));
/// result.push(Sequence::from("TGAAA"));
/// result.push(Sequence::from("TTGA"));
/// result.push(Sequence::from("TAT"));
/// result.push(Sequence::from("ATCAA"));
/// 
/// // Sort
/// let sorted = sort_lexicographically(&reads, &alphabet).unwrap();
/// assert_eq!(sorted, result);
/// ```
pub fn sort_lexicographically(sequences: &Tile, alphabet: &[u8]) -> Result<Tile> {
    // Build Trie data structure on which search is performed
    let mut trie_builder = Trie::new(alphabet);
    let trie = trie_builder.build(sequences)?;
    // Define recursive function that is used to perform a DFS on the Trie
    fn walk_trie_rec(trie: &Graph<TrieNode, u8>, node_id: u64, sorted: &mut Tile) {
        // Leaf node found -> terminate search path
        if trie.out_neighbors(node_id).count() == 0 {
            return;
        }
        // Iterate over children and recurse
        let node = trie.get_node(&node_id);
        for &c in node.data.children.iter() {
            if c != -1 {
                if trie.get_node(&(c as u64)).data.ending == true {
                    let substrings = &trie.get_node(&(c as u64)).data.substring;
                    substrings
                        .iter()
                        .for_each(|s| sorted.push(Sequence::from(s)));
                }
                walk_trie_rec(trie, c as u64, sorted);
            }
        }
    }
    // Sort
    let mut sorted = Tile::new();
    let root = trie.get_root().ok_or(BioError::ItemNotFound)?;
    walk_trie_rec(trie, root, &mut sorted);
    Ok(sorted)
}
