use biogarden::ds::sequence::Sequence;
use biogarden::ds::tile::Tile;

use biogarden::analysis::seq::*;
use biogarden::processing::patterns::*;
use biogarden::processing::transformers::*;

fn main() {
    
    let a = Sequence::from("TTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT");
    let b = Sequence::from("TCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT");

    // Get some properties for sequence A
    let gc_a = gc_content(&a);
    let lc_a = linguistic_complexity(&a).unwrap();
    println!("[A] GC Content: {}, Linguistic complexity: {}", gc_a, lc_a);
    
    // Get some properties for sequence B
    let gc_b = gc_content(&b);
    let lc_b = linguistic_complexity(&b).unwrap();
    println!("[B] GC Content: {}, Linguistic complexity: {}", gc_b, lc_b);

    // Comparative metrics
    let edit_dist = edit_distance(&a, &b).unwrap();
    let tt_ratio = transition_transversion_ratio(&a, &b).unwrap();
    println!("[A-B] Edit Distance: {}, TT Ratio: {}", edit_dist, tt_ratio);

    // Pattern finding
    let positions_tcg = find_motif(&a, &Sequence::from("TCG"));
    println!("[A] Positons ATA: {:?}", positions_tcg);
    let rev_cs = reverse_complement_substrings(&a, 4, 6);
    println!("[A] Reverse complement substrings: {:?}", rev_cs);
    
    // Pattern based compare
    let lcss = longest_common_subsequence(&a, &b);
    println!("[A-B] Longest common subsequence: {:?}", rev_cs);

    // Translation/Transcription
    let a_comp = complement_dna(a);
    println!("[A] Complement: {}", a_comp);
    
    let a_rna = transcribe_dna(a_comp);
    println!("[A] RNA: {}", a_rna);
    
    let a_prot = translate_rna(a_rna, Some(1));
    println!("[A] Prot: {}", a_prot);
}