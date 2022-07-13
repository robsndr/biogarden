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
    
    // Comparative metrics
    let edit_dist = edit_distance(&a, &b).unwrap();
    let tt_ratio = transition_transversion_ratio(&a, &b).unwrap();
    println!("[A-B] Edit Distance: {}, TT Ratio: {}", edit_dist, tt_ratio);

    // Pattern finding
    let positions_tcg = find_motif(&a, &Sequence::from("TCG"));
    println!("[A] Positons TCG: {:?}", positions_tcg);
    let rev_cs = reverse_complement_substrings(&a, 4, 6);
    println!("[A] Reverse complement substrings: {:?}", rev_cs);
    
    // Pattern based compare
    let lcss = longest_common_subsequence(&a, &b);
    println!("[A-B] Longest common subsequence: {}", lcss);

    // Transcribe    
    let b_rna = transcribe_dna(b);
    println!("[B] RNA: {}", b_rna);

}