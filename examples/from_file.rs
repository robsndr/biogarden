use biogarden::io::fasta::{FastaRead, Reader};
use std::path::Path;
use std::collections::HashSet;

use biogarden::alignment::*;
use biogarden::ds::tile::Tile;
use biogarden::processing::patterns::*;

fn main() {
    // Allocate container for sequences
    let mut tile = Tile::new();

    // Read FileB
    let mut path = Path::new(&"tests/data/input/semiglobal_alignment.fasta");
    let mut reader = Reader::from_file(path).unwrap();
    reader.read_all(&mut tile);

    // Sequence alignment
    let mut aligner = aligner::SequenceAligner::new();
    let gap_penalty_open = -1;
    let gap_penalty_enlarge = -2;

    // Semiglobal alignment
    let (align_score, a_align, b_align) = aligner
        .semiglobal_alignment(
            &tile[0],
            &tile[1],
            &score::blosum62,
            gap_penalty_open,
            gap_penalty_enlarge,
        )
        .unwrap();

    println!("[A-B] Align score: {}", align_score);
    println!("[A] Aligned: {}", a_align);
    println!("[B] Aligned: {}", b_align);

    let alphabet = HashSet::from([b'A', b'C', b'G', b'T']);
    // minimum number of times the common substring have to occur 
    // in all considered sequences
    let minimum_frequency = 6; 
    let lcs = longest_common_substring(&tile, &alphabet, minimum_frequency).unwrap();
    println!("[TILE] LCS: {}", lcs);
}
