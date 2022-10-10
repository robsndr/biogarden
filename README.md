# BioGarden

BioGarden is a collection of algorithms created as a project to learn about bioinformatics and rust.

### Installation

- `Cargo.toml`

```toml
[dependencies]
biogarden = "0.1.0"
```

- CLI application

```
$ cargo install biogarden
```

### Usage

For simple cases, sequences can be treated directly within the source:

```rust
use biogarden::ds::sequence::Sequence;

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
```

In cases where multiple long sequences are used, reading data from a file might be more practical.
For example, the semiglobal alignment for [`inputs`](https://github.com/symm3try/biogarden/blob/master/tests/data/input/semiglobal_alignment.fasta) 
of ~10000 base pairs can be computed as follows:

```rust
use biogarden::io::fasta::{FastaRead, Reader};
use std::path::Path;

use biogarden::alignment::*;
use biogarden::ds::tile::Tile;
use biogarden::processing::patterns::*;

fn main() {
    // Allocate container for sequences
    let mut tile = Tile::new();

    // Read input
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

    print!("[A-B] Align score: {}", align_score);
    print!("[A] Aligned: {}", a_align);
    print!("[B] Aligned: {}", b_align);

    let alphabet = HashSet::from([b'A', b'C', b'G', b'T']);
    // Minimum number of times the common substring have to occur
    let minimum_frequency = 6; 
    let lcs = longest_common_substring(&tile, &alphabet, minimum_frequency).unwrap();
    println!("[TILE] LCS: {}", lcs);
}
```
