// use crate biotech;
use biotech::algo;
use biotech::analysis;
use biotech::io::fasta::*;
use biotech::ds::sequence::{Sequence};
use biotech::ds::tile::{Tile};
use std::collections::HashMap;

#[cfg(test)]
mod integration {
    use super::*;

    fn read_sequences(name: &str) -> Tile {
        let x = format!("./tests/data/{}", name);
        let path = std::path::Path::new(&x);
        let mut reader = Reader::from_file(path).unwrap();
        let mut record = Record::new();
        let mut matrix = Tile::new();
        loop {
            reader
            .read(&mut record)
            .expect("fasta reader: got an io::Error or could not read_line()");
            if record.is_empty() {
                break;
            }
            matrix.push(Sequence::from(record.clone()));
        } 
        matrix
    }

    fn read_sequence(name: &str) -> Sequence {
        let x = read_sequences(name);
        x[0].clone()
    }

    #[test]
    fn count_nucleotides() {
        let input = read_sequence("input/count_nucleotides.fasta");
        let counts = analysis::seq::count_nucleotides(&input);
        assert_eq!(counts, HashMap::<u8, usize>::from([(b'A', 195), (b'C', 217), (b'G', 216), (b'T', 229)]));
    }

    #[test]
    fn gc_content() {
        let matrix = read_sequences("input/gc_content.fasta");
        let mut gcc = 0.0;
        // Get maximum GC content in the set of all input sequences
        for seq in &matrix {
            let temp = analysis::seq::gc_content(&seq);
            if temp > gcc {
                gcc = temp;
            }
        }
        assert_eq!(0.5273311897106109, gcc);
    }

    #[test]
    fn hamming_distance() {
        let input = read_sequences("input/hamming_distance.fasta");
        let hd = analysis::seq::hamming_distance(&input[0], &input[1]);
        assert_eq!(hd.unwrap(), 477);
    }

    #[test]
    fn edit_distance() {
        let input = read_sequences("input/edit_distance.fasta");
        let ed = analysis::seq::edit_distance(&input[0], &input[1]);
        assert_eq!(ed.unwrap(), 299);
    }
    
    #[test]
    fn transitions_transversions() {
        let input = read_sequences("input/transversions.fasta");
        let ratio = analysis::seq::transition_transversion_ratio(&input[0], &input[1]);
        assert_eq!(ratio.unwrap(), 2.032258064516129);
    }
}