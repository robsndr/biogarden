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

    #[test]
    fn linguistic_complexity() {
        let input = read_sequence("input/linguistic_complexity.fasta");
        let lc = analysis::seq::linguistic_complexity(&input);
        assert_eq!(lc.unwrap(), 0.9330378);
    }

    #[test]
    fn protein_mass() {
        let input = read_sequence("input/protein_mass.fasta");
        let mass = analysis::spectro::weighted_mass(&input);
        assert_eq!(mass.unwrap(), 114193.58444000047);
    }

    #[test]
    fn n_statistic() {
        let input = read_sequences("input/nxx_stat.fasta");
        let n75 = analysis::stat::n_statistic(&input, 75);
        assert_eq!(n75, 97);
        let n50 = analysis::stat::n_statistic(&input, 50);
        assert_eq!(n50, 125);
    }

    #[test]
    fn expected_restriction_sites() {
        // Inputs
        let recognition_seq = Sequence::from("AGCAAGGTG");
        let n = 888799;
        let gc = Vec::<f64>::from([
            0.000, 0.075, 0.114, 0.199, 0.205, 
            0.270, 0.321, 0.357, 0.422, 0.452, 
            0.523, 0.572, 0.629, 0.692, 0.712, 
            0.761, 0.835, 0.885, 0.922, 1.000
        ]);
        // Outputs
        let expected_result = Vec::<f64>::from([
            0.000, 0.004, 0.021, 0.224, 0.252, 
            0.708, 1.258, 1.721, 2.594, 2.954, 
            3.517, 3.567, 3.239, 2.479, 2.186, 
            1.446, 0.523, 0.165, 0.043, 0.000
        ]);
        // Compute
        let result = analysis::stat::expected_restriction_sites(&recognition_seq, n, &gc);
        assert_eq!(result, expected_result);
    }
}