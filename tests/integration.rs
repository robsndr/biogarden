// use crate biotech;
use biotech::algo;
use biotech::io::fasta::*;
use biotech::ds::sequence::{Sequence};
use biotech::ds::tile::{Tile};


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

    #[test]
    fn hamming_distance() {
        let input = read_sequences("input/hamming.fasta");
        if input.len() != 2 { panic!("Wrong input size.")};
        let hd = algo::seq::hamming_distance(&input[0], &input[1]);
        assert_eq!(hd, 481);
    }
}