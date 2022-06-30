// use crate biotech;
use biotech::analysis;
use biotech::ds::sequence::Sequence;
use biotech::ds::tile::Tile;
use biotech::io::fasta::*;
use biotech::processing;
use std::collections::HashMap;
use std::collections::HashSet;

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
        assert_eq!(
            counts,
            HashMap::<u8, usize>::from([(b'A', 195), (b'C', 217), (b'G', 216), (b'T', 229)])
        );
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
            0.000, 0.075, 0.114, 0.199, 0.205, 0.270, 0.321, 0.357, 0.422, 0.452, 0.523, 0.572,
            0.629, 0.692, 0.712, 0.761, 0.835, 0.885, 0.922, 1.000,
        ]);
        // Outputs
        let expected_result = Vec::<f64>::from([
            0.000, 0.004, 0.021, 0.224, 0.252, 0.708, 1.258, 1.721, 2.594, 2.954, 3.517, 3.567,
            3.239, 2.479, 2.186, 1.446, 0.523, 0.165, 0.043, 0.000,
        ]);
        // Compute
        let result = analysis::stat::expected_restriction_sites(&recognition_seq, n, &gc);
        assert_eq!(result, expected_result);
    }

    #[test]
    fn find_motif() {
        let input = read_sequences("input/find_motif.fasta");
        let positions = processing::patterns::find_motif(&input[0], &input[1]);
        let expected_pos = Vec::from([
            22, 65, 72, 138, 177, 254, 261, 268, 315, 348, 368, 454, 470, 489, 496, 550, 557, 601,
            608, 648, 687, 694, 711, 718, 777,
        ]);
        assert_eq!(positions, expected_pos);
    }

    #[test]
    fn longest_common_subsequence() {
        let input = read_sequences("input/longest_common_subseq.fasta");
        let result = read_sequence("output/longest_common_subseq.fasta");
        let lcss = processing::patterns::longest_common_subsequence(&input[0], &input[1]);
        assert_eq!(lcss, result);
    }

    #[test]
    fn shortest_common_supersequence() {
        let input = read_sequences("input/shortest_common_superseq.fasta");
        let result = read_sequence("output/shortest_common_superseq.fasta");
        let scss = processing::patterns::shortest_common_supersequence(&input[0], &input[1]);
        assert_eq!(scss, result);
    }

    #[test]
    fn longest_common_substring() {
        let input = read_sequences("input/longest_common_substring.fasta");
        let result = read_sequence("output/longest_common_substring.fasta");
        let alphabet = HashSet::<u8>::from([b'A', b'C', b'T', b'G']);
        let bound = 0;
        let lcs = processing::patterns::longest_common_substring(&input, &alphabet, bound);
        assert_eq!(lcs.unwrap(), result);
    }
    #[test]
    fn transcribe_dna() {
        let input = read_sequence("input/transcribe_dna.fasta");
        let output = read_sequence("output/transcribe_dna.fasta");
        let rna = processing::transformers::transcribe_dna(input);
        assert_eq!(rna, output);
    }
    #[test]
    fn complement_dna() {
        let input = read_sequence("input/complement_dna.fasta");
        let complement = read_sequence("output/complement_dna.fasta");
        assert_eq!(processing::transformers::complement_dna(input), complement);
    }

    #[test]
    fn translate_rna() {
        let input = read_sequence("input/translate_rna.fasta");
        let complement = read_sequences("output/translate_rna.fasta");
        assert_eq!(
            processing::transformers::translate_rna(input, Some(1)),
            complement
        );
    }

    #[test]
    fn orf() {
        let input = read_sequence("input/orf.fasta");
        let orfs = read_sequences("output/orf.fasta");
        assert_eq!(processing::transformers::open_reading_frames(&input), orfs);
    }
    #[test]
    fn rna_splice() {
        let mut seq = read_sequence("input/rna_splice_seq.fasta");
        let introns = read_sequences("input/rna_splice_introns.fasta");
        // Remove introns from preRNA
        seq = processing::transformers::splice_introns(seq, &introns);
        // Transcribe into RNA
        seq = processing::transformers::transcribe_dna(seq);
        // Translate RNA into protein and compare expected output
        let protein = read_sequences("output/splice_rna.fasta");
        assert_eq!(
            processing::transformers::translate_rna(seq, Some(1)),
            protein
        );
    }

    #[test]
    fn kmer_composition() {
        let seq = read_sequence("input/kmer_composition.fasta");
        let result = Vec::from([
            366, 390, 359, 391, 387, 395, 397, 362, 362, 377, 363, 376, 393, 393, 390, 369, 377,
            381, 406, 389, 390, 395, 386, 396, 361, 395, 422, 359, 389, 386, 407, 346, 389, 400,
            391, 369, 377, 391, 406, 374, 398, 400, 391, 364, 408, 376, 421, 373, 392, 375, 405,
            418, 410, 384, 384, 426, 406, 374, 362, 380, 395, 374, 384, 360, 363, 374, 396, 392,
            393, 401, 357, 397, 409, 395, 400, 378, 375, 405, 359, 399, 375, 390, 404, 370, 386,
            404, 381, 419, 375, 379, 388, 390, 408, 395, 385, 417, 349, 410, 392, 368, 398, 384,
            360, 400, 406, 396, 358, 392, 366, 417, 395, 380, 393, 371, 404, 409, 363, 401, 396,
            370, 418, 360, 406, 404, 379, 404, 390, 362, 377, 353, 381, 384, 404, 400, 417, 407,
            398, 416, 393, 390, 419, 405, 360, 345, 391, 380, 389, 401, 360, 397, 379, 405, 364,
            401, 366, 396, 374, 393, 397, 362, 364, 422, 404, 367, 385, 411, 414, 384, 386, 415,
            387, 359, 368, 388, 376, 380, 393, 368, 382, 397, 379, 389, 398, 403, 425, 378, 380,
            421, 362, 385, 391, 397, 400, 424, 342, 378, 369, 371, 366, 362, 380, 360, 397, 434,
            403, 402, 413, 400, 382, 397, 383, 378, 403, 394, 386, 385, 419, 367, 376, 413, 406,
            356, 399, 410, 393, 396, 410, 425, 401, 355, 347, 368, 367, 383, 411, 397, 398, 388,
            412, 402, 366, 354, 380, 394, 388, 394, 397, 371, 375, 359, 410, 395, 358, 387, 374,
            367,
        ]);
        assert_eq!(
            processing::transformers::k_mer_composition(&seq, 4, &[b'A', b'C', b'G', b'T'])
                .unwrap(),
            result
        );
    }
}
