// #![warn(missing_debug_implementations, missing_docs)]

mod sequence;
mod ds;
mod algo;
mod io;

use io::fasta::Record;
use io::fasta::{FastaRead, Reader};

use ds::tile::Tile;
use sequence::Sequence;
use ndarray::prelude::*;


fn main() {
    
    let mut test_string : String = String::from("TCAATGCATGCGGGTCTATATGCAT");
    let mut test_string2 : String = String::from("CATCGTAATGACGGCCT");

    let mut test_string : Sequence = Sequence::from(test_string);
    let mut test_string2 : Sequence = Sequence::from(test_string2);

    // print!("Count: {:?}\n", algo::count(&test_string));
    
    // test_string = algo::complement_dna(test_string);
    // print!("Complement: {}\n", test_string);
    
    // test_string = algo::transcribe_dna(test_string);
    // print!("Transcribe: {}", test_string);

    // let gc = algo::gc_content(&test_string);
    // print!("GC: {}\n", gc);

    // let rna = algo::transcribe_dna(test_string);
    // print!("RNA: {}\n", rna);
    
    // let protein = algo::tranlate_rna(test_string);
    // print!("{}\n", protein);

    // let mendel = algo::mendel_first_law(18, 19, 23);
    // print!("{}", mendel);

    // let offspring =  algo::expected_offspring(18137, 16426, 18904, 18674, 18160, 18728);
    // print!("{}", offspring);
    
    // print!("KMP: {:?}", algo::knuth_morris_pratt(&test_string, &test_string2));
    
    // print!("Hamming: {:?}", algo::hamming_distance(&test_string, &test_string2));

    // Reader based 
    // let mut reader = Reader::from_file(std::path::Path::new(r"C:\Users\Robert\Desktop\biotech\src\in.fasta")).unwrap();
    // let mut record = Record::new();
    // let mut matrix = Tile::new();
    // loop {
    //     reader
    //     .read(&mut record)
    //     .expect("fasta reader: got an io::Error or could not read_line()");
    //     if record.is_empty() {
    //         break;
    //     }
    //     matrix.push(Sequence::from(record.clone()));
    // } 
    // print!("{}", matrix);
    // let a  = matrix.into_array3();
    // let profile = algo::calc_profile(&a);
    // print!("Profile: {:#?}", profile);
    // print!("{}", Sequence::from(algo::calc_consensus(&profile)));


    // Permutations
    // let mut permutes : Vec<Vec<u32>> = vec![];
    // let mut x = vec![1,2,3,4,5]; 
    // algo::permutations(5, &mut x, &mut permutes);
    // for x in &permutes {
    //     print!("{:?}\n", x);
    // }
    // algo::overlap_graph(&matrix);


    // search_motifs(my_str1, my_str2);
    // println!("{:#?}", search_motifs(my_str1, my_str2));
    // let mut arr = array![[b'A', b'T', b'C', b'C', b'A', b'G', b'C', b'T'],
    //                      [b'G', b'G', b'G', b'C', b'A', b'A', b'C', b'T'],
    //                      [b'A', b'T', b'G', b'G', b'A', b'T', b'C', b'T'],
    //                      [b'A', b'A', b'G', b'C', b'A', b'A', b'C', b'C'],
    //                      [b'T', b'T', b'G', b'G', b'A', b'A', b'C', b'T'],
    //                      [b'A', b'T', b'G', b'C', b'C', b'A', b'T', b'T'],
    //                      [b'A', b'T', b'G', b'G', b'C', b'A', b'C', b'T']];
    // print!("{:#?}", encode_output(&calc_consensus(&calc_profile(&encode_input(&arr)))));


    // Splice RNA
    // let mut reader = Reader::from_file(std::path::Path::new(r"C:\Users\Robert\Desktop\biotech\src\in.fasta")).unwrap();
    // let mut record = Record::new();
    // let mut matrix = Tile::new();
    //
    // .read(&mut record)
    // .expect("fasta reader: got an io::Error or could not read_line()");
    // let mut pre_rna = Sequence::from(record.clone());
    // loop {
    //     reader
    //     .read(&mut record)
    //     .expect("fasta reader: got an io::Error or could not read_line()");
    //     if record.is_empty() {
    //         break;
    //     }
    //     matrix.push(Sequence::from(record.clone()));
    // } 
    // algo::rna_splice(&mut pre_rna, &matrix);
    // pre_rna = algo::transcribe_dna(pre_rna);
    // pre_rna = algo::tranlate_rna(pre_rna);
    // print!("{}", pre_rna);


    // Open Reading Frames
    let mut reader = Reader::from_file(std::path::Path::new(r"C:\Users\Robert\Desktop\biotech\src\in.fasta")).unwrap();
    let mut record = Record::new();
    let mut matrix = Tile::new();
    
    reader
    .read(&mut record)
    .expect("fasta reader: got an io::Error or could not read_line()");
    let mut dna = Sequence::from(record.clone());
    
    let frames = algo::open_reading_frames(&dna);

    for f in &frames {
        print!("{}\n", f);
    }


    // let a  = matrix.into_array3();
    // let profile = algo::calc_profile(&a);
    // print!("Profile: {:#?}", profile);
    // print!("{}", Sequence::from(algo::calc_consensus(&profile)))


}
