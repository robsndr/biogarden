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
    
    let mut test_string : String = String::from("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA");
    // let mut my_str2 : String = String::from("N{P}[ST]{P}");
    
    
    let mut test_string : Sequence = Sequence::from(test_string);
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

    // print!("{}", algo::mendel_first_law(18, 19, 23));

    // print!("{}", algo::expected_offspring(18137, 16426, 18904, 18674, 18160, 18728));
    
    // TODO:
    // print!("{:#?}", knuth_morris_pratt(my_str1, my_str2));
    // match hamming_distance(&my_str1, &my_str2) {
    //     Ok(value) => {
    //         print!("{}", value);
    //     }
    //     Err(err) => {
    //         println!("Error calculating the hamming distance: {}", err);
    //     }
    // }
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

    let mut reader = Reader::from_file(std::path::Path::new(r"C:\Users\Robert\Desktop\biotech\src\in.fasta")).unwrap();
    let mut record = Record::new();
    let mut matrix = Tile::new();

    // // // Check for errors parsing the record
    loop {
        reader
        .read(&mut record)
        .expect("fasta reader: got an io::Error or could not read_line()");
        if record.is_empty() {
            break;
        }
        matrix.push(Sequence::from(record.clone()));
    } 
    print!("{}", matrix);
    let a  = matrix.into_array3();
    let profile = algo::calc_profile(&a);
    print!("Profile: {:#?}", profile);
    print!("{}", Sequence::from(algo::calc_consensus(&profile)));

    // let mut permutes : Vec<Vec<u32>> = vec![];
    // let mut x = vec![1,2,3,4,5]; 
    // algo::permutations(5, &mut x, &mut permutes);
    // for x in &permutes {
    //     print!("{:?}\n", x);
    // }

    // algo::overlap_graph(&matrix);
}
