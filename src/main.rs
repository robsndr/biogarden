mod dna;
mod rna;
mod ds;
mod protein;
mod algo;
mod io;
use dna::DNA;
use ds::tile;
use ndarray::prelude::*;

use io::fasta::Record;
use io::fasta::{FastaRead, Reader};

fn main() {
    
    let mut my_str1 : String = String::from("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA");
    let mut my_str2 : String = String::from("N{P}[ST]{P}");
    
    // let mut d : DNA = DNA::from(my_str1);
    // print!("{}", &d);
    // print!("{:?}\n", &d.count());
    // d.complement();
    // print!("{}", d);
    // print!("{}", rna::RNA::from(&d));
    // print!("{}\n", d.gc_content());
    // let rna = rna::RNA::from(&d);
    // print!("{}\n", protein::Protein::from(&rna));
    // print!("{}", mendel_first_law(15, 17, 19));
    // print!("{}", expected_offspring(18137, 16426, 18904, 18674, 18160, 18728));
    // print!("{}", fibo_die(6, 3));
    // print!("{}", gc_content(&my_str));
    // print!("{:#?}", find_repeats(&my_str1, &my_str2));
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

    // for ((x, y, z), val) in a.indexed_iter() {
    //     println!("{:?} {:?} {:?}", x, y, z);
    // }

    // A == [1,0,0,0]
    // C == [0,1,0,0]
    // G == [0,0,1,0]
    // T == [0,0,0,1]

    // let mut arr = array![[b'A', b'T', b'C', b'C', b'A', b'G', b'C', b'T'],
    //                      [b'G', b'G', b'G', b'C', b'A', b'A', b'C', b'T'],
    //                      [b'A', b'T', b'G', b'G', b'A', b'T', b'C', b'T'],
    //                      [b'A', b'A', b'G', b'C', b'A', b'A', b'C', b'C'],
    //                      [b'T', b'T', b'G', b'G', b'A', b'A', b'C', b'T'],
    //                      [b'A', b'T', b'G', b'C', b'C', b'A', b'T', b'T'],
    //                      [b'A', b'T', b'G', b'G', b'C', b'A', b'C', b'T']];
    // print!("{:#?}", encode_output(&calc_consensus(&calc_profile(&encode_input(&arr)))));



    // const fasta_file: &'static [u8] = b">id desc";
    let mut reader = Reader::from_file(std::path::Path::new(r"C:\Users\Robert\Desktop\biotech\src\in.fasta")).unwrap();
    let mut record = Record::new();
    let mut matrix = tile::Tile::<dna::DNA>::new();

    // // Check for errors parsing the record
    loop {
        reader
        .read(&mut record)
        .expect("fasta reader: got an io::Error or could not read_line()");
        if record.is_empty() {
            break;
        }
        matrix.push(dna::DNA::from(record.clone()));
    } 
    // print!("{}", matrix);
    // print!("{}", a);
    // print!("{}", tile::Tile::<protein::Protein>::from(a));
    // let a  = matrix.into_array3();
    // let profile = algo::calc_profile(&a);
    // print!("Profile: {:#?}", profile);
    // print!("{}", dna::DNA::from(algo::calc_consensus(&profile)));

    // let mut permutes : Vec<Vec<u32>> = vec![];
    // let mut x = vec![1,2,3,4,5]; 
    // algo::permutations(5, &mut x, &mut permutes);
    // for x in &permutes {
    //     print!("{:?}\n", x);
    // }

    algo::overlap_graph(matrix);
}
