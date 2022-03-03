// #![warn(missing_debug_implementations, missing_docs)]

mod sequence;
mod ds;
mod algo;
mod io;

use io::fasta::Record;
use io::fasta::{FastaRead, Reader};

use ds::tile::Tile;
use std::fs::File;
use std::io::prelude::*;
use sequence::Sequence;
use ndarray::prelude::*;


fn main() {

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
    
    // reader
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
    // pre_rna = algo::rna_splice(pre_rna, &matrix);
    // pre_rna = algo::transcribe_dna(pre_rna);
    // // TODO: Refactor translate_rna() to take a bound on the number of proteins decoded
    // // TODO: Make translate_rna() return like Variant<Sequence, Vec<Sequence>>
    // let a = algo::translate_rna(pre_rna); 
    // print!("{}", a.first().unwrap());


    // Open Reading Frames
    // let mut reader = Reader::from_file(std::path::Path::new(r"C:\Users\Robert\Desktop\biotech\src\in.fasta")).unwrap();
    // let mut record = Record::new();
    // let mut matrix = Tile::new();
    
    // reader
    // .read(&mut record)
    // .expect("fasta reader: got an io::Error or could not read_line()");
    // let mut dna = Sequence::from(record.clone());
    
    // let frames = algo::open_reading_frames(&dna);

    // for f in &frames {
    //     print!("{}\n", f);
    // }



    // Infer RNA 
    // let mut reader = Reader::from_file(std::path::Path::new(r"C:\Users\Robert\Desktop\biotech\src\in.fasta")).unwrap();
    // let mut record = Record::new();
    // let mut matrix = Tile::new();
    
    // reader
    // .read(&mut record)
    // .expect("fasta reader: got an io::Error or could not read_line()");
    // let protein = Sequence::from(record.clone());
    
    // let num = algo::infer_number_rna(&protein);

    // print!("{}", num);



    // Infer RNA 
    let mut reader = Reader::from_file(std::path::Path::new(r"C:\Users\Robert\Desktop\biotech\src\in.fasta")).unwrap();
    let mut record = Record::new();
    let mut matrix = Tile::new();
    
    reader
    .read(&mut record)
    .expect("fasta reader: got an io::Error or could not read_line()");
    let protein = Sequence::from(record.clone());
    
    // let num : f64 = algo::weighted_mass(&protein);

    // print!("{}", num);

    // let a  = matrix.into_array3();
    // let profile = algo::calc_profile(&a);
    // print!("{}", Sequence::from(algo::calc_consensus(&profile)))


    // CAGCATGGTATCACAGCAGAG

    let seq = Sequence::from("A");
    let pat = Sequence::from("CAGCATGGTATCACAGCAGAG");

    let x  = algo::knuth_morris_pratt(&seq, &protein);


    let mut f = File::create("output.txt").expect("Unable to create file");                                                                                                          
    
    for i in &x{            
        write!(f, "{} ", i);                                                                                                                                                                                                                                                                              
    }   
}
