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

// N{P}[ST]{P}
pub fn generate_motifs(target: &Vec<u8>, result: &mut Vec<String>, i: usize, temp: &mut Vec<u8>) {

    if i == target.len() {
        result.push(std::str::from_utf8(temp).unwrap().to_string());
        return;
    }

    let alphabet = String::from("FLSYCWPHQRITMNKVADEG").into_bytes();
    let mut alternatives = Vec::<u8>::new();
    let mut j = i;
    if target[j] == b'[' {
        j += 1;
        while j < target.len() && target[j] != b']' {
            alternatives.push(target[j]);
            j += 1;
        }
        for elem in alternatives {
            temp.push(elem);
            generate_motifs(target, result, j + 1, temp);
            temp.pop();
        }
    }
    else if target[j] == b'{' {
        j += 1;
        for elem in alphabet {
            if elem != target[j] {
                temp.push(elem);
                generate_motifs(target, result, j + 2, temp);
                temp.pop();
            }
        }
        j += 2;
    }
    else{
        temp.push(target[i]);
        generate_motifs(target, result, i + 1, temp);
        temp.pop();
    }
}

fn search_motifs(st: String, pat: String) -> Vec<usize> {

    let st = st.into_bytes();
    let pat = pat.into_bytes();
    let mut res = Vec::<String>::new();
    let mut temp = Vec::<u8>::new();

    generate_motifs(&pat, &mut res, 0, &mut temp);

    let mut pos = Vec::<usize>::new();
    for elem in &res {
        pos.append(&mut algo::knuth_morris_pratt(&st, elem.as_bytes()));    
    }
    pos
}


fn main() {
    
    let mut my_str1 : String = String::from("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA");
    let mut my_str2 : String = String::from("N{P}[ST]{P}");
    
    let mut d : DNA = DNA::from(my_str1);
    print!("{}", &d);
    // print!("{:?}\n", &d.count());
    d.complement();
    print!("{}", d);
    print!("{}", rna::RNA::from(&d));
    print!("{}\n", d.gc_content());
    let rna = rna::RNA::from(&d);
    // print!("{}\n", protein::Protein::from(&rna));

    // for x in d {
    //     print!("{} ", x as char);
    // }

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
    let mut matrix = tile::Tile::<DNA>::new();
    //
    // // Check for errors parsing the record
    loop {
        reader
        .read(&mut record)
        .expect("fasta reader: got an io::Error or could not read_line()");
        if record.is_empty() {
            break;
        }
        matrix.push(dna::DNA::from(record.seq()));
        // print!("{}\n", Array2::<u16>::from(&dna::DNA::from(record.seq())));
    } 

    print!("{}", matrix);
}
