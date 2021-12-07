mod dna;
mod rna;
mod protein;
mod algo;
use dna::DNA;
use ndarray::prelude::*;

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

// Consensus and Profile
pub fn encode_input(arr: &Array2<u8>) -> Array3::<u16> {
    // Tensor with result encoding
    let mut converted = Array3::<u16>::zeros((arr.len_of(Axis(0)), arr.len_of(Axis(1)), 4));
    // Encode 2D input char arry into output tensor
    // A => [1,0,0,0], B => [0,1,0,0]
    // C => [0,0,1,0], D => [0,0,0,1]
    for ((i, j), value) in arr.indexed_iter() {
            match value {
                b'A' => {
                    converted[(i,j,0)] = 1;
                },
                b'C' => {
                    converted[(i,j,1)] = 1
                },
                b'G' => {
                    converted[(i,j,2)] = 1
                },
                b'T' => {
                    converted[(i,j,3)] = 1
                },
                _ => ()
            }
    }
    converted
}

pub fn calc_profile(arr: &Array3<u16>) -> Array2::<u16>  {
    // Squash tensor into 2d array 
    // Sum the occurs of letters different letters
    // Transpose at the end to get expected shape
    arr.sum_axis(Axis(0)).reversed_axes()
}


pub fn calc_consensus(arr: &Array2<u16>) -> Array1<u16> {
    // Allocate container for result
    let mut consensus = Array1::<u16>::zeros((arr.len_of(Axis(1))));
    // Calculate maximum index for every dimension in array
    for (i, ax) in arr.axis_iter(Axis(1)).enumerate() {
        let mut max : u16 = 0;
        let mut index : usize = 0;
        for (j, it) in ax.indexed_iter() {
            if *it > max {
                max = *it;
                index = j;
            }
        }
        consensus[i] = index as u16;
    }
    consensus
}

pub fn encode_output(arr: &Array1<u16>) -> Vec::<char> {
    // Tensor with result encoding
    let mut out : Vec<char> = vec![];

    for (i, value) in arr.indexed_iter() {
            match value {
                0 => {
                    out.push('A');
                },
                1 => {
                    out.push('C');
                },
                2 => {
                    out.push('G');
                },
                3 => {
                    out.push('T');
                },
                _ => ()
            }
    }
    out
}

// pub fn permute(arr: &, res: &Vec<Array1<u8>>) -> {
//     void perm(char* s, int n, int i){
//         if i >= n-1 {
//             print(s);
//         }
//         else {
//           perm(s, n, i+1);
//           for (int j = i+1; j<n; j++){
//             swap(s[i], s[j]);
//             perm(s, n, i+1);
//             swap(s[i], s[j]);
//           }
//         }
//       }
      
//       perm("ABC", 3, 0);
// }

fn main() {
    
    let mut my_str1 : String = String::from("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA");
    let mut my_str2 : String = String::from("N{P}[ST]{P}");
    
    let mut d : DNA = DNA::from(my_str1);
    print!("{}", d);
    print!("{:?}\n", d.count());
    // d.complement();
    print!("{}", d);
    print!("{}", rna::RNA::from(&d));
    print!("{}\n", d.gc_content());
    let rna = rna::RNA::from(&d);
    print!("{}\n", protein::Protein::from(&rna));

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

}
