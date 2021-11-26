use std::io;

// count number of chars in dna sequence
fn count_nuclea(dna : &String) -> [u16; 24] {
    let mut count: [u16; 24] = [0; 24];
    for c in dna.chars() {
        count[(c as usize) - 65] += 1;
    }
    return count;
}

// count number of chars in dna sequence
fn transcribe_rna(dna : &String) -> String {
    dna.replace("T", "U")
}

fn complement_dna(dna : &String) -> String {
    let mut t: String = dna.chars().rev().collect();
    // A <-> T
    t = t.replace("A", "AR");
    t = t.replace("T", "A");
    t = t.replace("AR", "T");
    // G <-> C
    t = t.replace("G", "GR");
    t = t.replace("C", "G");
    t = t.replace("GR", "C");
    return t;
}


fn mendel_first_law(k: u16, m: u16, n: u16 ) -> f32 {
    
    let k_f32 = k as f32;
    let m_f32 = m as f32;
    let n_f32 = n as f32;
    let total : f32 = k_f32 + m_f32 + n_f32;

    // calculate the probability for a homozygous recessive trait
    // heterozygous organisms mating with each other
    let m_p : f32 = m_f32/total * (m_f32 - 1.0 )/(total - 1.0) * 0.25;
    // homozygous organisms mating with each other
    let n_p : f32 = n_f32/total * (n_f32 - 1.0 )/(total - 1.0);
    // heterozygous with homozygous mating
    let mut mn_p : f32 = m_f32/total * n_f32/(total - 1.0) * 0.5;
    mn_p += m_f32/(total - 1.0) * n_f32/total * 0.5;
    
    // probability of a dominant allele is 1 - prob_homozygous
    1.0 - (m_p + n_p + mn_p)
}


fn expected_offspring(x: u16, y: u16, z: u16, q: u16, p: u16, r: u16 ) -> f32 {
    // The probability of AA-AA, AA-Aa, AA-aa couples
    // having dominant phenotype offspring is 1
    let xyz_e : f32 = 2.0 * (x as f32 + y as f32 + z as f32);
    // P[dominant offspring] = 0.75
    let q_e : f32 = 2.0 * q as f32 * 0.75;
    // P[dominant offspring] = 0.5
    let p_e : f32 = 2.0 * p as f32 *  0.5;
    xyz_e + q_e + p_e
}

fn fibo(n: usize, k: u64) -> u64 {
    // F1 == F2 == 1
    // Initialize the vector with 1's 
    let mut memory: Vec<u64> = vec![1; n];
    // Members of (i-1) generation are too young too breed
    // Members of (i-2) generation can breed k new members
    for i in 2..n {
        memory[i] = memory[i-1] + memory[i-2]*k;
    }
    memory[n -1]
}


// fn fibo_die(n: usize, m: usize) -> u64 {
//     // F1 == F2 == 1
//     // Initialize the vector with 1's 
//     let mut memory: Vec<u64> = vec![1; n];
//     for i in 2..n {
//         if i < m {
//             memory[i] = memory[i-1] + memory[i-2];
//         }
//         else{
//             memory[i] = memory[i-1] + memory[i-2] - memory[i-m];
//         }
//     }
//     memory[n -1]
// }

fn gc_content(s: &String) -> f64 {
    let mut gc_count : u32 = 0;
    for c in s.chars() {
        if c == 'G' ||  c == 'C' {
            gc_count += 1;
        }
    }
    gc_count as f64 / s.chars().count() as f64
}

fn hamming_distance(s1: &String, s2: &String) -> Result<u32, io::Error> {
    let both = s1.chars().zip(s2.chars());
    let mut hamming : u32 = 0;
    for pair in both {
        if pair.0 != pair.1 {
            hamming += 1;
        }
    }
    Ok(hamming)
}

fn find_repeats(s: &String, t: &String) -> Vec<usize> {
    s.match_indices(t).map(|(i, _)|i).collect()
}

fn main() {
    // >Rosalind_2484 -0.48629442

    let mut my_str1 : String = String::from("GATATATGCATATACTT");

    let mut my_str2 : String = String::from("ATAT");

    // print!("{}", complement_dna(&my_str));
    
    // print!("{}", mendel_first_law(15, 17, 19));

    // print!("{}", expected_offspring(18137, 16426, 18904, 18674, 18160, 18728));

    // print!("{}", fibo_die(6, 3));

    // print!("{}", gc_content(&my_str));

    // match hamming_distance(&my_str1, &my_str2) {
    //     Ok(value) => {
    //         print!("{}", value);
    //     }
    //     Err(err) => {
    //         println!("Error calculating the hamming distance: {}", err);
    //     }
    // }

    print!("{:#?}", find_repeats(&my_str1, &my_str2));
}
