use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::hash_map::Entry;
use std::cmp;

use crate::ds::sequence::Sequence;
use crate::ds::tile::Tile;

/// Obtain tuple containing edit-distance and edit-alignment of two genetic sequences.
///   
/// The Hamming distance provides a way to model point mutations transforming one genetic string into another. 
/// In practice, point mutations include insertions and deletions in addition to replacements only.
/// This can produce genetic sequence that vary in length, and cannot be compared using hamming distance.
/// In such scenarios, a measure of the minimum number of replacements / insertions / deletions between two sequences, is provided by edit distance.
/// The edit distance alignment provides additional information about the type and location where mutations have occurred.
/// 
/// # Arguments
///
/// * `seq1` - first sequence to calculate edit alignment
/// * `seq2` - second sequence to calculate edit alignment
/// 
pub fn edit_distance_alignment(seq1: &Sequence, seq2: &Sequence) -> (Sequence, Sequence, usize) {

    // Data containers 
    let mut edit_distance = 0_usize;
    let mut memo = vec![vec![0_usize; seq2.len() + 1 ]; seq1.len() + 1];
    let mut action_matrix = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];

    // initialize table
    for i in 0..(seq1.len() + 1) {
        memo[i][0] = i;
    }
    for j in 0..(seq2.len() + 1) {
        memo[0][j] = j;
    }

    // Calculate edit-distance dp table
    for i in 1..seq1.len()+1 {
        for j in 1..seq2.len()+1 {
            if seq1[i-1] == seq2[j-1] {
                memo[i][j] = memo[i-1][j-1];
                // No operation
                action_matrix[i][j] = b'N';
            } else {
                let minimum = cmp::min(memo[i-1][j-1], cmp::min(memo[i][j-1], memo[i-1][j]));
                // Evaluate whether edit, insert, replace
                if minimum == memo[i-1][j] {
                    action_matrix[i][j] = b'D';
                } else if  minimum == memo[i][j-1] {
                    action_matrix[i][j] = b'I';
                } else {
                    action_matrix[i][j] = b'R';
                }
                // Update edit distance in memoization table
                memo[i][j] = minimum + 1;
            }
        }
    }

    edit_distance = memo[seq1.len()][seq2.len()];

    // Backtrace to find edit alignment
    let mut s1_aligned = Sequence::new();
    let mut s2_aligned = Sequence::new();
    let mut k = seq1.len();
    let mut l = seq2.len();
    // Backtrack
    while k != 0 || l != 0 {
        match action_matrix[k][l] {
            b'N' | b'R' => {
                // No-op | Replace
                s1_aligned.push(seq1[k-1]);
                s2_aligned.push(seq2[l-1]);
                k -= 1;
                l -= 1;
            }
            b'D' => {
                // Delete
                s1_aligned.push(seq1[k-1]);
                s2_aligned.push(b'-');
                k -= 1;
            }
            b'I' => {
                // Insert
                s1_aligned.push(b'-');
                s2_aligned.push(seq2[l-1]);
                l -= 1;
            }
            _ => ()
        }
    }

    s1_aligned.reverse();
    s2_aligned.reverse();

    (s1_aligned, s2_aligned, edit_distance)
}

pub fn global_alignment(seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> i32  
{
    // Data containers 
    let mut alignement_score = 0_i32;
    let mut m = vec![vec![0_i32; seq2.len() + 1 ]; seq1.len() + 1];
    let mut x = vec![vec![i32::MIN; seq2.len() + 1 ]; seq1.len() + 1];
    let mut y = vec![vec![i32::MIN; seq2.len() + 1 ]; seq1.len() + 1];

    let mut m_trace = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];
    let mut x_trace = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];
    let mut y_trace = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];


    m[1][0] = a;
    m[0][1] = a;
    m_trace[1][0] = b'X';
    m_trace[0][1] = b'Y';
    
    x_trace[1][0] = b'I';
    y_trace[0][1] = b'I';

    for i in 2..(seq1.len() + 1) {
        m[i][0] = m[i-1][0] + b;
        m_trace[i][0] = b'X';
        x_trace[i][0] = b'I';
    }
    for j in 2..(seq2.len() + 1) {
        m[0][j] = m[0][j-1] + b;
        m_trace[0][j] = b'Y';
        y_trace[0][j] = b'I';
    }

    for i in 1..seq1.len()+1 {
        for j in 1..seq2.len()+1 {

            x[i][j] = cmp::max( m[i-1][j] + a, x[i-1][j].saturating_add(b));
            x_trace[i][j] = if x[i][j] == m[i-1][j] + a { b'M' } else { b'I' };


            y[i][j] = cmp::max( m[i][j-1] + a, y[i][j-1].saturating_add(b));
            y_trace[i][j] = if y[i][j] == m[i][j-1] + a { b'M' } else { b'I' };


            let maximum = cmp::max( m[i-1][j-1] + score(&seq1[i-1], &seq2[j-1]), 
                                    cmp::max( x[i][j], y[i][j] )
                                );

            if maximum ==  y[i][j] {
                m_trace[i][j] = b'Y';
            } 
            else if maximum == x[i][j] {
                m_trace[i][j] = b'X';
            } 
            else {
                m_trace[i][j] = b'R';
            }
            
            // Update edit distance in memoization table
            m[i][j] = maximum;
        }
    }

    for i in 0..seq1.len()+1 {
        for j in 0..seq2.len()+1 {
            // if  m_trace[i][j] == 0 {
            print!("{} ", m_trace[i][j] as char);
            // }
        }
        print!("\n");
    }

    print!("\n\n");

    for i in 0..seq1.len()+1 {
        for j in 0..seq2.len()+1 {
            // if  m_trace[i][j] == 0 {
            print!("{} ", x_trace[i][j] as char);
            // }
        }
        print!("\n");
    }

    // Backtrace to find edit alignment
    let mut s1_aligned = Sequence::new();
    let mut s2_aligned = Sequence::new();
    let mut k = seq1.len();
    let mut l = seq2.len();

    let mut curtrace : u8 = b'M';
    // Backtrack
    while k != 0 || l != 0 {
        // print!("TRACE: {}", m_trace[k][l] as char);
        match curtrace {

            b'M' => {

                match m_trace[k][l] {
                    b'R' => {
                        // No-op | Replace
                        s1_aligned.push(seq1[k-1]);
                        s2_aligned.push(seq2[l-1]);
                        k -= 1;
                        l -= 1;
                    }   

                    b'X' => {
                        curtrace = b'X';
                        s1_aligned.push(seq1[k-1]);
                        s2_aligned.push(b'-');
                        k -= 1;
                    }

                    b'Y' => {
                        curtrace = b'Y';
                        s1_aligned.push(b'-');
                        s2_aligned.push(seq2[l-1]);
                        l -= 1;
                    }

                    _ => {
                        ()
                    }
                }
                             
            }


            b'X' => {
                // Delete
                match x_trace[k][l] {
                    b'M' => {
                        curtrace = b'M';
                    }
                    _ => {
                        s1_aligned.push(seq1[k-1]);
                        s2_aligned.push(b'-');
                        k -= 1;
                    }
                }
            }
            b'Y' => {
                // Delete
                match y_trace[k][l] {
                    b'M' => {
                        curtrace = b'M';
                    }
                    _ => {
                        s1_aligned.push(b'-');
                        s2_aligned.push(seq2[l-1]);
                        l -= 1;
                    }
                }
            }
            _ => {
                ()
            }
        }
    }

    s1_aligned.reverse();
    s2_aligned.reverse();


    print!("{}\n", s1_aligned);
    print!("{}\n", s2_aligned);




    alignement_score = m[seq1.len()][seq2.len()];
    alignement_score
}



pub fn local_alignment(seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> (Sequence, Sequence, i32)
{
    // Data containers 
    let mut alignement_score = 0_i32;
    let mut m = vec![vec![0_i32; seq2.len() + 1 ]; seq1.len() + 1];
    let mut x = vec![vec![0; seq2.len() + 1 ]; seq1.len() + 1];
    let mut y = vec![vec![0; seq2.len() + 1 ]; seq1.len() + 1];

    let mut m_trace = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];
    let mut x_trace = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];
    let mut y_trace = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];

    let mut max_pos = (0_usize, 0_usize);

    m[1][0] = 0;
    m[0][1] = 0;
    m_trace[1][0] = b'X';
    m_trace[0][1] = b'Y';
    
    x_trace[1][0] = b'I';
    y_trace[0][1] = b'I';

    for i in 2..(seq1.len() + 1) {
        m[i][0] = 0;
        m_trace[i][0] = b'X';
        x_trace[i][0] = b'I';
    }
    for j in 2..(seq2.len() + 1) {
        m[0][j] = 0;
        m_trace[0][j] = b'Y';
        y_trace[0][j] = b'I';
    }

    for i in 1..seq1.len()+1 {
        for j in 1..seq2.len()+1 {

            
            x[i][j] = cmp::max( m[i-1][j] + a, x[i-1][j].saturating_add(b));
            x_trace[i][j] = if x[i][j] == m[i-1][j] + a { b'M' } else { b'I' };
            x[i][j] = if x[i][j] < 0 { 0 } else { x[i][j] };

            y[i][j] = cmp::max( m[i][j-1] + a, y[i][j-1].saturating_add(b));
            y_trace[i][j] = if y[i][j] == m[i][j-1] + a { b'M' } else { b'I' };
            y[i][j] = if y[i][j] < 0 { 0 } else { y[i][j] };


            let maximum = cmp::max( m[i-1][j-1] + score(&seq1[i-1], &seq2[j-1]), 
                                    cmp::max( x[i][j], y[i][j])
                                );

            if maximum ==  y[i][j] {
                m_trace[i][j] = b'Y';
            } 
            else if maximum == x[i][j] {
                m_trace[i][j] = b'X';
            } 
            else {
                m_trace[i][j] = b'R';
            }
            
            // Update edit distance in memoization table
            m[i][j] = if maximum < 0 { 0 } else { maximum };

            max_pos = if maximum > m[max_pos.0][max_pos.1] { (i, j) } else { max_pos };
        }
    }

    print!("MAXPOS: {}   {} \n", max_pos.0, max_pos.1);
   
    // for i in 1..seq1.len()+1 {
    //     for j in 1..seq2.len()+1 {
    //         print!("{} ", m[i][j]);
    //     }
    //     print!("\n");
    // }

    // Backtrace to find edit alignment
    let mut s1_aligned = Sequence::new();
    let mut s2_aligned = Sequence::new();
    let mut k = max_pos.0;
    let mut l = max_pos.1;

    let mut curtrace : u8 = b'M';
    // Backtrack
    while (k != 0 || l != 0) && m[k][l] > 0 {
        // print!("TRACE: {}", m_trace[k][l] as char);
        match curtrace {

            b'M' => {

                match m_trace[k][l] {
                    b'R' => {
                        // No-op | Replace
                        s1_aligned.push(seq1[k-1]);
                        s2_aligned.push(seq2[l-1]);
                        k -= 1;
                        l -= 1;
                    }   

                    b'X' => {
                        curtrace = b'X';
                        s1_aligned.push(seq1[k-1]);
                        s2_aligned.push(b'-');
                        k -= 1;
                    }

                    b'Y' => {
                        curtrace = b'Y';
                        s1_aligned.push(b'-');
                        s2_aligned.push(seq2[l-1]);
                        l -= 1;
                    }

                    _ => {
                        ()
                    }
                }
                             
            }


            b'X' => {
                // Delete
                match x_trace[k][l] {
                    b'M' => {
                        curtrace = b'M';
                    }
                    _ => {
                        s1_aligned.push(seq1[k-1]);
                        s2_aligned.push(b'-');
                        k -= 1;
                    }
                }
            }
            b'Y' => {
                // Delete
                match y_trace[k][l] {
                    b'M' => {
                        curtrace = b'M';
                    }
                    _ => {
                        s1_aligned.push(b'-');
                        s2_aligned.push(seq2[l-1]);
                        l -= 1;
                    }
                }
            }
            _ => {
                ()
            }
        }
    }

    s1_aligned.reverse();
    s2_aligned.reverse();

    alignement_score = m[max_pos.0][max_pos.1];
    (s1_aligned, s2_aligned, alignement_score)
}