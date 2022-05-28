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

pub fn global_alignment(seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, penalty: &dyn Fn(&usize) -> i32) -> i32  
{
    // Data containers 
    let mut alignement_score = 0_i32;
    let mut memo = vec![vec![0_i32; seq2.len() + 1 ]; seq1.len() + 1];
    let mut action_matrix = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];

    // initialize table
    for i in 1..(seq1.len() + 1) {
        memo[i][0] = memo[i-1][0] + penalty(&0);
    }
    for j in 1..(seq2.len() + 1) {
        memo[0][j] = memo[0][j-1] + penalty(&0);
    }

    // Calculate edit-distance dp table
    for i in 1..seq1.len()+1 {
        for j in 1..seq2.len()+1 {
            if seq1[i-1] == seq2[j-1] {
                memo[i][j] = memo[i-1][j-1] + score(&seq2[j-1], &seq1[i-1]);
                // No operation
                // action_matrix[i][j] = b'N';
            } else {
                let maximum = cmp::max(memo[i-1][j-1] + score(&seq1[i-1], &seq2[j-1]), cmp::max(memo[i][j-1] + penalty(&0), memo[i-1][j] + penalty(&0)));
                // Evaluate whether edit, insert, replace
                // if maximum == memo[i-1][j] - 5 {
                //     action_matrix[i][j] = b'D';
                // } else if  maximum == memo[i][j-1] - 5 {
                //     action_matrix[i][j] = b'I';
                // } else {
                //     action_matrix[i][j] = b'R';
                // }
                // Update edit distance in memoization table
                memo[i][j] = maximum;
            }
        }
    }

    alignement_score = memo[seq1.len()][seq2.len()];
    alignement_score
}
