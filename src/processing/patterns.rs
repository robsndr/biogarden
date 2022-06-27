use std::collections::HashMap;
use std::collections::HashSet;

use crate::ds::builders::suffix_tree::SuffixTreeBuilder;
use crate::ds::sequence::Sequence;
use crate::ds::tile::Tile;

pub fn knuth_morris_pratt(seq: &Sequence, pat: &Sequence) -> Vec<usize> {
    let seq = seq.to_string().into_bytes();
    let pat = pat.to_string().into_bytes();

    // Build the partial match table
    let mut partial = vec![0];
    for i in 1..pat.len() {
        let mut j = partial[i - 1];
        while j > 0 && pat[j] != pat[i] {
            j = partial[j - 1];
        }
        partial.push(if pat[j] == pat[i] { j + 1 } else { j });
    }

    // Read 'string' to find 'pattern'
    let mut ret = vec![];
    let mut j = 0;

    for (i, &c) in seq.iter().enumerate() {
        while j > 0 && c != pat[j] {
            j = partial[j - 1];
        }
        if c == pat[j] {
            j += 1;
        }
        if j == pat.len() {
            ret.push(i + 1 - j);
            j = partial[j - 1];
        }
    }
    ret
}

pub fn longest_increasing_subsequence(seq: &[u64]) -> Vec<u64> {
    let mut lis = vec![0; seq.len()];
    let mut pointers = vec![0; seq.len()];

    let mut max_idx: usize = 0;
    let mut max_len: usize = 0;

    lis[0] = 1;
    max_idx = 0;

    for i in 1..lis.len() {
        lis[i] = 1;
        pointers[i] = i;

        for j in 0..i {
            // TODO: Make more generic
            // Pass comparator as parameter
            if seq[i] > seq[j] && lis[i] < lis[j] + 1 {
                lis[i] = lis[j] + 1;
                pointers[i] = j;
                if lis[i] > max_len {
                    max_idx = i;
                    max_len = lis[i];
                }
            }
        }
    }

    let mut result: Vec<u64> = vec![];

    while max_len > 0 {
        result.push(seq[max_idx]);
        max_idx = pointers[max_idx];
        max_len -= 1;
    }

    result.reverse();
    result
}

pub fn longest_common_subsequence(seq1: &Sequence, seq2: &Sequence) -> Sequence {
    let mut match_table = vec![vec![0_usize; seq2.len() + 1]; seq1.len() + 1];
    let mut prev_table = vec![vec![(0_usize, 0_usize); seq2.len() + 1]; seq1.len() + 1];

    for i in 1..(seq1.len() + 1) {
        for j in 1..(seq2.len() + 1) {
            if seq1[i - 1] == seq2[j - 1] {
                match_table[i][j] = match_table[i - 1][j - 1] + 1;
                prev_table[i][j] = (i - 1, j - 1);
            } else {
                if match_table[i - 1][j] > match_table[i][j - 1] {
                    match_table[i][j] = match_table[i - 1][j];
                    prev_table[i][j] = (i - 1, j);
                } else {
                    match_table[i][j] = match_table[i][j - 1];
                    prev_table[i][j] = (i, j - 1);
                }
            }
        }
    }

    let mut lcs = Sequence::new();

    let mut i = seq1.len();
    let mut j = seq2.len();

    while match_table[i][j] != 0 {
        let i_next = prev_table[i][j].0;
        let j_next = prev_table[i][j].1;
        if i_next == i - 1 && j_next == j - 1 {
            lcs.push(seq1[i_next]);
        }
        i = i_next;
        j = j_next
    }

    lcs.reverse();
    lcs
}

pub fn longest_common_supersequence(seq1: &Sequence, seq2: &Sequence) -> Sequence {
    let lcs = longest_common_subsequence(seq1, seq2);
    let mut superseq = Vec::<u8>::new();

    let mut s1 = seq1.into_iter();
    let mut s2 = seq2.into_iter();

    for c in lcs {
        loop {
            match s1.next() {
                Some(x) if *x != c => {
                    superseq.push(*x);
                }
                _ => {
                    break;
                }
            }
        }

        loop {
            match s2.next() {
                Some(x) if *x != c => {
                    superseq.push(*x);
                }
                _ => {
                    break;
                }
            }
        }

        superseq.push(c);
    }

    superseq.extend(s1);
    superseq.extend(s2);

    Sequence::from(superseq.as_slice())
}

pub fn longest_common_substring(matrix: &Tile, bound: usize) -> Sequence {
    let alphabet = HashSet::<u8>::from([b'A', b'C', b'T', b'G']);
    let mut suffix_sequence = Sequence::new();

    // Transform set of sequences into one global search string
    // Separate words using any characters that are not members of the alphabet
    let mut separator = 0;
    for a in matrix {
        suffix_sequence.extend(a.clone());
        suffix_sequence.push(separator);
        separator += 1;
        while alphabet.contains(&separator) {
            separator += 1;
        }
    }

    // Build suffix tre using ukonnen's algorithm in O(n)
    let mut ukonnen_builder = SuffixTreeBuilder::new(alphabet);
    let graph = ukonnen_builder.build(&suffix_sequence);

    // DFS
    // Finds the longest common substring
    let mut lcs = Sequence::new();
    let mut max_len = 0;
    let mut stack = Vec::<(u64, Sequence)>::new();
    stack.push((graph.get_root().unwrap(), Sequence::new()));

    while !stack.is_empty() {
        let (cur_node_id, cur_sequence) = stack.pop().unwrap();

        if cur_sequence.len() > max_len {
            max_len = cur_sequence.len();
            lcs = cur_sequence.clone();
        }

        for eid in graph.out_edges(cur_node_id) {
            let rs_i = graph
                .get_node(&graph.get_edge(eid).end)
                .data
                .reachable_suffixes
                .iter();

            if rs_i.clone().filter(|x| **x == 0).count() == 0 && rs_i.sum::<u64>() >= bound as u64 {
                let start = graph.get_edge(eid).data.as_ref().unwrap().suffix_start;
                let stop = graph.get_edge(eid).data.as_ref().unwrap().suffix_stop;

                let mut t = cur_sequence.clone();
                for i in start..stop + 1 {
                    t.push(suffix_sequence[i as usize]);
                }
                stack.push((graph.get_edge(eid).end, t));
            }
        }
    }

    lcs
}

pub fn subsequences(a: &Sequence, b: &Sequence, limit: Option<usize>) -> Vec<Vec<usize>> {
    let mut result = vec![];
    let mut temp = Vec::<usize>::new();
    let a_idx: usize = 0;
    let b_idx: usize = 0;

    pub fn subsequences_recursive(
        a: &Sequence,
        a_idx: usize,
        b: &Sequence,
        b_idx: usize,
        temp: &mut Vec<usize>,
        result: &mut Vec<Vec<usize>>,
        limit: Option<usize>,
    ) {
        if b_idx == b.len() {
            result.push(temp.clone());
            return;
        }
        for i in a_idx..a.len() {
            if limit.is_some() && result.len() == limit.unwrap() {
                return;
            }
            if b[b_idx] == a[i] {
                temp.push(i);
                subsequences_recursive(a, i + 1, b, b_idx + 1, temp, result, limit);
                temp.pop();
            }
        }
    }

    subsequences_recursive(a, a_idx, b, b_idx, &mut temp, &mut result, limit);
    return result;
}

// Find all reverse-palindromes within seq of n <= length <= m
// Return tuples containing position and length of each palindrome O(n^3)?
pub fn reverse_palindromes(seq: &Sequence, n: usize, m: usize) -> Vec<(usize, usize)> {
    let mut palindromes: Vec<(usize, usize)> = vec![];
    let complements = HashMap::from([(b'A', b'T'), (b'T', b'A'), (b'G', b'C'), (b'C', b'G')]);

    // iterate over every offset within the initial string
    for i in 0..seq.len() {
        // iterate over possible lengths of palindromic substrings
        for j in n..(m + 1) {
            // break if potential substring cannot fit
            if i + j > seq.len() {
                break;
            }
            // check if substring with length `j` at offset `i`
            // is a reverse palindrome
            let mut is_palindrome = true;
            for k in 0..j {
                if seq.chain[i + k] != complements[&seq.chain[i + j - 1 - k]] {
                    is_palindrome = false;
                    break;
                }
            }
            // append (offset, length) into result set
            if is_palindrome {
                palindromes.push((i + 1, j));
            }
        }
    }
    palindromes
}
