use std::collections::HashMap;
use std::collections::HashSet;

use crate::ds::builders::suffix_tree::SuffixTreeBuilder;
use crate::ds::sequence::Sequence;
use crate::ds::tile::Tile;
use crate::error::{BioError, Result};

/// Return all positions at which a predefined motif occurs in a genetic sequence
///
/// Complexity: O(n)
///
/// It is a common task in biology to search the genome of an organism for a pattern of biological importance.
/// Such patterns are denoted as `motifs` and may encode specific traits or even disorders.
///  
/// # Arguments
/// * `seq` - genetic sequence to search for motifs
/// * `pat` - the pattern / motif that has to be found
///
/// # Example
/// ```
/// use biogarden::processing::patterns::find_motif;
/// use biogarden::ds::sequence::Sequence;
///
/// let genome = Sequence::from("GATATATGCATATACTT");
/// let motif = Sequence::from("ATAT");
///
/// let pos = find_motif(&genome, &motif);
/// assert_eq!(pos, [1, 3, 9]);
/// ```
pub fn find_motif(seq: &Sequence, pat: &Sequence) -> Vec<usize> {
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

    for (i, &c) in seq.into_iter().enumerate() {
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

/// Return the longest shared non-contiguous motif of genes in two genetic strings.
///
/// Complexity: O(mn)
///
/// Coding regions of DNA are very often interleaved by introns, which are intervals of
/// DNA that do not influence the translation into protein. Introns do not affect the expression of a gene.
///  
/// # Arguments
/// * `seq1`, `seq2` - sequences for which the longest common non-contiguous motif has to be found
///
/// # Example
/// ```
/// use biogarden::processing::patterns::longest_common_subsequence;
/// use biogarden::ds::sequence::Sequence;
///
/// let a = Sequence::from("AACCTTGG");
/// let b = Sequence::from("ACACTGTGA");
///
/// let lcss = longest_common_subsequence(&a, &b);
/// assert_eq!(lcss, Sequence::from("ACCTGG"));
/// ```
pub fn longest_common_subsequence(seq1: &Sequence, seq2: &Sequence) -> Sequence {
    let mut match_table = vec![vec![0_usize; seq2.len() + 1]; seq1.len() + 1];
    let mut prev_table = vec![vec![(0_usize, 0_usize); seq2.len() + 1]; seq1.len() + 1];

    for i in 1..(seq1.len() + 1) {
        for j in 1..(seq2.len() + 1) {
            if seq1[i - 1] == seq2[j - 1] {
                match_table[i][j] = match_table[i - 1][j - 1] + 1;
                prev_table[i][j] = (i - 1, j - 1);
            } else if match_table[i - 1][j] > match_table[i][j - 1] {
                    match_table[i][j] = match_table[i - 1][j];
                    prev_table[i][j] = (i - 1, j);
            } else {
                    match_table[i][j] = match_table[i][j - 1];
                    prev_table[i][j] = (i, j - 1);
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

/// Return the longest increasing subsequence of genes occurring in a genetic strings
///
/// Complexity: O(mn)
///
/// During mutations, genes taken from two different organisms can be moved around in the course of evolution.
/// Based on the initial ordering of the genes in a chromosome, subsequences that are contained
/// in permutations are of interest.
///  
/// # Arguments
/// * `seq` - genetic sequence to search for longest increasing subsequence
///
/// # Example
/// ```
/// use biogarden::processing::patterns::longest_increasing_subsequence;
/// use biogarden::ds::sequence::Sequence;
///
/// let seq = Vec::from([5, 1, 4, 2, 5, 3, 9]);
///
/// let pos = longest_increasing_subsequence(&seq);
/// assert_eq!(pos, [1, 4, 5, 9]);
/// ```
pub fn longest_increasing_subsequence(seq: &[u64]) -> Vec<u64> {
    // TODO: specify ordering operator such that
    // the function can be applied to actual genomes
    let mut lis = vec![0; seq.len()];
    let mut pointers = vec![0; seq.len()];

    let mut max_idx = 0;
    let mut max_len = 0;

    lis[0] = 1;
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

    // Backtrack
    let mut result: Vec<u64> = vec![];
    while max_len > 0 {
        result.push(seq[max_idx]);
        max_idx = pointers[max_idx];
        max_len -= 1;
    }

    result.reverse();
    result
}

/// Return the shortest super-sequence that can be obtained by interleaving two motifs
///
/// Complexity: O(mn)
///  
/// # Arguments
/// * `seq1`, `seq2` - motifs to be interleaved into a top-level super-sequence
///
/// # Example
/// ```
/// use biogarden::processing::patterns::shortest_common_supersequence;
/// use biogarden::ds::sequence::Sequence;
///
/// let a = Sequence::from("TGCATA");
/// let b = Sequence::from("ATCTGAT");
///
/// assert_eq!(shortest_common_supersequence(&a, &b), Sequence::from("ATGCATGAT"));
/// ```
pub fn shortest_common_supersequence(seq1: &Sequence, seq2: &Sequence) -> Sequence {
    let lcs = longest_common_subsequence(seq1, seq2);
    let mut superseq = Vec::<u8>::new();

    let mut s1 = seq1.into_iter();
    let mut s2 = seq2.into_iter();

    for c in lcs {
        loop {
            match s1.next() {
                Some(&x) if x != c => {
                    superseq.push(x);
                }
                _ => {
                    break;
                }
            }
        }

        loop {
            match s2.next() {
                Some(&x) if x != c => {
                    superseq.push(x);
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

    Sequence::from(superseq)
}

/// Return the longest common substring of multiple genetic sequences
///
/// Complexity: O(mn), with `m` - number of sequences, `n` - length of longest sequence
///
/// The motif shared by multiple organisms is often not known in advance.
/// Therefore, it might be required to find the longest interval of shared genes.
/// The implemented function performs the search using a suffix tree, that is built and traversed in `O(mn)`.
///  
/// # Arguments
/// * `tile` - container holding `m` genetic sequences
/// 
/// * `alphabet` - alphabet of `asci` characters that should be considered 
///                as part of the sequences in `tile`, eg. [`A`, `C`, `T`, `G`] for DNA.
///                Characters that are not part of the alphabet will be considered as inter-string separators.
///                
/// * `bound` - **minimal** number of times the given substring has to occur in total.
///             (eg. it might occur more then once in a sequence, or we might search for longest repeat in single genome)
///
/// # Example
/// ```
/// use biogarden::processing::patterns::longest_common_substring;
/// use biogarden::ds::tile::Tile;
/// use biogarden::ds::sequence::Sequence;
/// use std::collections::HashSet;
/// 
/// let mut sequences = Tile::new();
/// sequences.push(Sequence::from("GATTACA"));
/// sequences.push(Sequence::from("TAGACCA"));
/// sequences.push(Sequence::from("ATACA"));
/// 
/// let alphabet = HashSet::<u8>::from([b'A', b'C', b'T', b'G']);
/// let bound = 0;
/// 
/// let lcs = longest_common_substring(&sequences, &alphabet, bound);
/// assert_eq!(lcs.unwrap(), Sequence::from("CA"));
/// ```
pub fn longest_common_substring(tile: &Tile, alphabet: &HashSet::<u8>, bound: usize) -> Result<Sequence> {
    // let alphabet = HashSet::<u8>::from([b'A', b'C', b'T', b'G']);
    let mut suffix_sequence = Sequence::new();

    // Transform set of sequences into one global search string (for Ukkonnen's algo)
    // Separate words using any characters that are not members of the alphabet
    let mut separator = 0;
    for a in tile {
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

    // Initialize data for LCS search
    let mut lcs = Sequence::new();
    let mut max_len = 0;
    let mut stack = vec![(
        graph.get_root().ok_or(BioError::ItemNotFound)?,
        Sequence::new(),
    )];

    // DFS
    while !stack.is_empty() {

        // Fetch next node and suffix candidate from stack
        let (cur_node_id, cur_sequence) = stack.pop().ok_or(BioError::ItemNotFound)?;

        // Check if new candidate for longest substring is found
        if cur_sequence.len() > max_len {
            max_len = cur_sequence.len();
            lcs = cur_sequence.clone();
        }
        // Process outgoing neighbour edges for current node
        for eid in graph.out_edges(cur_node_id) {
            // Get iterator over all suffixes that are reachable from currend node over `eid` edge
            let rs_i = graph.get_node(&graph.get_edge(eid).end).data.reachable_suffixes.iter();
            // Consider only nodes that can reach leaf nodes for suffixes of all strings that are searched
            if rs_i.clone().filter(|&&x| x == 0).count() == 0 && rs_i.sum::<u64>() >= bound as u64 {
                // Id of start node of an edge
                let start = graph.get_edge(eid).data.as_ref().ok_or(BioError::ItemNotFound)?.suffix_start;
                // Id of start node of an edge
                let stop = graph.get_edge(eid).data.as_ref().ok_or(BioError::ItemNotFound)?.suffix_stop;
                // Infer sequence that is enceded for start <= idx <= stop
                let mut t = cur_sequence.clone();
                for i in start..stop + 1 {
                    t.push(suffix_sequence[i as usize]);
                }
                // Put discovered node, and corresponding substring on stack
                stack.push((graph.get_edge(eid).end, t));
            }
        }
    }
    Ok(lcs)
}

/// Get start-positions and lengths of all reverse-complement substrings within a dna strand
///
/// When an enzyme disarms a virus gene, it cuts the virus DNA at spots denoted as restriction sites.
/// This process is most effective, when the two target strands appear directly across from each other along the viral DNA.
/// Something that occurs precisely, when the target is equal to its own reverse complement.
///  
/// # Arguments
/// * `dna` - sequence to search for substrings that are equivalent to their reverse complement
/// * `min_len` - substring length lower bound
/// * `max_len` - substring length higher bound
/// 
/// # Example
/// ```
/// use biogarden::processing::patterns::reverse_complement_substrings;
/// use biogarden::ds::sequence::Sequence;
///
/// let a = Sequence::from("TCAATGCATGCGGGTCTATATGCAT");
/// let min_bound = 4;
/// let max_bound = 12;
///
/// let reverses = [(3, 6), (4, 4), (5, 6), (6, 4), (16, 4), (17, 4), (19, 6), (20, 4)];
/// assert_eq!(reverse_complement_substrings(&a, min_bound, max_bound), reverses);
/// ```
pub fn reverse_complement_substrings(dna: &Sequence, min_len: usize, max_len: usize) -> Vec<(usize, usize)> {
    let mut palindromes: Vec<(usize, usize)> = vec![];
    let complements = HashMap::from([(b'A', b'T'), (b'T', b'A'), (b'G', b'C'), (b'C', b'G')]);

    // iterate over every offset within the initial string
    for i in 0..dna.len() {
        // iterate over possible lengths of palindromic substrings
        for j in min_len..(max_len + 1) {
            // break if potential substring cannot fit
            if i + j > dna.len() {
                break;
            }
            // check if substring with length `j` at offset `i`
            // is a reverse palindrome
            let mut is_palindrome = true;
            for k in 0..j {
                if dna.chain[i + k] != complements[&dna.chain[i + j - 1 - k]] {
                    is_palindrome = false;
                    break;
                }
            }
            // append (offset, length) into result set
            if is_palindrome {
                palindromes.push((i, j));
            }
        }
    }
    palindromes
}

/// Get positions where bases from one genetic string occur in another genome as a subsequence
///
/// Within a DNA string, regions of significance are often scattered and interleaved by introns.
/// It is therefore useful to obtain the positions, where relevant nuclea can be found. 
/// 
/// # Arguments
/// * `a` - top-level sequence to search 
/// * `b` - subsequence to find in top-level sequence
/// * `limit` - maximum number of position sets to generate
/// 
/// # Example
/// ```
/// use biogarden::processing::patterns::subsequences;
/// use biogarden::ds::sequence::Sequence;
///
/// let a = Sequence::from("ACGTACGTGACG");
/// let b = Sequence::from("GTA");
/// let limit = Some(4);
///
/// let reverses : Vec<Vec<usize>> = vec![
///     vec![2, 3, 4],
///     vec![2, 3, 9],
///     vec![2, 7, 9],
///     vec![6, 7, 9]
/// ];
/// assert_eq!(subsequences(&a, &b, limit), reverses);
/// ```
pub fn subsequences(a: &Sequence, b: &Sequence, limit: Option<usize>) -> Vec<Vec<usize>> {
    // TODO: Refactor into Iterator, such that limit is not needed
    let mut result = vec![];
    let mut temp = Vec::<usize>::new();
    let a_idx: usize = 0;
    let b_idx: usize = 0;

    pub fn subsequences_recursive(
        a: &Sequence,
        a_idx: usize,
        b: &Sequence,
        b_idx: usize,
        limit: Option<usize>,
        temp: &mut Vec<usize>,
        result: &mut Vec<Vec<usize>>,
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
                subsequences_recursive(a, i + 1, b, b_idx + 1, limit, temp, result);
                temp.pop();
            }
        }
    }

    subsequences_recursive(a, a_idx, b, b_idx, limit, &mut temp, &mut result);
    result
}
