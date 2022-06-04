use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::hash_map::Entry;
use std::cmp;

use crate::ds::sequence::Sequence;
use crate::ds::tile::Tile;

pub struct SequenceAligner {
    // Size of cost and size buffers
    buffer_size: (usize, usize),
    // Cost buffers
    m: Vec<Vec<i32>>,
    x: Vec<Vec<i32>>,
    y: Vec<Vec<i32>>,    
    // Buffers with traceback data
    m_trace: Vec<Vec<u8>>,
    x_trace: Vec<Vec<u8>>,
    y_trace: Vec<Vec<u8>>,
}

impl SequenceAligner  {

    pub fn new() -> SequenceAligner {
        let buffer_dim = (1024_usize, 1024_usize);
        SequenceAligner {
            buffer_size: buffer_dim,
            m: vec![vec![0_i32; buffer_dim.0 + 1 ]; buffer_dim.1 + 1],
            x: vec![vec![i32::MIN; buffer_dim.0 + 1 ]; buffer_dim.1 + 1],
            y: vec![vec![i32::MIN; buffer_dim.0 + 1 ]; buffer_dim.1 + 1],
            m_trace: vec![vec![0_u8; buffer_dim.0 + 1 ]; buffer_dim.1 + 1],
            x_trace: vec![vec![0_u8; buffer_dim.0 + 1 ]; buffer_dim.1 + 1],
            y_trace: vec![vec![0_u8; buffer_dim.0 + 1 ]; buffer_dim.1 + 1],
        }
    }

    pub fn global_alignment(&mut self, seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> (i32, Sequence, Sequence) {

        if seq1.len() > self.buffer_size.0 || seq2.len() > self.buffer_size.1 {
            self.resize_buffers(seq1.len(), seq2.len());
        }
        
        self.m[1][0] = a;
        self.m[0][1] = a;
        self.m_trace[1][0] = b'X';
        self.m_trace[0][1] = b'Y';
        self.x_trace[1][0] = b'I';
        self.y_trace[0][1] = b'I';
    
        for i in 2..(seq1.len() + 1) {
            self.m[i][0] = self.m[i-1][0] + b;
            self.m_trace[i][0] = b'X';
            self.x_trace[i][0] = b'I';
        }
        for j in 2..(seq2.len() + 1) {
            self.m[0][j] = self.m[0][j-1] + b;
            self.m_trace[0][j] = b'Y';
            self.y_trace[0][j] = b'I';
        }
    
        self.compute_scores_global(seq1, seq2, score, a, b);

        let k = seq1.len();
        let l = seq1.len();
        let trace_valid = |k: usize, l: usize| -> bool { k != 0 || l != 0 };

        let (s1_aligned, s2_aligned) = self.backtrack(seq1, seq2, k, l, &trace_valid);

        (self.m[seq1.len()][seq2.len()], s1_aligned, s2_aligned)
    }

    pub fn local_alignment(&mut self, seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> (i32, Sequence, Sequence) {

        if seq1.len() > self.buffer_size.0 || seq2.len() > self.buffer_size.1 {
            self.resize_buffers(seq1.len(), seq2.len());
        }

        self.m[1][0] = 0;
        self.m[0][1] = 0;
        self.m_trace[1][0] = b'X';
        self.m_trace[0][1] = b'Y';
        
        self.x_trace[1][0] = b'I';
        self.y_trace[0][1] = b'I';
    
        for i in 2..(seq1.len() + 1) {
            self.m[i][0] = 0;
            self.m_trace[i][0] = b'X';
            self.x_trace[i][0] = b'I';
        }
        for j in 2..(seq2.len() + 1) {
            self.m[0][j] = 0;
            self.m_trace[0][j] = b'Y';
            self.y_trace[0][j] = b'I';
        }
    

        self.compute_scores_local(seq1, seq2, score, a, b);

        // Copmute position of maximum score
        let mut max_pos = (0_usize, 0_usize);
        let mut maximum = i32::MIN;
        for i in 1..seq1.len()+1 {
            for j in 1..seq2.len()+1 {
                if maximum < self.m[i][j] { 
                    max_pos = (i, j);
                    maximum = self.m[i][j];
                }
            }
        }

        // Backtrace to find 
        let mut k = max_pos.0;
        let mut l = max_pos.1;
        let trace_valid = |k: usize, l: usize| -> bool { (k != 0 || l != 0) && self.m[k][l] > 0 };
        let (s1_aligned, s2_aligned) = self.backtrack(seq1, seq2, k, l, &trace_valid);

        let align_score = self.m[max_pos.0][max_pos.1];
        (align_score, s1_aligned, s2_aligned)

    }

    pub fn fitting_alignment(&mut self, seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> (i32, Sequence, Sequence) {

        // TODO: Refactor into proper Result<> return value
        if seq1.len() < seq2.len() {
            panic!("LEN(SEQ1) !< LEN(SEQ2)");
        }

        if seq1.len() > self.buffer_size.0 || seq2.len() > self.buffer_size.1 {
            self.resize_buffers(seq1.len(), seq2.len());
        }

        // m[1][0] = a;
        self.m[0][1] = a;
        self.m_trace[1][0] = b'X';
        self.m_trace[0][1] = b'Y';
        
        self.x_trace[1][0] = b'I';
        self.y_trace[0][1] = b'I';

        for i in 2..(seq1.len() + 1) {
            // m[i][0] = m[i-1][0] + b;
            self.m_trace[i][0] = b'X';
            self.x_trace[i][0] = b'I';
        }
        for j in 2..(seq2.len() + 1) {
            self.m[0][j] = self.m[0][j-1] + b;
            self.m_trace[0][j] = b'Y';
            self.y_trace[0][j] = b'I';
        }

        self.compute_scores_global(seq1, seq2, score, a, b);

        let mut mx = i32::MIN;
        let mut max_pos = (0, 0);
        for i in 0..seq1.len()+1 {
            if mx < self.m[i][seq2.len()] {
                mx = self.m[i][seq2.len()];
                max_pos = (i, seq2.len());
            }
        }

        let mut k = max_pos.0;
        let mut l = max_pos.1;
        let trace_valid = |k: usize, l: usize| -> bool { l != 0 };
        let (s1_aligned, s2_aligned) = self.backtrack(seq1, seq2, k, l, &trace_valid);

        let align_score = self.m[max_pos.0][max_pos.1];
        (align_score, s1_aligned, s2_aligned)
    }

    fn compute_scores_global(&mut self, seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) {
        for i in 1..seq1.len()+1 {
            for j in 1..seq2.len()+1 {
    
                self.x[i][j] = cmp::max( self.m[i-1][j] + a, self.x[i-1][j].saturating_add(b));
                self.x_trace[i][j] = if self.x[i][j] == self.m[i-1][j] + a { b'M' } else { b'I' };
    
    
                self.y[i][j] = cmp::max( self.m[i][j-1] + a, self.y[i][j-1].saturating_add(b));
                self.y_trace[i][j] = if self.y[i][j] == self.m[i][j-1] + a { b'M' } else { b'I' };
    
    
                let maximum = cmp::max( self.m[i-1][j-1] + score(&seq1[i-1], &seq2[j-1]), 
                                        cmp::max( self.x[i][j], self.y[i][j] )
                                    );
    
                if maximum ==  self.y[i][j] {
                    self.m_trace[i][j] = b'Y';
                } 
                else if maximum == self.x[i][j] {
                    self.m_trace[i][j] = b'X';
                } 
                else {
                    self.m_trace[i][j] = b'R';
                }

                // Update edit distance in memoization table
                self.m[i][j] = maximum;
            }
        }
    }

    fn compute_scores_local(&mut self, seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) {
        for i in 1..seq1.len()+1 {
            for j in 1..seq2.len()+1 {
    
                self.x[i][j] = cmp::max( self.m[i-1][j] + a, self.x[i-1][j].saturating_add(b));
                self.x_trace[i][j] = if self.x[i][j] == self.m[i-1][j] + a { b'M' } else { b'I' };
                self.x[i][j] = if self.x[i][j] < 0 { 0 } else { self.x[i][j] };
    
                self.y[i][j] = cmp::max( self.m[i][j-1] + a, self.y[i][j-1].saturating_add(b));
                self.y_trace[i][j] = if self.y[i][j] == self.m[i][j-1] + a { b'M' } else { b'I' };
                self.y[i][j] = if self.y[i][j] < 0 { 0 } else { self.y[i][j] };
    
    
                let maximum = cmp::max( self.m[i-1][j-1] + score(&seq1[i-1], &seq2[j-1]), 
                                        cmp::max( self.x[i][j], self.y[i][j])
                                    );
    
                if maximum ==  self.y[i][j] {
                    self.m_trace[i][j] = b'Y';
                } 
                else if maximum == self.x[i][j] {
                    self.m_trace[i][j] = b'X';
                } 
                else {
                    self.m_trace[i][j] = b'R';
                }
                
                // Update edit distance in memoization table
                self.m[i][j] = if maximum < 0 { 0 } else { maximum };    
            }
        }
    }

    fn backtrack(&self, seq1: &Sequence, seq2: &Sequence, mut k: usize, mut l: usize, 
                    trace_valid: &dyn Fn(usize, usize) -> bool) -> (Sequence, Sequence) {

        // backtrack to find optimal alignment based on obtained costs
        let mut s1_aligned = Sequence::new();
        let mut s2_aligned = Sequence::new();
        let mut curtrace : u8 = b'M';

        // Backtrack
        while trace_valid(k, l) {
            match curtrace {

                b'M' => {
                    match self.m_trace[k][l] {
                        b'R' => {
                            // Replace/Match
                            s1_aligned.push(seq1[k-1]);
                            s2_aligned.push(seq2[l-1]);
                            k -= 1;
                            l -= 1;
                        }   

                        b'X' => {
                            // Insert/Delete in X, Gap in Y
                            curtrace = b'X';
                            s1_aligned.push(seq1[k-1]);
                            s2_aligned.push(b'-');
                            k -= 1;
                        }

                        b'Y' => {
                            // Insert/Delete in Y, Gap in X
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
                    match self.x_trace[k][l] {
                        b'M' => {
                            // Replace/Match, Return to M
                            curtrace = b'M';
                        }
                        _ => {
                            // Insert/Delete in X, Gap in Y
                            s1_aligned.push(seq1[k-1]);
                            s2_aligned.push(b'-');
                            k -= 1;
                        }
                    }
                }

                b'Y' => {
                    // Delete
                    match self.y_trace[k][l] {
                        b'M' => {
                            // Replace/Match, Return to M
                            curtrace = b'M';
                        }
                        _ => {
                            // Insert/Delete in Y, Gap in X
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

        (s1_aligned, s2_aligned)
    }

    fn resize_buffers(&mut self, dim1: usize, dim2: usize) {
        self.buffer_size = (dim1, dim2);
        self.m = vec![vec![0_i32; dim2 + 1 ]; dim1 + 1];
        self.x = vec![vec![i32::MIN; dim2 + 1 ]; dim1 + 1];
        self.y = vec![vec![i32::MIN; dim2 + 1 ]; dim1 + 1];
        self.m_trace = vec![vec![0_u8; dim2 + 1 ]; dim1 + 1];
        self.x_trace = vec![vec![0_u8; dim2 + 1 ]; dim1 + 1];
        self.y_trace = vec![vec![0_u8; dim2 + 1 ]; dim1 + 1];
    }
}


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
pub fn edit_distance_alignment(seq1: &Sequence, seq2: &Sequence) -> (Sequence, Sequence, u128) {

    // Data containers 
    let mut edit_distance = 0_u128;
    let mut memo = vec![vec![0_u128; seq2.len() + 1 ]; seq1.len() + 1];
    let mut count = vec![vec![0_u128; seq2.len() + 1 ]; seq1.len() + 1];
    let mut action_matrix = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];

    // initialize table
    for i in 0..(seq1.len() + 1) {
        memo[i][0] = i as u128;
        count[i][0] = 1;
    }
    for j in 0..(seq2.len() + 1) {
        memo[0][j] = j as u128;
        count[0][j] = 1;
    }

    // Calculate edit-distance dp table
    for i in 1..seq1.len()+1 {
        for j in 1..seq2.len()+1 {
            let minimum = cmp::min(memo[i-1][j-1] + ((seq1[i-1] != seq2[j-1]) as u128), cmp::min(memo[i][j-1] + 1, memo[i-1][j] + 1));
            // Evaluate whether edit, insert, replace
            if minimum == memo[i-1][j-1] + ((seq1[i-1] != seq2[j-1]) as u128) {
                action_matrix[i][j] = b'R';
                count[i][j] += count[i-1][j-1] % 134_217_727;
            }
            if minimum == memo[i-1][j] + 1 {
                action_matrix[i][j] = b'D';
                count[i][j] += count[i-1][j] % 134_217_727;
            } 
            if  minimum == memo[i][j-1] + 1 {
                action_matrix[i][j] = b'I';
                count[i][j] += count[i][j-1] % 134_217_727 ;
            } 
            
            // Update edit distance in memoization table
            memo[i][j] = minimum % 134_217_727;
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

    print!("COUNT: {} ", count[seq1.len()][seq2.len()]);

    (s1_aligned, s2_aligned, edit_distance)
}


pub fn overlap_alignment(seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> i32
{
    // Data containers 
    let mut alignement_score = 0_i32;
    let mut m = vec![vec![0_i32; seq2.len() + 1 ]; seq1.len() + 1];
    let mut x = vec![vec![0; seq2.len() + 1 ]; seq1.len() + 1];
    let mut y = vec![vec![0; seq2.len() + 1 ]; seq1.len() + 1];

    let mut m_trace = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];
    let mut x_trace = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];
    let mut y_trace = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];

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

            y[i][j] = cmp::max( m[i][j-1] + a, y[i][j-1].saturating_add(b));
            y_trace[i][j] = if y[i][j] == m[i][j-1] + a { b'M' } else { b'I' };

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
            m[i][j] = maximum;
        }
    }

    
    let mut mx = i32::MIN;
    let mut max_pos = (0, 0);
    
    for j in 0..seq2.len()+1 {
        if mx <= m[seq1.len()][j] {
            mx = m[seq1.len()][j];
            max_pos = (seq1.len(), j);
        }
    }

    // Backtrace to find edit alignment
    let mut s1_aligned = Sequence::new();
    let mut s2_aligned = Sequence::new();
    let mut k = max_pos.0;
    let mut l = max_pos.1;

    let mut curtrace : u8 = b'M';
    // Backtrack
    while l != 0 {
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
    alignement_score = m[max_pos.0 as usize][max_pos.1 as usize];
    
    print!("SCORE: {} \n", alignement_score);
    print!("{}\n", s1_aligned);
    print!("{}\n", s2_aligned);

    alignement_score
}

pub fn semiglobal_alignment(seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> i32
{
    // Data containers 
    let mut max_pos = (0, 0);
    let mut alignement_score = 0_i32;
    let mut m = vec![vec![0_i32; seq2.len() + 1 ]; seq1.len() + 1];
    let mut x = vec![vec![0; seq2.len() + 1 ]; seq1.len() + 1];
    let mut y = vec![vec![0; seq2.len() + 1 ]; seq1.len() + 1];

    let mut m_trace = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];
    let mut x_trace = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];
    let mut y_trace = vec![vec![0_u8; seq2.len() + 1 ]; seq1.len() + 1];

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
            // x[i][j] = if x[i][j] < 0 { 0 } else { x[i][j] };

            y[i][j] = cmp::max( m[i][j-1] + a, y[i][j-1].saturating_add(b));
            y_trace[i][j] = if y[i][j] == m[i][j-1] + a { b'M' } else { b'I' };
            // y[i][j] = if y[i][j] < 0 { 0 } else { y[i][j] };


            let maximum = cmp::max( m[i-1][j-1] + score(&seq1[i-1], &seq2[j-1]), 
                                    cmp::max( x[i][j], y[i][j])
                                );


            if maximum == m[i-1][j-1] + score(&seq1[i-1], &seq2[j-1]) {
                m_trace[i][j] = b'R';
            }
            else if maximum ==  y[i][j] {
                m_trace[i][j] = b'Y';
            } 
            else {
                m_trace[i][j] = b'X';
            } 
            
            // Update edit distance in memoization table
            m[i][j] = if maximum < 0 { maximum } else { maximum };

            max_pos = if maximum > m[max_pos.0][max_pos.1] { (i, j) } else { max_pos };
        }
    }

    
    let mut max_col = i32::MIN;
    let mut max_row = i32::MIN;
    let mut max_pos_col = (0, 0);
    let mut max_pos_row = (0, 0);
    let mut max_pos = (0, 0);
    
    for j in 0..seq2.len()+1 {
        if max_row <= m[seq1.len()][j] {
            max_row = m[seq1.len()][j];
            max_pos_row = (seq1.len(), j);
        }
    }

    for i in 0..seq1.len()+1 {
        if max_col <= m[i][seq2.len()] {
            max_col = m[i][seq2.len()];
            max_pos_col = (i, seq2.len());
        }
    }


    // Backtrace to find edit alignment
    let mut s1_aligned = Sequence::new();
    let mut s2_aligned = Sequence::new();


    if max_col > max_row {
        max_pos = max_pos_col;
        for i in (max_pos_col.0+1..seq1.len()+1).rev() {
            s1_aligned.push(seq1[i-1]);
            s2_aligned.push(b'-');
        }
    }
    else {
        max_pos = max_pos_row;
        for i in (max_pos_row.1+1..seq2.len()+1).rev() {
            s1_aligned.push(b'-');
            s2_aligned.push(seq2[i-1]);
        }
    }

    let mut k = max_pos.0;
    let mut l = max_pos.1;

    let mut curtrace : u8 = b'M';
    // Backtrack
    while l*k != 0 {
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

    if max_col > max_row {
        max_pos = max_pos_col;
        for i in (0..k).rev() {
            s1_aligned.push(seq1[i]);
            s2_aligned.push(b'-');
        }
    }
    else {
        max_pos = max_pos_row;
        for i in (0..l).rev() {
            s1_aligned.push(b'-');
            s2_aligned.push(seq2[i]);
        }
    }

    s1_aligned.reverse();
    s2_aligned.reverse();
    alignement_score = m[max_pos.0 as usize][max_pos.1 as usize];
    
    print!("SCORE: {} \n", alignement_score);
    print!("{}\n", s1_aligned);
    print!("{}\n", s2_aligned);

    alignement_score
}

pub fn counting_alignments(seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> i32  
{

    // TODO: Refactor into proper Result<> return value
    if seq1.len() < seq2.len() {
        panic!("LEN(SEQ1) !< LEN(SEQ2)");
    }

    // Data containers 
    let mut alignement_score = 0_i32;
    let mut m = vec![vec![0_i32; seq2.len() + 1 ]; seq1.len() + 1];
    let mut x = vec![vec![i32::MIN; seq2.len() + 1 ]; seq1.len() + 1];
    let mut y = vec![vec![i32::MIN; seq2.len() + 1 ]; seq1.len() + 1];

    let mut scorings = vec![vec![0; seq2.len() ]; seq1.len()];

    m[1][0] = a;
    m[0][1] = a;

    for i in 2..(seq1.len() + 1) {
        m[i][0] = m[i-1][0] + b;
    }
    for j in 2..(seq2.len() + 1) {
        m[0][j] = m[0][j-1] + b;
    }

    for i in 1..seq1.len()+1 {
        for j in 1..seq2.len()+1 {

            x[i][j] = cmp::max( m[i-1][j] + a, x[i-1][j].saturating_add(b));

            y[i][j] = cmp::max( m[i][j-1] + a, y[i][j-1].saturating_add(b));

            m[i][j] = cmp::max( m[i-1][j-1] + score(&seq1[i-1], &seq2[j-1]), 
                                    cmp::max( x[i][j], y[i][j] )
                                );
        }
    }

    for i in 0..seq1.len() {
        for j in 0..seq2.len() {
            
        }
    }

    alignement_score = m[seq1.len()][seq2.len()];
    alignement_score
}