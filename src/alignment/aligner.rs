use std::cmp;
use ndarray::Array2;

use crate::error::{BioError, Result};
use crate::ds::sequence::Sequence;

/// Compute the alignment of genetic strings against each other
/// 
/// Sequence alignment is the process of comparing and detecting similarities between biological sequences. 
/// It imposes a way of arranging sequences to identify regions of similarity, that may be a consequence of evolutionary relationships.
/// 
/// Four operations are performed to align genetic strings: insertion, deletion, substitution, and match.
/// A scoring function can be used to quantify the likelihood of different bases being matched against each other. 
/// Depending on the end-goal and given scenario, different types of alignment can be performed.
/// An hypothetical alignment scenario might look as follows:
/// ```no_compile
/// // No alignment:
/// CGATGCTAGCGTATCGTAGTCTATCGTAC
/// ACGATGCTAGCGTTTCGTATCATCGTA
///
/// // Aligned:
/// -CGATGCTAGCGTATCGTAGTCTATCGTAC
/// ACGATGCTAGCGTTTCGTA-TC-ATCGTA-
/// ```
/// 
/// The implemented aligner supports affine gap penalties, which can differentiate between the cost `a` of starting a gap
/// and the cost of extending of an already started one `b`.
pub struct SequenceAligner {
    // Size of cost and size buffers
    buffer_size: (usize, usize),
    // Cost buffers
    m: Array2::<i32>,
    x: Array2::<i32>,
    y: Array2::<i32>,
    // Buffers with traceback data
    m_trace: Array2::<u8>,
    x_trace: Array2::<u8>,
    y_trace: Array2::<u8>,
}

impl SequenceAligner  {

    /// Construct new aligner object
    pub fn new() -> SequenceAligner {
        let buffer_dim = (1024_usize, 1024_usize);
        SequenceAligner {
            buffer_size: buffer_dim,
            m: Array2::<i32>::zeros(buffer_dim),
            x: Array2::from_elem(buffer_dim, i32::MIN),
            y: Array2::from_elem(buffer_dim, i32::MIN),
            m_trace: Array2::<u8>::zeros(buffer_dim),
            x_trace: Array2::from_elem(buffer_dim, b'I'),
            y_trace: Array2::from_elem(buffer_dim, b'I'),
        }
    }

    /// Compute the global alignment of two genetic strings
    ///
    /// Maximize the cumulative alignment score between two sequences
    /// 
    /// # Arguments
    /// * `seq1`, `seq2` - sequences to align
    /// * `score` - scoring function to be used during alignment 
    /// * `a` - penalty for opening a gap
    /// * `b` - penalty for extending a gap
    ///
    /// # Example
    /// ```
    /// use biotech::alignment::aligner::SequenceAligner;
    /// use biotech::alignment::score::blosum62;
    /// use biotech::ds::sequence::Sequence;
    /// 
    /// let mut aligner = SequenceAligner::new();
    ///
    /// let s1 = Sequence::from("PRTEINS");
    /// let s2 = Sequence::from("PRTWPSEIN");
    /// let (align_score, s1_aligned, s2_aligned) = 
    ///     aligner.global_alignment(&s1, &s2, &blosum62, -11, -1).unwrap();
    /// 
    /// assert_eq!(align_score, 8);
    /// assert_eq!(s1_aligned, Sequence::from("PRT---EINS"));
    /// assert_eq!(s2_aligned, Sequence::from("PRTWPSEIN-"));
    /// ```
    pub fn global_alignment(&mut self, seq1: &Sequence, seq2: &Sequence,
                            score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> Result<(i32, Sequence, Sequence)> {

        if a > 0 || b > 0 {
            return Err(BioError::InvalidArgumentRange);
        }

        // Allocate a larger buffer if sequences cannot fit
        if seq1.len() > self.buffer_size.0 || seq2.len() > self.buffer_size.1 {
            self.resize_buffers(seq1.len()+1, seq2.len()+1);
        }

        // Initialize score buffer boundaries
        // Alongside uppermost row
        let mut row0 = self.m.row_mut(0);
        row0[1] = a;
        for j in 2..(seq2.len() + 1) { row0[j] = row0[j-1] + b; }
        // Alongside uppermost column
        let mut col0 = self.m.column_mut(0);
        col0[1] = a;
        for i in 2..(seq1.len() + 1) { col0[i] = col0[i-1] + b; }

        // Initialize trace buffer boundaries (used for backtracking)
        self.m_trace.column_mut(0).fill(b'X');
        self.m_trace.row_mut(0).fill(b'Y');

        // Fill global alignment score buffer and trace
        self.compute_scores_global(seq1, seq2, score, a, b);
        let alignment_score = self.m[[seq1.len(), seq2.len()]];

        // Backtrack optimal solution
        let mut k = seq1.len();
        let mut l = seq2.len();
        let trace_valid = |x: &usize, y: &usize| -> bool { *x != 0 || *y != 0 };
        let (s1_aligned, s2_aligned) = self.backtrack(seq1, seq2, &mut k, &mut l, &trace_valid);

        Ok((alignment_score, s1_aligned, s2_aligned))
    }

    /// Compute the local alignment of two genetic strings
    ///
    /// Maximize the alignment score across all substrings of two genetic sequences
    /// 
    /// # Arguments
    /// * `seq1`, `seq2` - sequences to align
    /// * `score` - scoring function to be used during alignment 
    /// * `a` - penalty for opening a gap
    /// * `b` - penalty for extending a gap
    ///
    /// # Example
    /// ```
    /// use biotech::alignment::aligner::SequenceAligner;
    /// use biotech::alignment::score::blosum62;
    /// use biotech::ds::sequence::Sequence;
    /// 
    /// let mut aligner = SequenceAligner::new();
    ///
    /// let s1 = Sequence::from("PLEASANTLY");
    /// let s2 = Sequence::from("MEANLY");
    /// let (align_score, s1_aligned, s2_aligned) = 
    ///     aligner.local_alignment(&s1, &s2, &blosum62, -11, -1).unwrap();
    /// 
    /// assert_eq!(align_score, 12);
    /// assert_eq!(s1_aligned, Sequence::from("LEAS"));
    /// assert_eq!(s2_aligned, Sequence::from("MEAN"));
    /// ```
    pub fn local_alignment(&mut self, seq1: &Sequence, seq2: &Sequence,
                            score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> Result<(i32, Sequence, Sequence)> {

        if a > 0 || b > 0 {
            return Err(BioError::InvalidArgumentRange);
        }

        // Allocate a larger buffer if sequences cannot fit
        if seq1.len() > self.buffer_size.0 || seq2.len() > self.buffer_size.1 {
            self.resize_buffers(seq1.len()+1, seq2.len()+1);
        }

        // Initialize score buffer
        self.m.fill(0);

        // Initialize trace buffer boundaries (used for backtracking)
        self.m_trace.column_mut(0).fill(b'X');
        self.m_trace.row_mut(0).fill(b'Y');

        // Fill alignment score buffer and trace
        self.compute_scores_local(seq1, seq2, score, a, b);

        // Compute position and value of maximum alignment score
        let maximum = self.m.indexed_iter()
                            .fold(((0,0),i32::MIN), |max, x| if *x.1 > max.1 { (x.0, *x.1) } else { max });
        let max_pos = maximum.0;
        let align_score = maximum.1;

        // Backtrace to find optimal alignment
        let mut k = max_pos.0;
        let mut l = max_pos.1;
        let trace_valid = |x: &usize, y: &usize| -> bool { (*x != 0 || *y != 0) && self.m[[*x,*y]] > 0 };
        let (s1_aligned, s2_aligned) = self.backtrack(seq1, seq2, &mut k, &mut l, &trace_valid);

        Ok((align_score, s1_aligned, s2_aligned))
    }

    /// Compute the fitting alignment of two genetic strings
    ///
    /// Maximize the alignment score of a string `seq1` against another string `seq2`,
    /// where `seq2` is aligned with a substring `seq1′` of `seq1`.
    /// 
    /// # Arguments
    /// * `seq1`, `seq2` - sequences to align, len(seq1) >= len(seq2)
    /// * `score` - scoring function to be used during alignment 
    /// * `a` - penalty for opening a gap
    /// * `b` - penalty for extending a gap
    ///
    /// # Example
    /// ```
    /// use biotech::alignment::aligner::SequenceAligner;
    /// use biotech::alignment::score::unit;
    /// use biotech::ds::sequence::Sequence;
    /// 
    /// let mut aligner = SequenceAligner::new();
    ///
    /// let s1 = Sequence::from("GCAAACCATAAGCCCTACGTGCCGCCTGTTTAAACTCGCGAACTGAAT\
    ///                          CTTCTGCTTCACGGTGAAAGTACCACAATGGTATCACACCCCAAGGAAAC");
    /// let s2 = Sequence::from("GCCGTCAGGCTGGTGTCCG");
    /// let (align_score, s1_aligned, s2_aligned) = 
    ///     aligner.fitting_alignment(&s1, &s2, &unit, -1, -1).unwrap();
    /// 
    /// assert_eq!(align_score, 5);
    /// assert_eq!(s1_aligned, Sequence::from("GCCCT-A--C-G-TG-CCG"));
    /// assert_eq!(s2_aligned, Sequence::from("GCCGTCAGGCTGGTGTCCG"));
    /// ```
    pub fn fitting_alignment(&mut self, seq1: &Sequence, seq2: &Sequence,
                                score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> Result<(i32, Sequence, Sequence)> {

        if a > 0 || b > 0 {
            return Err(BioError::InvalidArgumentRange);
        }

        if seq1.len() < seq2.len() {
            return Err(BioError::InvalidInputSize);
        }

        // Allocate a larger buffer if sequences cannot fit
        if seq1.len() > self.buffer_size.0 || seq2.len() > self.buffer_size.1 {
            self.resize_buffers(seq1.len()+1, seq2.len()+1);
        }

        // Initialize score buffer
        self.m.fill(0);
        // Initialize score buffer alongside uppermost row
        let mut row0 = self.m.row_mut(0);
        row0[1] = a;
        for j in 2..(seq2.len() + 1) { row0[j] = row0[j-1] + b; }

        // Initialize trace buffer boundaries (used for backtracking)
        self.m_trace.column_mut(0).fill(b'X');
        self.m_trace.row_mut(0).fill(b'Y');

        // Fill alignment score buffer and trace
        self.compute_scores_global(seq1, seq2, score, a, b);

        // Compute position and value of maximum alignment score alongside last column
        let maximum = self.m.column(seq2.len())
                            .indexed_iter()
                            .fold((0,i32::MIN), |max, x| if *x.1 > max.1 { (x.0, *x.1) } else { max });
        let max_pos = (maximum.0, seq2.len());
        let alignment_score = maximum.1;

        // Backtrace to find optimal alignment
        let mut k = max_pos.0;
        let mut l = max_pos.1;
        let trace_valid = |_: &usize, y: &usize| -> bool { *y != 0 };
        let (s1_aligned, s2_aligned) = self.backtrack(seq1, seq2, &mut k, &mut l, &trace_valid);

        Ok((alignment_score, s1_aligned, s2_aligned))
    }

    /// Compute the overlap alignment of two genetic strings
    ///
    /// Maximize the local alignment score over all substrings of `seq1` and `seq2`, 
    /// where a suffix of `seq1` is aligned with a prefix of `seq2`. 
    /// 
    /// # Arguments
    /// * `seq1`, `seq2` - sequences to align
    /// * `score` - scoring function to be used during alignment 
    /// * `a` - penalty for opening a gap
    /// * `b` - penalty for extending a gap
    ///
    /// # Example
    /// ```
    /// use biotech::alignment::aligner::SequenceAligner;
    /// use biotech::alignment::score::unit;
    /// use biotech::ds::sequence::Sequence;
    /// 
    /// let mut aligner = SequenceAligner::new();
    ///
    /// let s1 = Sequence::from("CTAAGGGATTCCGGTAATTAGACAG");
    /// let s2 = Sequence::from("ATAGACCATATGTCAGTGACTGTGTAA");
    /// let (align_score, s1_aligned, s2_aligned) = 
    ///     aligner.overlap_alignment(&s1, &s2, &unit, -2, -2).unwrap();
    /// 
    /// assert_eq!(align_score, 2);
    /// assert_eq!(s1_aligned, Sequence::from("ATTAGAC-AG"));
    /// assert_eq!(s2_aligned, Sequence::from("AT-AGACCAT"));
    /// ```
    pub fn overlap_alignment(&mut self, seq1: &Sequence, seq2: &Sequence,
                                score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> Result<(i32, Sequence, Sequence)> {

        // Allocate a larger buffer if sequences cannot fit
        if seq1.len() > self.buffer_size.0 || seq2.len() > self.buffer_size.1 {
            self.resize_buffers(seq1.len()+1, seq2.len()+1);
        }

        // Initialize score buffer
        self.m.fill(0);
        // Initialize trace buffer boundaries (used for backtracking)
        self.m_trace.column_mut(0).fill(b'X');
        self.m_trace.row_mut(0).fill(b'Y');

        // Fill alignment score buffer and trace
        self.compute_scores_global(seq1, seq2, score, a, b);

        // Compute position and value of maximum alignment score alongside last row
        let maximum = self.m.row(seq1.len())
                            .indexed_iter()
                            .fold((0,i32::MIN), |max, x| if *x.1 >= max.1 { (x.0, *x.1) } else { max });
        let max_pos = (seq1.len(), maximum.0);
        let alignment_score = maximum.1;

        // Backtrace to find optimal alignment
        let mut k = max_pos.0;
        let mut l = max_pos.1;
        let trace_valid = |_: &usize, l: &usize| -> bool { *l != 0 };
        let (s1_aligned, s2_aligned) = self.backtrack(seq1, seq2, &mut k, &mut l, &trace_valid);

        Ok((alignment_score, s1_aligned, s2_aligned))
    }

    /// Compute the semiglobal alignment of two genetic strings
    ///
    /// Maximize the alignment score over `seq1` and `seq2` in a scenario, where any gaps appearing as prefixes 
    /// or suffixes of `seq1` and `seq2` do not contribute to the score of the alignment.
    /// 
    /// # Arguments
    /// * `seq1`, `seq2` - sequences to align
    /// * `score` - scoring function to be used during alignment 
    /// * `a` - penalty for opening a gap
    /// * `b` - penalty for extending a gap
    ///
    /// # Example
    /// ```
    /// use biotech::alignment::aligner::SequenceAligner;
    /// use biotech::alignment::score::unit;
    /// use biotech::ds::sequence::Sequence;
    /// 
    /// let mut aligner = SequenceAligner::new();
    ///
    /// let s1 = Sequence::from("TAGCACTTGGATTCTCGG");
    /// let s2 = Sequence::from("CAGCGTGG");
    /// let (align_score, s1_aligned, s2_aligned) = 
    ///     aligner.semiglobal_alignment(&s1, &s2, &unit, -1, -1).unwrap();
    /// 
    /// assert_eq!(align_score, 4);
    /// assert_eq!(s1_aligned, Sequence::from("TAGCA-CTTGGATTCTCGG"));
    /// assert_eq!(s2_aligned, Sequence::from("---CAGCGTGG--------"));
    /// ```
    pub fn semiglobal_alignment(&mut self, seq1: &Sequence, seq2: &Sequence, 
                                    score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) -> Result<(i32, Sequence, Sequence)> {

        // Allocate a larger buffer if sequences cannot fit
        if seq1.len() > self.buffer_size.0 || seq2.len() > self.buffer_size.1 {
            self.resize_buffers(seq1.len()+1, seq2.len()+1);
        }

        // Initialize score buffer
        self.m.fill(0);
        // Initialize trace buffer boundaries (used for backtracking)
        self.m_trace.column_mut(0).fill(b'X');
        self.m_trace.row_mut(0).fill(b'Y');

        // Fill alignment score buffer and trace
        self.compute_scores_global(seq1, seq2, score, a, b);

        // Compute position and value of maximum alignment score alongside last row
        let mr = self.m.row(seq1.len())
                                .indexed_iter()
                                .fold((0,i32::MIN), |max, x| if *x.1 >= max.1 { (x.0, *x.1) } else { max });
        let max_pos_row = (seq1.len(), mr.0);
        let max_score_row = mr.1;

        // Compute position and value of maximum alignment score alongside last column
        let mc = self.m.column(seq2.len())
                                .indexed_iter()
                                .fold((0,i32::MIN), |max, x| if *x.1 > max.1 { (x.0, *x.1) } else { max });
        let max_pos_col = (mc.0, seq2.len());
        let max_score_col = mc.1;

        // Traceback to find edit alignment
        let mut s1_aligned = Sequence::new();
        let mut s2_aligned = Sequence::new();
        let mut max_pos = (0, 0);
        let mut alignment_score = 0_i32;

        // Traceback gaps at tail of alignment
        if max_score_col > max_score_row {
            max_pos = max_pos_col;
            alignment_score = max_score_col;
            for i in (max_pos_col.0+1..seq1.len()+1).rev() {
                s1_aligned.push(seq1[i-1]);
                s2_aligned.push(b'-');
            }
        }
        else {
            max_pos = max_pos_row;
            alignment_score = max_score_row;
            for i in (max_pos_row.1+1..seq2.len()+1).rev() {
                s1_aligned.push(b'-');
                s2_aligned.push(seq2[i-1]);
            }
        }

        // Traceback overlapping local alignment 
        let mut k = max_pos.0;
        let mut l = max_pos.1;
        let trace_valid = |x: &usize, y: &usize| -> bool { x*y != 0 };
        let (mut t1, mut t2) = self.backtrack(seq1, seq2, &mut k, &mut l, &trace_valid);
        t1.reverse();
        t2.reverse();
        s1_aligned.extend(t1);
        s2_aligned.extend(t2);

        // Traceback prefix gaps of alignment 
        if max_score_col > max_score_row {
            for i in (0..k).rev() {
                s1_aligned.push(seq1[i]);
                s2_aligned.push(b'-');
            }
        }
        else {
            for i in (0..l).rev() {
                s1_aligned.push(b'-');
                s2_aligned.push(seq2[i]);
            }
        }

        // Reverse obtained alignment to get requested order
        s1_aligned.reverse();
        s2_aligned.reverse();

        Ok((alignment_score, s1_aligned, s2_aligned))
    }

    fn compute_scores_global(&mut self, seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) {

        for i in 1..(seq1.len() + 1) {
            for j in 1..(seq2.len() + 1) {

                // Establish optimal action to perform with regards to gaps in x (seq1)
                self.x[[i,j]] = cmp::max(self.m[[i-1,j]] + a, self.x[[i-1,j]].saturating_add(b));
                self.x_trace[[i,j]] = if self.x[[i,j]] == self.m[[i-1,j]] + a { b'M' } else { b'I' };

                // Establish optimal action to perform with regards to gaps in y (seq2)
                self.y[[i,j]] = cmp::max(self.m[[i,j-1]] + a, self.y[[i,j-1]].saturating_add(b));
                self.y_trace[[i,j]] = if self.y[[i,j]] == self.m[[i,j-1]] + a { b'M' } else { b'I' };

                // Find optimal action from: Replace/Match | Gap in X | Gap in Y
                let maximum = cmp::max(self.m[[i-1,j-1]] + score(&seq1[i-1], &seq2[j-1]),
                                       cmp::max(self.x[[i,j]], self.y[[i,j]]));

                // Update main trace depending on optimal action
                if maximum ==  self.y[[i,j]] {
                    self.m_trace[[i,j]] = b'Y';
                }
                else if maximum == self.x[[i,j]] {
                    self.m_trace[[i,j]] = b'X';
                }
                else {
                    self.m_trace[[i,j]] = b'R';
                }

                // Update main score buffer with optimal cost
                self.m[[i,j]] = maximum;
            }
        }
    }

    fn compute_scores_local(&mut self, seq1: &Sequence, seq2: &Sequence, score: &dyn Fn(&u8, &u8) -> i32, a: i32, b: i32) {

        for i in 1..(seq1.len() + 1) {
            for j in 1..(seq2.len() + 1) {

                // Establish optimal action to perform with regards to gaps in x (seq1)
                self.x[[i,j]] = cmp::max( self.m[[i-1,j]] + a, self.x[[i-1,j]].saturating_add(b));
                self.x_trace[[i,j]] = if self.x[[i,j]] == self.m[[i-1,j]] + a { b'M' } else { b'I' };
                // Local alignment -> zero out score if < 0 as alignment can start and end anywhere
                self.x[[i,j]] = if self.x[[i,j]] < 0 { 0 } else { self.x[[i,j]] };

                // Establish optimal action to perform with regards to gaps in y (seq2)
                self.y[[i,j]] = cmp::max( self.m[[i,j-1]] + a, self.y[[i,j-1]].saturating_add(b));
                self.y_trace[[i,j]] = if self.y[[i,j]] == self.m[[i,j-1]] + a { b'M' } else { b'I' };
                // Local alignment -> zero out score if < 0 as alignment can start and end anywhere
                self.y[[i,j]] = if self.y[[i,j]] < 0 { 0 } else { self.y[[i,j]] };

                // Find optimal action from: Replace/Match | Gap in X | Gap in Y
                let maximum = cmp::max( self.m[[i-1,j-1]] + score(&seq1[i-1], &seq2[j-1]),
                                        cmp::max( self.x[[i,j]], self.y[[i,j]])
                                    );

                // Update main trace depending on optimal action
                if maximum ==  self.y[[i,j]] {
                    self.m_trace[[i,j]] = b'Y';
                }
                else if maximum == self.x[[i,j]] {
                    self.m_trace[[i,j]] = b'X';
                }
                else {
                    self.m_trace[[i,j]] = b'R';
                }

                // Update edit distance in memoization table
                // Local alignment -> zero out score if < 0 as alignment can start and end anywhere
                self.m[[i,j]] = if maximum < 0 { 0 } else { maximum };
            }
        }
    }

    fn backtrack(&self, seq1: &Sequence, seq2: &Sequence, k: &mut usize, l: &mut usize,
                    trace_valid: &dyn Fn(&usize, &usize) -> bool) -> (Sequence, Sequence) {

        // backtrack to find optimal alignment based on obtained costs
        let mut s1_aligned = Sequence::new();
        let mut s2_aligned = Sequence::new();
        let mut curtrace : u8 = b'M';

        // Backtrack
        while trace_valid(k, l) {
            match curtrace {

                b'M' => {
                    match self.m_trace[[*k,*l]] {
                        b'R' => {
                            // Replace/Match
                            s1_aligned.push(seq1[*k-1]);
                            s2_aligned.push(seq2[*l-1]);
                            *k -= 1;
                            *l -= 1;
                        }

                        b'X' => {
                            // Insert/Delete in X, Gap in Y
                            curtrace = b'X';
                            s1_aligned.push(seq1[*k-1]);
                            s2_aligned.push(b'-');
                            *k -= 1;
                        }

                        b'Y' => {
                            // Insert/Delete in Y, Gap in X
                            curtrace = b'Y';
                            s1_aligned.push(b'-');
                            s2_aligned.push(seq2[*l-1]);
                            *l -= 1;
                        }

                        _ => {}
                    }

                }

                b'X' => {
                    match self.x_trace[[*k,*l]] {
                        b'M' => {
                            // Replace/Match, Return to M
                            curtrace = b'M';
                        }
                        _ => {
                            // Insert/Delete in X, Gap in Y
                            s1_aligned.push(seq1[*k-1]);
                            s2_aligned.push(b'-');
                            *k -= 1;
                        }
                    }
                }

                b'Y' => {
                    // Delete
                    match self.y_trace[[*k,*l]] {
                        b'M' => {
                            // Replace/Match, Return to M
                            curtrace = b'M';
                        }
                        _ => {
                            // Insert/Delete in Y, Gap in X
                            s1_aligned.push(b'-');
                            s2_aligned.push(seq2[*l-1]);
                            *l -= 1;
                        }
                    }
                }
                _ => {}
            }
        }

        s1_aligned.reverse();
        s2_aligned.reverse();

        (s1_aligned, s2_aligned)
    }

    fn resize_buffers(&mut self, dim1: usize, dim2: usize) {
        self.buffer_size = (dim1, dim2);
        self.m = Array2::zeros(self.buffer_size);
        self.x = Array2::from_elem(self.buffer_size, i32::MIN);
        self.y = Array2::from_elem(self.buffer_size, i32::MIN);
        self.m_trace = Array2::zeros(self.buffer_size);
        self.x_trace = Array2::from_elem(self.buffer_size, b'I');
        self.y_trace = Array2::from_elem(self.buffer_size, b'I');
    }
}

impl Default for SequenceAligner {
    fn default() -> Self {
        Self::new()
    }
}
