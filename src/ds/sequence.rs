#![warn(missing_debug_implementations)]
use ndarray::prelude::*;
use std::hash::{Hash, Hasher};
use crate::io::fasta;
use std::ops::Add;
use std::fmt; 

#[derive(Debug, Clone)]
/// Data structure representing genetic sequences
pub struct Sequence {
    pub chain: Vec<u8>,
    pub id: Option<String>,
}

impl Sequence {

    pub fn new() -> Sequence {
        Sequence { chain: Vec::<u8>::new(), id: None}
    }

    pub fn push(&mut self, x: u8) {
        self.chain.push(x);
    }

    pub fn pop(&mut self) -> Option<u8> {
        self.chain.pop()
    }

    pub fn extend(&mut self, b: Sequence) {
        self.chain.extend(b.chain);
    }

    pub fn back(&self) -> Option<&u8> {
        self.chain.last()
    }

    pub fn len(&self) -> usize {
        self.chain.len()
    }

    pub fn is_empty(&self) -> bool {
        self.chain.is_empty()
    }

    pub fn reverse(&mut self) {
        self.chain.reverse();
    }

    pub fn starts_with(&self, prefix: &Sequence) -> bool {

        if prefix.len() > self.chain.len() {
            return false;
        }

        let mut flag = true;
        for (index, x) in prefix.into_iter().enumerate() {
            if *x != self.chain[index] {
                flag = false;
                break;
            }   
        }
        flag
    }

    pub fn ends_with(&self, suffix: &Sequence) -> bool {

        if suffix.len() > self.chain.len() {
            return false;
        }

        let chain_len = self.chain.len();
        let mut flag = true;
        for (index, x) in suffix.into_iter().enumerate() {
            if *x != self.chain[chain_len - index - 1] {
                flag = false;
                break;
            }   
        }
        flag
    }
}

impl Default for Sequence {
    fn default() -> Self {
        Self::new()
    }
}

impl Add for Sequence {
    
    type Output = Self;

    fn add(self, other: Self) -> Self {

        let temp = [self.chain, other.chain].concat();

        Self {
            id: other.id,
            chain: temp
        }
    }
}

impl Hash for Sequence {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // self.id.hash(state);
        self.chain.hash(state);
    }
}

impl Eq for Sequence {}

impl PartialEq for Sequence {
    fn eq(&self, other: &Self) -> bool {
        self.chain == other.chain
    }
}

impl<Idx> std::ops::Index<Idx> for Sequence
where
    Idx: std::slice::SliceIndex<[u8]>,
{
    type Output = Idx::Output;

    fn index(&self, index: Idx) -> &Self::Output {
        &self.chain[index]
    }
}

impl<'a,Idx> std::ops::Index<Idx> for &'a Sequence
where
    Idx: std::slice::SliceIndex<[u8]>,
{
    type Output = Idx::Output;

    fn index(&self, index: Idx) -> &Self::Output {
        &self.chain[index]
    }
}

/*** Type-Conversion Traits ***/ 
// String -> Sequence
impl From<String> for Sequence {
    fn from(s: String) -> Self {
        Sequence { chain: s.into_bytes(), id: None }
    }
}

/*** Type-Conversion Traits ***/ 
// String -> Sequence
impl From<&str> for Sequence {
    fn from(s: &str) -> Self {
        Sequence { chain: Vec::from(s.as_bytes()), id: None }
    }

}

// &[u8] -> Sequence
impl From<&[u8]> for Sequence {
    fn from(s: &[u8]) -> Self {
        Sequence { chain: s.to_vec(), id: None }
    }
}


// &[u8] -> Sequence
impl From<&Vec<u8>> for Sequence {
    fn from(s: &Vec<u8>) -> Self {
        Sequence { chain: s.clone(), id: None }
    }
}

// [u8] -> Sequence
impl From<Vec<u8>> for Sequence {
    fn from(s: Vec<u8>) -> Self {
        Sequence { chain: s, id: None }
    }
}

// Array1 -> Sequence
impl From<Array1<u8>> for Sequence {
    fn from(a: Array1<u8>) -> Self {
        Sequence { chain: a.to_vec(), id: None }
    }
}
// fasta::Record -> Sequence
impl From<fasta::Record> for Sequence {
    fn from(r: fasta::Record) -> Self {
        Sequence { chain: r.seq().to_vec(), id: Some(r.id().to_string()) }
    }
}
// String <- Sequence
impl From<&Sequence> for String {
    fn from(seq: &Sequence) -> Self {
        seq.chain.iter().map(|&c| c as char).collect::<String>()
    }
}

impl From<Sequence> for String {
    fn from(seq: Sequence) -> Self {
        seq.chain.iter().map(|&c| c as char).collect::<String>()
    }
}

/*** Utility Traits ***/ 
// Iterator Trait
impl<'a> IntoIterator for &'a Sequence {

    type Item = <std::slice::Iter<'a, u8> as Iterator>::Item;
    type IntoIter = std::slice::Iter<'a, u8>;

    fn into_iter(self) -> Self::IntoIter {
        self.chain.iter()
    }
}

/*** Utility Traits ***/ 
// Iterator Trait
impl<'a> IntoIterator for &'a mut Sequence {

    type Item = <std::slice::IterMut<'a, u8> as Iterator>::Item;
    type IntoIter = std::slice::IterMut<'a, u8>;

    fn into_iter(self) -> Self::IntoIter {
        self.chain.iter_mut()
    }
}

impl IntoIterator for Sequence {
    type Item = u8;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.chain.into_iter()
    }
}


impl FromIterator<u8> for Sequence {
    fn from_iter<I: IntoIterator<Item=u8>>(iter: I) -> Self {
        let mut s = Sequence::new();
        for i in iter {
            s.push(i);
        }
        s
    }
}

/*** Debug Traits ***/ 
// Display 
impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut temp : String = String::from("");
        temp += std::str::from_utf8(&self.chain).unwrap();
        write!(f, "{}", temp)
    }
}
