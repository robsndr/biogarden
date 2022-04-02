#![warn(missing_debug_implementations)]
use ndarray::prelude::*;
use std::hash::{Hash, Hasher};
use std::ops::{Index, IndexMut};
use std::fmt; 
use crate::io::fasta;

#[derive(Debug, Clone)]
pub struct Sequence {
    pub chain: Vec<u8>,
    pub id: Option<String>,
}

impl Sequence {

    pub fn new() -> Sequence {
        Sequence { chain: Vec::<u8>::new(), id: None}
    }

    pub fn push(&mut self, x: u8) -> () {
        self.chain.push(x);
    }

    pub fn len(&self) -> usize {
        self.chain.len()
    }

    pub fn reverse(&mut self) {
        self.chain.reverse();
    }

}

impl Hash for Sequence {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
        self.chain.hash(state);
    }
}

impl Eq for Sequence {}

impl PartialEq for Sequence {
    fn eq(&self, other: &Self) -> bool {
        self.chain == other.chain
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
        (&self.chain).iter()
    }
}


/*** Utility Traits ***/ 
// Iterator Trait
impl<'a> IntoIterator for &'a mut Sequence {

    type Item = <std::slice::Iter<'a, u8> as Iterator>::Item;
    type IntoIter = std::slice::Iter<'a, u8>;

    fn into_iter(self) -> Self::IntoIter {
        (&mut self.chain).iter()
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

impl Index<usize> for Sequence {
    type Output = u8;
    fn index<'a>(&'a self, i: usize) -> &'a u8 {
        &self.chain[i]
    }
}

// impl<'a> Index<usize> for &'a Sequence {
//     type Output = u8;
//     fn index<'b>(&'b self, i: usize) -> &'b u8 {
//         &self.chain[i]
//     }
// }

// impl IndexMut<usize> for Sequence {
//     fn index_mut<'a>(&'a mut self, i: usize) -> &'a mut Sequence {
//         &mut self.data[i]
//     }
// }