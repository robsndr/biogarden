use std::fmt; // Import `fmt`
use ndarray::prelude::*;
use super::io::fasta;

#[derive(Debug, Clone)]
pub struct DNA {
    nuclea: Vec<u8>,
    curr: usize,
    id: Option<String>,
}

impl DNA {

    pub fn new() -> DNA {
        DNA { nuclea: Vec::<u8>::new(), curr: 0, id: None}
    }

    // Complement the DNA string by reversing in the first step.
    // Swap: 'A' <-> 'T' and 'G' <-> 'C'
    pub fn complement(&mut self) -> () {
        let mut t: String = self.nuclea.iter().rev().map(|&c| c as char).collect::<String>();
        // A <-> T
        t = t.replace("A", "X");
        t = t.replace("T", "A");
        t = t.replace("X", "T");
        // G <-> C
        t = t.replace("G", "X");
        t = t.replace("C", "G");
        t = t.replace("X", "C");
        self.nuclea = t.into_bytes();
    }

    // Count number of chars in dna sequence
    // Return array with numbers representing #occur of given char
    // count[0] == count ['A'] and count[23] == count['Z']
    pub fn count(&self) -> [u16; 24] {
        let mut count: [u16; 24] = [0; 24];
        for c in &self.nuclea {
            count[(*c as usize) - 65] += 1;
        }
        return count;
    }

    pub fn gc_content(&self) -> f64 {
        let mut gc_count : u32 = 0;
        for c in &self.nuclea {
            if *c == b'G' ||  *c == b'C' {
                gc_count += 1;
            }
        }
        gc_count as f64 / self.nuclea.len() as f64
    }

    pub fn len(&self) -> usize {
        self.nuclea.len()
    }
}

/*** Type-Conversion Traits ***/ 
// String -> DNA
impl From<String> for DNA {
    fn from(s: String) -> Self {
        DNA { nuclea: s.into_bytes(), curr: 0, id: None }
    }
}
// &[u8] -> DNA
impl From<&[u8]> for DNA {
    fn from(s: &[u8]) -> Self {
        DNA { nuclea: s.to_vec(), curr: 0, id: None }
    }
}
// Array1 -> DNA
impl From<Array1<u8>> for DNA {
    fn from(a: Array1<u8>) -> Self {
        DNA { nuclea: a.to_vec(), curr: 0, id: None }
    }
}
// fasta::Record -> DNA
impl From<fasta::Record> for DNA {
    fn from(r: fasta::Record) -> Self {
        DNA { nuclea: r.seq().to_vec(),  curr: 0, id: Some(r.id().to_string()) }
    }
}
// String <- DNA
impl From<&DNA> for String {
    fn from(dna: &DNA) -> Self {
        dna.nuclea.iter().map(|&c| c as char).collect::<String>()
    }
}

/*** Utility Traits ***/ 
// Iterator Trait
impl Iterator for DNA {
    // We can refer to this type using Self::Item
    type Item = u8;
    
    fn next(&mut self) -> Option<Self::Item> {
        if (self.curr+1) > self.nuclea.len() {
            return None
        }
        let idx = self.curr;
        self.curr+=1;
        Some(self.nuclea[idx])
    }
}
// ExactSizeIterator : Iterator
impl ExactSizeIterator for DNA{
    fn len(&self) -> usize { 
        self.nuclea.len()
    }
}

/*** Debug Traits ***/ 
// Display 
impl fmt::Display for DNA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut temp : String;
        match &self.id {
            Some(p) => temp = p.to_string() + "\n",
            None => temp = "".to_string(),
        }
        temp += std::str::from_utf8(&self.nuclea).unwrap();
        write!(f, "{}", temp)
    }
}