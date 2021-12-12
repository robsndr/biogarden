use std::fmt; // Import `fmt`
use super::dna::DNA;
use super::io::fasta;

#[derive(Debug)]
pub struct RNA {
    nuclea: Vec<u8>,
    curr: usize,
    id: Option<String>,
}

impl RNA {
    pub fn new() -> RNA {
        RNA { nuclea: Vec::<u8>::new(), curr: 0, id: None }
    }
}

/*** Type-Conversion Traits ***/ 
// String -> RNA
impl From<String> for RNA {
    fn from(s: String) -> Self {
        RNA { nuclea: s.into_bytes(), curr: 0, id: None }
    }
}
// &[u8] <-> RNA
impl From<&[u8]> for RNA {
    fn from(s: &[u8]) -> Self {
        RNA { nuclea: s.to_vec(), curr: 0, id: None }
    }
}
// DNA -> RNA
impl From<&DNA> for RNA {
    fn from(dna: &DNA) -> Self {
        let temp = String::from(dna);
        RNA { nuclea: temp.replace("T", "U").into_bytes(), curr: 0, id: None }
    }
}
// fasta::Record -> DNA
impl From<fasta::Record> for RNA {
    fn from(r: fasta::Record) -> Self {
        RNA { nuclea: r.seq().to_vec(),  curr: 0, id: Some(r.id().to_string()) }
    }
}
// String <- RNA
impl From<&RNA> for String {
    fn from(rna: &RNA) -> Self {
        rna.nuclea.iter().map(|&c| c as char).collect::<String>()
    }
}

/*** Utility Traits ***/
// Iterator Trait
impl Iterator for RNA {
    // We can refer to this type using Self::Item
    type Item = u8;
    
    fn next(&mut self) -> Option<Self::Item> {
        if (self.curr+1) > self.nuclea.len() {
            self.curr = 0;
            return None
        }
        let idx = self.curr;
        self.curr+=1;
        Some(self.nuclea[idx])
    }
}
// ExactSizeIterator : Iterator
impl ExactSizeIterator for RNA{
    fn len(&self) -> usize { 
        self.nuclea.len()
    }
}

/*** Debug Traits ***/
impl fmt::Display for RNA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let temp: String = self.nuclea.iter().map(|&c| c as char).collect::<String>();
        write!(f, "----------\nRNA: \n{}\n----------\n", temp)
    }
}

