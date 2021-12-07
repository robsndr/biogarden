use std::fmt; // Import `fmt`
use super::dna::DNA;

#[derive(Debug)]
pub struct RNA {
    nuclea: Vec<u8>,
}

impl RNA {
    pub fn new() -> RNA {
        RNA { nuclea: Vec::<u8>::new() }
    }
}

impl From<String> for RNA {
    fn from(s: String) -> Self {
        RNA { nuclea: s.into_bytes() }
    }
}

impl From<&RNA> for String {
    fn from(rna: &RNA) -> Self {
        rna.nuclea.iter().map(|&c| c as char).collect::<String>()
    }
}

impl From<&DNA> for RNA {
    fn from(dna: &DNA) -> Self {
        let temp = String::from(dna);
        RNA { nuclea: temp.replace("T", "U").into_bytes() }
    }
}

impl fmt::Display for RNA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let temp: String = self.nuclea.iter().map(|&c| c as char).collect::<String>();
        write!(f, "----------\nRNA: \n{}\n----------\n", temp)
    }
}

