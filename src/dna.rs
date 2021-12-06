use std::fmt; // Import `fmt`

#[derive(Debug)]
pub struct DNA {
    nuclea: Vec<u8>,
}


impl DNA {
    pub fn new() -> DNA {
        DNA { nuclea: Vec::<u8>::new() }
    }

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
}

impl From<String> for DNA {
    fn from(s: String) -> Self {
        DNA { nuclea: s.into_bytes() }
    }
}

impl fmt::Display for DNA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let temp :String = self.nuclea.iter().map(|&c| c as char).collect::<String>();
        write!(f, "DNA: \n{}\n", temp)
    }
}


