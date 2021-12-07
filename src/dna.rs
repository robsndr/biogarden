use std::fmt; // Import `fmt`

#[derive(Debug)]
pub struct DNA {
    nuclea: Vec<u8>,
}

impl DNA {
    pub fn new() -> DNA {
        DNA { nuclea: Vec::<u8>::new() }
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
}

impl From<String> for DNA {
    fn from(s: String) -> Self {
        DNA { nuclea: s.into_bytes() }
    }
}

impl From<&DNA> for String {
    fn from(dna: &DNA) -> Self {
        dna.nuclea.iter().map(|&c| c as char).collect::<String>()
    }
}

impl fmt::Display for DNA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let temp = String::from(self);
        write!(f, "----------\nDNA: \n{}\n----------\n", temp)
    }
}

