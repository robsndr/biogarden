use std::fmt; // Import `fmt`
use super::rna::RNA;
use std::collections::HashMap;

#[derive(Debug)]
pub struct Protein {
    amino: Vec<u8>,
}

impl Protein {
    pub fn new() -> Protein {
        Protein { amino: Vec::<u8>::new() }
    }
}

impl From<String> for Protein {
    fn from(s: String) -> Self {
        Protein { amino: s.into_bytes() }
    }
}

impl From<&RNA> for Protein {
    fn from(rna: &RNA) -> Self {
        // mRNA <-> amino-acid translation table (codon table)
        let codon_table: HashMap<&str, &str> = HashMap::from(
            [ ("UUU", "F"),    ("CUU", "L"),   ("AUU", "I"),   ("GUU", "V"),
                ("UUC", "F"),    ("CUC", "L"),   ("AUC", "I"),   ("GUC", "V"),
                ("UUA", "L"),    ("CUA", "L"),   ("AUA", "I"),   ("GUA", "V"),
                ("UUG", "L"),    ("CUG", "L"),   ("AUG", "M"),   ("GUG", "V"),
                ("UCU", "S"),    ("CCU", "P"),   ("ACU", "T"),   ("GCU", "A"),
                ("UCC", "S"),    ("CCC", "P"),   ("ACC", "T"),   ("GCC", "A"),
                ("UCA", "S"),    ("CCA", "P"),   ("ACA", "T"),   ("GCA", "A"),
                ("UCG", "S"),    ("CCG", "P"),   ("ACG", "T"),   ("GCG", "A"),
                ("UAU", "Y"),    ("CAU", "H"),   ("AAU", "N"),   ("GAU", "D"),
                ("UAC", "Y"),    ("CAC", "H"),   ("AAC", "N"),   ("GAC", "D"),
                ("UAA", "Stop"), ("CAA", "Q"),   ("AAA", "K"),   ("GAA", "E"),
                ("UAG", "Stop"), ("CAG", "Q"),   ("AAG", "K"),   ("GAG", "E"),
                ("UGU", "C"),    ("CGU", "R"),   ("AGU", "S"),   ("GGU", "G"),
                ("UGC", "C"),    ("CGC", "R"),   ("AGC", "S"),   ("GGC", "G"),
                ("UGA", "Stop"), ("CGA", "R"),   ("AGA", "R"),   ("GGA", "G"),
                ("UGG", "W"),    ("CGG", "R"),   ("AGG", "R"),   ("GGG", "G") ]);
        // Container for final result of transcription
        let mut amino_acid: String = String::from("");
        // Run the translation 
        let s = String::from(rna);
        let mut z = s.chars().peekable();
        while z.peek().is_some() {
            let chunk: String = z.by_ref().take(3).collect();
            match codon_table.get(&chunk as &str) {
                Some(value) => {
                    if value == &"Stop"{
                        break;
                    }
                    amino_acid.push_str(value);
                },
                None => println!("Codon not found in codon table.")
            }
        }
        Protein { amino: amino_acid.into_bytes() } 
    }
}

impl fmt::Display for Protein {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let temp: String = self.amino.iter().map(|&c| c as char).collect::<String>();
        write!(f, "----------\nProtein: \n{}\n----------\n", temp)
    }
}


