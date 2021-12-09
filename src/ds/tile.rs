use crate::dna;
use std::fmt; // Import `fmt`
use ndarray::prelude::*;


#[derive(Debug)]
pub struct Tile<T: IntoIterator + fmt::Display> {
    data: Vec<T>,
}

impl<T: IntoIterator + fmt::Display> Tile<T> {
    pub fn new() -> Tile<T> {
        Tile::<T>{ data: Vec::new() }
    }

    pub fn push(&mut self, value: T) {
        self.data.push(value);
    }

    pub fn pop(&mut self) -> Option<T> {
        self.data.pop()
    }
}

impl<T: IntoIterator + fmt::Display> fmt::Display for Tile<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for x in &self.data {
            write!(f,"{}\n", x);
        }
        Ok(())
    }
}

// Tile <-> Array3<u16>
impl<T: IntoIterator + fmt::Display> From<&Tile<T>> for Array2<u16> {
    fn from(tile: &Tile<T>) -> Self {
        let mut converted = Array2::<u16>::zeros((tile.data.len(), tile.data.first().iter().len()));
        for x in tile.data.iter() {
            print!("{}", x);
            // for value in x.into_iter() {
            //     // converted[(i,j)] = value;
            // }
        }  
        converted
    }
}
// Array2<u16> <-> DNA
// impl From<Array1<u16>> for DNA {
//     fn from(arr: Array1<u16>) -> Self {
//         let mut temp = DNA::new();
//         for (i, value) in arr.indexed_iter() {
//                 match value {
//                     0 => temp.nuclea.push(b'A'),
//                     1 => temp.nuclea.push(b'C'),
//                     2 => temp.nuclea.push(b'G'),
//                     3 => temp.nuclea.push(b'T'),
//                     _ => ()
//                 }
//         }
//         temp
//     }
// }