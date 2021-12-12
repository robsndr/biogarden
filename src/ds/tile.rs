use crate::dna;
use std::fmt; // Import `fmt`
use ndarray::prelude::*;


#[derive(Debug)]
<<<<<<< HEAD
pub struct Tile<T> 
    where T:  Iterator<Item=u8> + ExactSizeIterator + fmt::Display + From<String> 
{
    data: Vec<T>,
    curr: usize
}

impl<T> Tile<T>
    where T:  Iterator<Item=u8> + ExactSizeIterator + fmt::Display + From<String>  
{
    pub fn new() -> Tile<T> {
        Tile::<T>{ data: Vec::new(), curr: 0 }
=======
pub struct Tile<T: ExactSizeIterator + Iterator + fmt::Display> {
    data: Vec<T>,
}

impl<T:ExactSizeIterator + Iterator + fmt::Display> Tile<T> {
    pub fn new() -> Tile<T> {
        Tile::<T>{ data: Vec::new() }
>>>>>>> f2d6b30b5dd36c6aab6feff01d6f2848ac118c25
    }

    pub fn push(&mut self, value: T) {
        self.data.push(value);
    }

    pub fn pop(&mut self) -> Option<T> {
        self.data.pop()
    }
<<<<<<< HEAD

    // Convert tile into Array2 
    // TODO: move into From method
    pub fn into_array2(self) -> Array2<u16> {
        let mut converted = Array2::<u16>::zeros((self.data.len(), self.data.first().unwrap().len()));
        for (i,x) in self.data.into_iter().enumerate() {
            for (j, value) in x.into_iter().enumerate() {
                converted[(i,j)] = value as u16;
            }
        }  
        converted
    }

    // Convert tile into Array2 
    // TODO: refactor into From method
    // TODO: remove constant value used to code 3rd dimension of Array3
    // The 3rd dimension is used to encode the character in alphabet
    // Should be reduced by providing a custom alphabet in DNA/RNA/Protein
    pub fn into_array3(self) -> Array3<u16> {
        let mut converted = Array3::<u16>::zeros((self.data.len(), self.data.first().unwrap().len(), 120));
        for (i,x) in self.data.into_iter().enumerate() {
            for (j, value) in x.into_iter().enumerate() {
                converted[(i,j, value as usize)] = 1;
            }
        }  
        converted
    }
}

// Array2<u16> <-> DNA
impl<T> From<Array2<u16>> for Tile<T> 
    where T:  Iterator<Item=u8> + ExactSizeIterator + fmt::Display + From<String> 
{
    fn from(arr: Array2<u16>) -> Self {
        let mut temp = Tile::<T>::new();
        for (i, ax) in arr.axis_iter(Axis(0)).enumerate() {
            let mut s : String = "".to_string();
            for (j, value) in ax.indexed_iter() {
                s.push((*value as u8) as char)
            }
            temp.push(T::from(s));
        }
        temp
    }
}

// Array3<u16> <-> DNA
impl<T> From<Array3<u16>> for Tile<T> 
    where T:  Iterator<Item=u8> + ExactSizeIterator + fmt::Display + From<String> 
{
    fn from(arr: Array3<u16>) -> Self {
        let mut temp = Tile::<T>::new();
        for (i, ax1) in arr.axis_iter(Axis(0)).enumerate() {
            let mut s : String = "".to_string();
            for (j, ax2) in ax1.axis_iter(Axis(0)).enumerate() {
                for (k, val) in ax2.indexed_iter() {
                    if *val == 1 {
                        s.push((k as u8) as char);
                        break;
                    }
                }
            }
            temp.push(T::from(s));
        }
        temp
    }
}

/*** Utility Traits ***/ 
// Iterator Trait
// TODO: Refactor into streaming iterator such that Clone not needed
// Note: Streaming iterator is currently not supported -> check rust::futures
impl<T> Iterator for Tile<T> 
    where T: Clone + Iterator<Item=u8> + ExactSizeIterator + fmt::Display + From<String> 
{
    // We can refer to this type using Self::Item
    type Item = T;
    
    fn next(&mut self) -> Option<Self::Item> {
        if (self.curr+1) > self.data.len() {
            self.curr = 0;
            return None
        }
        let idx = self.curr;
        self.curr+=1;
        Some(self.data[idx].clone())
    }
}


impl<T> fmt::Display for Tile<T> 
    where T:  Iterator<Item=u8> + ExactSizeIterator + fmt::Display + From<String>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for x in &self.data {
            write!(f,"{}\n", x);
        }
        Ok(())
    }
}
=======
}

impl<T: ExactSizeIterator + Iterator + fmt::Display> fmt::Display for Tile<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for x in &self.data {
            write!(f,"{}\n", x);
        }
        Ok(())
    }
}

// Tile <-> Array3<u16>
impl<T: ExactSizeIterator + Iterator + fmt::Display> From<&Tile<T>> for Array2<u16> {
    fn from(tile: &Tile<T>) -> Self {
        let mut converted = Array2::<u16>::zeros((tile.data.len(), tile.data.first().iter().len()));
        for x in tile.data.iter() {
            for value in x.iter() {
                // converted[(i,j)] = value;
            }
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
>>>>>>> f2d6b30b5dd36c6aab6feff01d6f2848ac118c25
