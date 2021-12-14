use crate::dna;
use std::fmt; // Import `fmt`
use ndarray::prelude::*;


#[derive(Debug)]
pub struct Tile<T> 
    where T:  Iterator<Item=u8> + ExactSizeIterator + fmt::Display + From<String> 
{
    data: Vec<T>,
    curr: usize
}

impl<'a, T> Tile<T>
    where T: 'a + Iterator<Item=u8> + ExactSizeIterator + fmt::Display + From<String>  
{
    pub fn new() -> Tile<T> {
        Tile::<T>{ data: Vec::new(), curr: 0 }
    }

    pub fn push(&mut self, value: T) {
        self.data.push(value);
    }

    pub fn pop(&mut self) -> Option<T> {
        self.data.pop()
    }

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
impl<'a, T> From<Array2<u16>> for Tile<T> 
    where T: 'a + Iterator<Item=u8> + ExactSizeIterator + fmt::Display + From<String> 
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
impl<'a, T> From<Array3<u16>> for Tile<T> 
    where T: 'a + Iterator<Item=u8> + ExactSizeIterator + fmt::Display + From<String> 
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
impl<'a, T> IntoIterator for &'a Tile<T> 
    where T: 'a + Iterator<Item=u8> + ExactSizeIterator + fmt::Display + From<String> 
{
    type Item = <std::slice::Iter<'a, T> as Iterator>::Item;
    type IntoIter = std::slice::Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.as_slice().into_iter()
    }
}


impl<'a, T> fmt::Display for Tile<T> 
    where T: 'a + Iterator<Item=u8> + ExactSizeIterator + fmt::Display + From<String>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for x in &self.data {
            write!(f,"{}\n", x);
        }
        Ok(())
    }
}
