use std::fmt; // Import `fmt`
use ndarray::prelude::*;
use std::ops::{Index, IndexMut};

use crate::sequence::Sequence;

#[derive(Debug)]
pub struct Tile {
    data: Vec<Sequence>,
}

impl Tile
{
    pub fn new() -> Tile {
        Tile{ data: Vec::new() }
    }

    pub fn push(&mut self, value: Sequence) {
        self.data.push(value);
    }

    pub fn pop(&mut self) -> Option<Sequence> {
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

/*** Type-Conversion Traits ***/ 
// Array2<u16> <-> DNA
impl From<Array2<u16>> for Tile 
{
    fn from(arr: Array2<u16>) -> Self {
        let mut temp = Tile::new();
        for (i, ax) in arr.axis_iter(Axis(0)).enumerate() {
            let mut s : String = "".to_string();
            for (j, value) in ax.indexed_iter() {
                s.push((*value as u8) as char)
            }
            temp.push(Sequence::from(s));
        }
        temp
    }
}
// Array3<u16> <-> DNA
impl From<Array3<u16>> for Tile 
{
    fn from(arr: Array3<u16>) -> Self {
        let mut temp = Tile::new();
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
            temp.push(Sequence::from(s));
        }
        temp
    }
}

/*** Utility Traits ***/ 
// Iterator Trait
impl IntoIterator for Tile 
{
    type Item = Sequence;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl<'a> IntoIterator for &'a Tile {

    type Item = <std::slice::Iter<'a, Sequence> as Iterator>::Item;
    type IntoIter = std::slice::Iter<'a, Sequence>;

    fn into_iter(self) -> Self::IntoIter {
        (&self.data).into_iter() // can also be .iter()
    }
}

/*** Debug Traits */
// Display
impl fmt::Display for Tile 
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for x in &self.data {
            write!(f,"{}\n", x);
        }
        Ok(())
    }
}

impl Index<usize> for Tile {
    type Output = Sequence;
    fn index<'a>(&'a self, i: usize) -> &'a Sequence {
        &self.data[i]
    }
}

impl IndexMut<usize> for Tile {
    fn index_mut<'a>(&'a mut self, i: usize) -> &'a mut Sequence {
        &mut self.data[i]
    }
}