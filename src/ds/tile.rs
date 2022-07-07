use std::fmt; // Import `fmt`
use ndarray::prelude::*;
use std::ops::{Index, IndexMut};

use super::sequence::Sequence;

#[derive(Debug, Clone)]
pub struct Tile {
    pub data: Vec<Sequence>,
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

    pub fn remove(&mut self, index: usize) ->  Sequence {
        self.data.remove(index)
    }

    pub fn size(& self) -> (usize, usize) {
        (self.data.len(), self.data[0].len())
    }

    pub fn len(& self) -> usize {
        self.data.len()
    }

    pub fn is_empty(& self) -> bool {
        self.data.is_empty()
    }

    pub fn extend(&mut self, b: Tile) {
        self.data.extend(b.data);
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

impl Default for Tile {
    fn default() -> Self {
        Self::new()
    }
}

/*** Type-Conversion Traits ***/ 
// Array2<u16> <-> DNA
impl From<Array2<u16>> for Tile 
{
    fn from(arr: Array2<u16>) -> Self {
        let mut temp = Tile::new();
        for ax in arr.axis_iter(Axis(0)) {
            let mut s : String = "".to_string();
            for &value in ax.iter() {
                s.push((value as u8) as char)
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
        for ax1 in arr.axis_iter(Axis(0)) {
            let mut s : String = "".to_string();
            for ax2 in ax1.axis_iter(Axis(0)) {
                for (k, &val) in ax2.indexed_iter() {
                    if val == 1 {
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

impl From<&[Sequence]> for Tile {
    fn from(s: &[Sequence]) -> Self {
        Tile { data: s.to_vec() }
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
        (&self.data).iter() // can also be .iter()
    }
}

/*** Debug Traits */
// Display
impl fmt::Display for Tile 
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for x in &self.data {
            writeln!(f,"{}", x)?;
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
    fn index_mut(&'_ mut self, i: usize) -> &'_ mut Sequence {
        &mut self.data[i]
    }
}

impl Eq for Tile {}

impl PartialEq for Tile {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}
