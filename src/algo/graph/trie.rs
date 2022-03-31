use std::fmt;

use std::collections::HashMap;
use std::ops::{Index, IndexMut};

use crate::ds::graph::Graph;
use crate::ds::graph::GraphErr;
use crate::ds::graph::GraphProperties;

/// Ukonen Node
pub struct TrieNode {
    pub substring: String, 
    pub children: Vec<i64>,
    pub ending: bool
}

impl TrieNode {
    pub fn new(alphabet_len: usize) -> TrieNode {
        TrieNode {
            substring: String::from("-"), 
            children: vec![-1; alphabet_len],
            ending: false
        }
    }
}

impl fmt::Display for TrieNode  {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // let mut temp : String = String::from("");
        // temp += std::str::from_utf8(&self.substring).unwrap();
        write!(f, "S: {} | E: {}", self.substring, self.ending)
    }
}

// TODO: make generic
pub struct Trie {
    graph: Graph<TrieNode, u8>,
    // map letters to positional representation
    alphabet: HashMap<u8, usize>
}

impl Trie {
    
    pub fn new(alphabet: &[u8]) -> Trie {

        let mut temp = HashMap::<u8, usize>::new();
        for (i, c) in alphabet.iter().enumerate() {
            temp.insert(*c, i);
        }

        let mut graph = Graph::<TrieNode, u8>::new(GraphProperties{directed: true});
        let root_id = graph.add_node(TrieNode::new(temp.len()));
        graph.set_root(root_id);

        Trie {
            graph: graph,
            alphabet: temp
        }
    }

    pub fn insert_word<'a, T>(&mut self, word: &'a T) -> Result<&Graph::<TrieNode, u8>, GraphErr> 
        where &'a T: 'a + Clone + IntoIterator<Item=&'a u8>
    {
        
        let mut cur_node_id = match self.graph.get_root() {
            Some(root) => root, 
            None => return Err(GraphErr::NoSuchVertex)
        };

        let mut substring = String::from("");

        for c in word.into_iter() { 
            
            // TODO: refine error handline
            let idx = match self.alphabet.get(c) {
                Some(x) => x.to_owned(),
                None => return Err(GraphErr::NoSuchVertex) 
            };
            
            if self.graph.get_node(&cur_node_id).data.children[idx] == -1 {
                let nid = self.graph.add_node(TrieNode::new(self.alphabet.len()));
                self.graph.add_edge(&cur_node_id, &nid, None).ok();                
                self.graph.get_node_mut(&cur_node_id).data.children[idx] = nid as i64;
            }
            
            let cur_node = self.graph.get_node_mut(&cur_node_id);
            let next_node_id = cur_node.data.children[idx] as u64;

            let next_node = self.graph.get_node_mut(&next_node_id);
            substring.push(*c as char);

            cur_node_id = next_node_id;
        }
     
        self.graph.get_node_mut(&cur_node_id).data.ending = true;
        self.graph.get_node_mut(&cur_node_id).data.substring = substring;

        Ok(&self.graph)
    }

    pub fn build<'a, X, T>(&mut self, wordlist: &'a X) -> Result<&Graph::<TrieNode, u8>, GraphErr>
        where &'a X: 'a + IntoIterator<Item=&'a T>,
                &'a T: 'a + Clone + IntoIterator<Item=&'a u8>
    {
        for word in wordlist {
            if self.insert_word(&word).is_err() {
                return Err(GraphErr::CannotAddEdge)
            }
        }
        
        Ok(&self.graph)
    }
}
