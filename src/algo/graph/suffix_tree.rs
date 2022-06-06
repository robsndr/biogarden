use std::fmt;
use std::collections::{HashMap, HashSet};
use std::ops::{Index, IndexMut};

use crate::ds::tile::Tile;
use crate::ds::sequence::Sequence;
use crate::ds::graph::Graph;
use crate::ds::graph::GraphProperties;

// SuffixTreeBuilder Algorithm
pub struct SuffixTreeBuilder {
    // Sequence to be processed
    // seq: T, 
    idx: usize,
    // State indicators
    remaining : usize,
    active_node_id : u64, 
    active_edge_seq_idx : usize,
    active_length : usize,
    end: u8,
    // Graph 
    root_id : u64,
    // graph: Graph<SuffixTreeNode, SuffixTreeEdge>,
    // Resolved
    previous_new_node : u64,
    // Postprocessing metadata
    alphabet: HashSet<u8>,
}

impl SuffixTreeBuilder {

    pub fn new(alphabet: HashSet<u8>) -> SuffixTreeBuilder {

        SuffixTreeBuilder {
            idx : 0,
            remaining : 0,
            active_node_id : 0, 
            active_edge_seq_idx : 0,
            active_length : 0,
            end : 0,
            root_id : 0,
            // graph :  graph,
            previous_new_node : 0,
            alphabet: alphabet,
        }
    }    

    pub fn build<'a, T>(&mut self, seq: &'a T) -> Graph<SuffixTreeNode, SuffixTreeEdge>
            where &'a T: 'a + fmt::Display + Clone + Index<usize, Output=u8> + IntoIterator<Item=&'a u8>
    {
        // Initialize Graph
        let mut graph = Graph::<SuffixTreeNode, SuffixTreeEdge>::new(GraphProperties{directed: true});
        let root_id = graph.add_node(SuffixTreeNode::new());
        graph.set_root(root_id);
        
        // Initialize variables
        self.active_node_id = root_id;
        self.previous_new_node = root_id;
        self.root_id = root_id;
        
        // Build SuffixTree
        seq.into_iter().for_each(|_| { self.step(&mut graph, seq); });

        // Endings are encoded by -1 in ukonnen's algorithm
        // Wordmap is required to decode suffix ending enumerations 
        let mut wordmap = Vec::<(usize, usize)>::new();
        let mut word_len = 0_usize;
        let mut word_idx = 0_usize;
        for (id, x) in seq.into_iter().enumerate() {
            word_len += 1;
            if !self.alphabet.contains(&x) {
                wordmap.extend(vec![(word_idx, id); word_len]);
                word_len = 0;
                word_idx += 1;
            }
        }
        
        // Decode each -1 at leaf into a valid index pointing into 'seq'
        SuffixTreeBuilder::postprocess(&mut graph, root_id, &wordmap);
        graph
    }

    fn split_suffix_edge<'a, T>(graph: &mut Graph<SuffixTreeNode, SuffixTreeEdge>, edge_id: u64, split_index: usize, value_index: usize, seq: &'a T) -> u64 
            where &'a T: 'a + fmt::Display + Clone + Index<usize, Output=u8> + IntoIterator<Item=&'a u8>    
    {
        
        // Get indizes of nodes in graph that are connected to edge
        let edge_start = graph.get_edge(&edge_id).start;
        let edge_end = graph.get_edge(&edge_id).end;

        // Get indizes of suffix that is encoded by the given edge
        let suffix_start = graph.get_edge(&edge_id).data.as_ref().unwrap().suffix_start;
        let suffix_stop = graph.get_edge(&edge_id).data.as_ref().unwrap().suffix_stop;

        // Split edge
        let split_node_id = graph.add_node(SuffixTreeNode::new());
        let split_edge_pre = graph.add_edge(
            &edge_start, 
            &split_node_id, 
            Some(SuffixTreeEdge{suffix_start: suffix_start, suffix_stop: (suffix_start + split_index - 1) as i64})
        ).unwrap();
        let split_edge_succ = graph.add_edge(
            &split_node_id, 
            &edge_end, 
            Some(SuffixTreeEdge{suffix_start: suffix_start + split_index, suffix_stop: suffix_stop as i64})
        ).unwrap();

        // Add value into node created at split
        let leaf_node = graph.add_node(SuffixTreeNode::new());
        let value_edge = graph.add_edge(
            &split_node_id, 
            &leaf_node, 
            Some(SuffixTreeEdge{suffix_start: value_index, suffix_stop: -1})
        ).unwrap();

        // Remove old edge that was split
        graph.remove_edge(&edge_id);
        graph.get_node_mut(&edge_start).data.suffix_edge_ids.remove(&seq[suffix_start]);

        // Update mapping for outgoing suffix edges for the involved nodes
        graph.get_node_mut(&edge_start).data.suffix_edge_ids.insert(seq[suffix_start], split_edge_pre);
        graph.get_node_mut(&split_node_id).data.suffix_edge_ids.insert(seq[suffix_start + split_index], split_edge_succ);
        graph.get_node_mut(&split_node_id).data.suffix_edge_ids.insert(seq[value_index], value_edge);

        split_node_id
    }

    pub fn step<'a, T>(&mut self, graph: &mut Graph<SuffixTreeNode, SuffixTreeEdge>, seq: &'a T) 
        where &'a T: 'a + fmt::Display + Clone + Index<usize, Output=u8> + IntoIterator<Item=&'a u8>
    {

        let mut cur_value = seq[self.idx];
        // Keeps track of the number of `suffixes` to be resolved after a new character is added 
        self.remaining += 1;
        // Set the current `end` character to new value
        self.end = cur_value;
        // Set previous new node for suffix links back to root 
        self.previous_new_node = self.root_id;

        // For given position `self.idx` in input sequence `self.seq`
        // Iterate as long as suffix after adding given letter is resolved
        while self.remaining > 0 { 
            
            if self.active_length == 0 {
                self.active_edge_seq_idx = self.idx;
            }

            let cur_node  = graph.get_node(&self.active_node_id);
            cur_value = seq[self.active_edge_seq_idx];

            // Check if there is a suffix corresponding to `cur_value` adjacent to current node, 
            // Add a new None/Edge to the graph, if there is no suffix corresponding
            // Activate the edge that corresponds to the suffix otherwise
            if !cur_node.data.suffix_edge_ids.contains_key(&cur_value) {

                // Create a new node and edge
                let nid = graph.add_node(SuffixTreeNode::new());
                let sfx = SuffixTreeEdge{suffix_start: self.idx, suffix_stop: -1};
                let eid = graph.add_edge(&self.active_node_id, &nid, Some(sfx)).unwrap();

                // Keep track of newly added edge (suffix) inside current node
                let cur_node_mut  = graph.get_node_mut(&self.active_node_id);
                cur_node_mut.data.suffix_edge_ids.insert(cur_value, eid);          
                
                if self.previous_new_node != self.root_id {
                    graph.get_node_mut(&self.previous_new_node).data.link = self.active_node_id;
                    self.previous_new_node = self.active_node_id;
                }

                self.remaining -= 1;
            }
            else {

                // Characters might match the suffix on the path from root downwards
                // This might introduce an update to the active node, if the current edge gets exhausted
                let cur_node  = graph.get_node(&self.active_node_id);
                let active_edge_identifier = cur_node.data.suffix_edge_ids[&seq[self.active_edge_seq_idx]];
                let active_edge = graph.get_edge(&active_edge_identifier).clone();
                let lookup_idx = self.active_length + active_edge.data.as_ref().unwrap().suffix_start;

                // Check if an internal node has been reach, and update state accordingly
                let suffix_stop = active_edge.data.as_ref().unwrap().suffix_stop;
                let suffix_start = active_edge.data.as_ref().unwrap().suffix_start;

                if suffix_stop != -1 && lookup_idx > suffix_stop as usize {
                    self.active_node_id = active_edge.end;
                    self.active_length = lookup_idx - suffix_stop as usize - 1;
                    self.active_edge_seq_idx += (suffix_stop as usize - suffix_start) + 1; // +1?
                    continue;
                }

                // Check if next character on active edge matches the current value
                // If matches, we proceed along the edge and go to next character
                if seq[lookup_idx] == seq[self.idx] {
                    self.active_length += 1;
                    if self.previous_new_node != self.root_id {
                        graph.get_node_mut(&self.previous_new_node).data.link = self.active_node_id;
                        self.previous_new_node = self.active_node_id;
                    }
                    break;
                }

                // Resolve all partial suffixes have been accumulated during transition along edge
                // Break the matching on the active edge and introduce new leaf node for new suffix
                let split_node_id = SuffixTreeBuilder::split_suffix_edge(
                    graph, 
                    active_edge_identifier, 
                    self.active_length, 
                    self.idx, 
                    seq
                );
                graph.get_node_mut(&split_node_id).data.link = self.root_id;

                // Update suffix links on previous node, if subsequent internal node created
                if self.previous_new_node != self.root_id {
                    graph.get_node_mut(&self.previous_new_node).data.link = split_node_id;
                }
                
                self.previous_new_node = split_node_id;
                self.remaining -= 1;
                
            }

            // Jump to next edge if are on an edge from root and length > 0, 
            // Otherwise we are at some internal node, jump to next internal node as per suffix link resolution rules
            if self.active_length > 0 && self.active_node_id == self.root_id {
                // Current character was processed, reduce length
                self.active_length -= 1;
                self.active_edge_seq_idx = self.idx - self.remaining + 1;
            }
            else if self.active_node_id != self.root_id {      
                self.active_node_id = graph.get_node(&self.active_node_id).data.link;
            }
        }  
        self.idx += 1;
    }

    fn postprocess(graph: &mut Graph<SuffixTreeNode, SuffixTreeEdge>, node_id: u64, wordmap: &Vec<(usize, usize)>)
    {
        let out_neighbors : Vec<u64> = graph.out_neighbors(node_id).cloned().collect();
        for t in out_neighbors{
            SuffixTreeBuilder::postprocess(graph, t, wordmap);
        }

        let mut reach = vec![0; wordmap.last().unwrap().0 + 1];
        let out_edges : Vec<u64> = graph.out_edges(node_id).cloned().collect();

        for eid in out_edges {

            let successor_node_id = graph.get_edge(&eid).end;
            let suffix_start = graph.get_edge(&eid).data.as_ref().unwrap().suffix_start;
            let suffix_stop = graph.get_edge(&eid).data.as_ref().unwrap().suffix_stop;

            if suffix_stop == -1 {
                reach[wordmap[suffix_start].0] = 1;
                graph.get_edge_mut(&eid).data.as_mut().unwrap().suffix_stop = wordmap[suffix_start].1 as i64;
            }
            else {

                for (i, elem) in graph.get_node(&successor_node_id).data.reachable_suffixes.iter().enumerate() {
                    if *elem == 1 {
                        reach[i] = 1;
                    }
                };
            }
        }

        graph.get_node_mut(&node_id).data.reachable_suffixes = reach;
    } 
}

/// SuffixTreeBuilder Node
#[derive(Clone)]
pub struct SuffixTreeNode {
    link: u64,
    suffix_edge_ids: HashMap<u8, u64>, 
    pub reachable_suffixes: Vec<u64>,
}

impl SuffixTreeNode {
    pub fn new() -> SuffixTreeNode {
        SuffixTreeNode {
            link: 0,
            suffix_edge_ids: HashMap::new(),
            reachable_suffixes: Vec::<u64>::new()
        }
    }
}

impl fmt::Display for SuffixTreeNode  {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "UN: Link {:?}", self.link)
    }
}

// SuffixTreeBuilder Edge
#[derive(Clone)]
pub struct SuffixTreeEdge {
    pub suffix_start: usize,
    pub suffix_stop: i64
}

impl SuffixTreeEdge {
    pub fn new() -> SuffixTreeEdge {
        SuffixTreeEdge {
            suffix_start: 0,
            suffix_stop: 0,
        }
    }
}

impl fmt::Display for SuffixTreeEdge  {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{},{}]", self.suffix_start, self.suffix_stop)
    }
}