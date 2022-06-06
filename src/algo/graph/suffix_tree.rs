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
    graph: Graph<SuffixTreeNode, SuffixTreeEdge>,
    // Resolved
    previous_new_node : u64,
    // Postprocessing metadata
    alphabet: HashSet<u8>,
}

impl SuffixTreeBuilder {

    pub fn new(alphabet: &[u8]) -> SuffixTreeBuilder {

        let mut graph = Graph::<SuffixTreeNode, SuffixTreeEdge>::new(GraphProperties{directed: true});
        let root_id = graph.add_node(SuffixTreeNode::new());
        graph.set_root(root_id);

        SuffixTreeBuilder {
            // seq : sequence, 
            idx : 0,
            remaining : 0,
            active_node_id : root_id, 
            active_edge_seq_idx : 0,
            active_length : 0,
            end : 0,
            root_id : root_id,
            graph :  graph,
            previous_new_node : root_id,
            alphabet: HashSet::from_iter(alphabet.iter().cloned()),
        }
    }    

    pub fn build<'a, X, T>(&mut self, wordlist: &'a X) -> &mut Graph<SuffixTreeNode, SuffixTreeEdge>
                where &'a X: 'a + IntoIterator<Item=&'a T>,
                        &'a T: 'a + fmt::Display + Clone + Index<usize, Output=u8> + IntoIterator<Item=&'a u8>
    {    
        let mut separator : u8 = '!' as u8; 
        let mut temp : Vec<u8> = vec![];
        let mut wordmap : Vec<(usize, usize)> = vec![];
        let mut tile_len = 0_usize;


        for (idx, a) in wordlist.into_iter().enumerate() {
            println!("{}", idx);
            let lenold = temp.len();
            temp.extend(a);
            temp.push(separator);
            wordmap.extend(vec![(idx, temp.len() - 1); temp.len()-lenold + 1]);
            separator += 1;
            while self.alphabet.contains(&separator) {
                separator += 1;
            }
            tile_len += 1;
        }
    
        let g = self.process(&temp);

        fn generate_reachbility_map(graph: &mut Graph<SuffixTreeNode, SuffixTreeEdge>, node_id: u64, wordmap: &Vec<(usize, usize)>)
        {
            let out_neighbors : Vec<u64> = graph.out_neighbors(node_id).cloned().collect();
            for t in out_neighbors{
                generate_reachbility_map(graph, t, wordmap);
            }
    
            let mut reach = vec![0; wordmap.last().unwrap().0 + 1];
            let out_edges : Vec<u64> = graph.out_edges(node_id).cloned().collect();
    
            for eid in out_edges{
    
                let successor_node_id = graph.get_edge(&eid).end;
                let suffix_start = graph.get_edge(&eid).data.as_ref().unwrap().suffix_start;
                let suffix_stop = graph.get_edge(&eid).data.as_ref().unwrap().suffix_stop;
    
                if suffix_stop == -1 {
                    reach[wordmap[suffix_start].0] = 1;
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
    
        let root = self.graph.get_root().unwrap();
        generate_reachbility_map(&mut self.graph, root, &wordmap);
    
        // // Function to perform DFS traversal on the graph
        fn resolve_suffix_endings(graph: &mut Graph<SuffixTreeNode, SuffixTreeEdge>, node_id: u64, 
                                     wordmap: &Vec<(usize, usize)>)
        {
            let out_neighbors : Vec<u64> = graph.out_neighbors(node_id).cloned().collect();
            for t in out_neighbors{
                resolve_suffix_endings(graph, t, wordmap);
            }
    
            let out_edges : Vec<u64> = graph.out_edges(node_id).cloned().collect();
            for eid in out_edges {            
                let suffix_start = graph.get_edge(&eid).data.as_ref().unwrap().suffix_start;
                let suffix_stop = graph.get_edge(&eid).data.as_ref().unwrap().suffix_stop;
                if suffix_stop == -1 {
                    graph.get_edge_mut(&eid).data.as_mut().unwrap().suffix_stop = wordmap[suffix_start].1 as i64;
                }
            }
        }
    
        // visited.clear();
        resolve_suffix_endings(&mut self.graph, root, &wordmap);
    
        let mut lcs : Vec<(usize, i64)> = vec![];
        let mut cur : Vec<(usize, i64)> = vec![];
    
        let mut longest = 0;
        let mut cur_len = 0;
    
        fn dfs_recursive_substrings(graph: &Graph<SuffixTreeNode, SuffixTreeEdge>, node_id: u64, 
                                            cur_suffix: &mut Vec<(usize, i64)>, cur_length: usize, 
                                                lcs: &mut Vec<(usize, i64)>, longest: &mut usize, len: u8)
        {
            // mark the current node as discovered
    
            if cur_length > *longest {
                *longest = cur_length;
                *lcs = cur_suffix.clone();
            } 
    
            // do for every edge (v, u)
            for e in graph.out_edges(node_id){
                if graph.get_node(&graph.get_edge(e).end).data.reachable_suffixes.iter().sum::<u64>() == len as u64{
                    let start = graph.get_edge(e).data.as_ref().unwrap().suffix_start;
                    let stop = graph.get_edge(e).data.as_ref().unwrap().suffix_stop;
                    cur_suffix.push((start, stop));
                    dfs_recursive_substrings(graph, graph.get_edge(e).end, cur_suffix, cur_length + stop as usize - start + 1, lcs, longest, len);
                    cur_suffix.pop();
                }
            }
    
    
        }
    
        // visited.clear();
        dfs_recursive_substrings(&mut self.graph, root, &mut cur, cur_len, &mut lcs, &mut longest, tile_len as u8);
    
    
    
        let mut result = vec![];
    
        println!("{:?}", lcs);
    
        for part in lcs {
    
            for i in part.0..(part.1 as usize + 1) {
                result.push(temp[i] as u8);
            }
        }
    
        let x =  Sequence::from(result.as_slice());
    
        print!("{}", x);
        &mut self.graph

    }

    fn process(&mut self, seq: &Vec<u8>) -> () {
        for s in seq {
            self.step(seq);
            self.idx += 1;
        }
    }

    fn split_suffix_edge(graph: &mut Graph<SuffixTreeNode, SuffixTreeEdge>, edge_id: u64, split_index: usize, value_index: usize, seq: &Vec<u8>) -> u64 {
        
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

    pub fn step(&mut self, seq: &Vec<u8>) {

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

            let cur_node  = self.graph.get_node(&self.active_node_id);
            cur_value = seq[self.active_edge_seq_idx];

            // Check if there is a suffix corresponding to `cur_value` adjacent to current node, 
            // Add a new None/Edge to the graph, if there is no suffix corresponding
            // Activate the edge that corresponds to the suffix otherwise
            if !cur_node.data.suffix_edge_ids.contains_key(&cur_value) {

                // Create a new node and edge
                let nid = self.graph.add_node(SuffixTreeNode::new());
                let sfx = SuffixTreeEdge{suffix_start: self.idx, suffix_stop: -1};
                let eid = self.graph.add_edge(&self.active_node_id, &nid, Some(sfx)).unwrap();

                // Keep track of newly added edge (suffix) inside current node
                let cur_node_mut  = self.graph.get_node_mut(&self.active_node_id);
                cur_node_mut.data.suffix_edge_ids.insert(cur_value, eid);          
                
                if self.previous_new_node != self.root_id {
                    self.graph.get_node_mut(&self.previous_new_node).data.link = self.active_node_id;
                    self.previous_new_node = self.active_node_id;
                }

                self.remaining -= 1;
            }
            else {

                // Characters might match the suffix on the path from root downwards
                // This might introduce an update to the active node, if the current edge gets exhausted
                let cur_node  = self.graph.get_node(&self.active_node_id);
                let active_edge_identifier = cur_node.data.suffix_edge_ids[&seq[self.active_edge_seq_idx]];
                let active_edge = self.graph.get_edge(&active_edge_identifier).clone();
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
                        self.graph.get_node_mut(&self.previous_new_node).data.link = self.active_node_id;
                        self.previous_new_node = self.active_node_id;
                    }
                    break;
                }

                // Resolve all partial suffixes have been accumulated during transition along edge
                // Break the matching on the active edge and introduce new leaf node for new suffix
                let split_node_id = SuffixTreeBuilder::split_suffix_edge(
                    &mut self.graph, 
                    active_edge_identifier, 
                    self.active_length, 
                    self.idx, 
                    seq
                );
                self.graph.get_node_mut(&split_node_id).data.link = self.root_id;

                // Update suffix links on previous node, if subsequent internal node created
                if self.previous_new_node != self.root_id {
                    self.graph.get_node_mut(&self.previous_new_node).data.link = split_node_id;
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
                self.active_node_id = self.graph.get_node(&self.active_node_id).data.link;
            }
        }  
    }
}

/// SuffixTreeBuilder Node
#[derive(Clone)]
pub struct SuffixTreeNode {
    link: u64,
    suffix_edge_ids: HashMap<u8, u64>, 
    reachable_suffixes: Vec<u64>,
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