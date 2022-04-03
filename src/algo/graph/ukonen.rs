use std::fmt;
use std::collections::HashMap;
use std::ops::{Index, IndexMut};

use crate::ds::tile::Tile;
use crate::ds::sequence::Sequence;
use crate::ds::graph::Graph;
use crate::ds::graph::GraphProperties;


/// Ukonen Node
pub struct UkonenNode {
    link: u64,
    suffix_edge_ids: HashMap<u8, u64>
}

impl UkonenNode {
    pub fn new() -> UkonenNode {
        UkonenNode {
            link: 0,
            suffix_edge_ids: HashMap::new(),
        }
    }
}

impl fmt::Display for UkonenNode  {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "UN: Link {:?}", self.link)
    }
}

// Ukonen Edge
#[derive(Clone)]
pub struct UkonenEdge {
    pub suffix_start: usize,
    pub suffix_stop: i64
}

impl UkonenEdge {
    pub fn new() -> UkonenEdge {
        UkonenEdge {
            suffix_start: 0,
            suffix_stop: 0,
        }
    }
}

impl fmt::Display for UkonenEdge  {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{},{}]", self.suffix_start, self.suffix_stop)
    }
}

// Ukonen Algorithm
pub struct Ukonen<T:fmt::Display + Index<usize, Output=u8> + IntoIterator> {
    // Sequence to be processed
    seq: T, 
    idx: usize,
    // State indicators
    remaining : usize,
    active_node_id : u64, 
    active_edge_seq_idx : usize,
    active_length : usize,
    end: u8,
    // Graph 
    root_id : u64,
    graph: Graph<UkonenNode, UkonenEdge>,
    // Resolved
    previous_new_node : u64,

}

impl<T: fmt::Display + Clone + Index<usize, Output=u8> + IntoIterator> Ukonen<T> {

    pub fn new(sequence: T) -> Ukonen<T> {

        let mut graph = Graph::<UkonenNode, UkonenEdge>::new(GraphProperties{directed: true});
        let root_id = graph.add_node(UkonenNode::new());
        graph.set_root(root_id);

        Ukonen {
            seq : sequence, 
            idx : 0,
            remaining : 0,
            active_node_id : root_id, 
            active_edge_seq_idx : 0,
            active_length : 0,
            end : 0,
            root_id : root_id,
            graph :  graph,
            previous_new_node : root_id,
        }
    }    

    pub fn process(&mut self) -> &mut Graph<UkonenNode, UkonenEdge>{
        for s in (self.seq.clone()).into_iter() {
            self.step();
            self.idx += 1;
        }
        &mut self.graph
    }

    fn split_suffix_edge(graph: &mut Graph<UkonenNode, UkonenEdge>, edge_id: u64, split_index: usize, value_index: usize, seq: &T) -> u64 {
        
        // Get indizes of nodes in graph that are connected to edge
        let edge_start = graph.get_edge(&edge_id).start;
        let edge_end = graph.get_edge(&edge_id).end;

        // Get indizes of suffix that is encoded by the given edge
        let suffix_start = graph.get_edge(&edge_id).data.as_ref().unwrap().suffix_start;
        let suffix_stop = graph.get_edge(&edge_id).data.as_ref().unwrap().suffix_stop;

        // Split edge
        let split_node_id = graph.add_node(UkonenNode::new());
        let split_edge_pre = graph.add_edge(
            &edge_start, 
            &split_node_id, 
            Some(UkonenEdge{suffix_start: suffix_start, suffix_stop: (suffix_start + split_index - 1) as i64})
        ).unwrap();
        let split_edge_succ = graph.add_edge(
            &split_node_id, 
            &edge_end, 
            Some(UkonenEdge{suffix_start: suffix_start + split_index, suffix_stop: suffix_stop as i64})
        ).unwrap();

        // Add value into node created at split
        let leaf_node = graph.add_node(UkonenNode::new());
        let value_edge = graph.add_edge(
            &split_node_id, 
            &leaf_node, 
            Some(UkonenEdge{suffix_start: value_index, suffix_stop: -1})
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

    pub fn step(&mut self) {

        let mut cur_value = self.seq[self.idx];
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
            cur_value = self.seq[self.active_edge_seq_idx];

            // Check if there is a suffix corresponding to `cur_value` adjacent to current node, 
            // Add a new None/Edge to the graph, if there is no suffix corresponding
            // Activate the edge that corresponds to the suffix otherwise
            if !cur_node.data.suffix_edge_ids.contains_key(&cur_value) {

                // Create a new node and edge
                let nid = self.graph.add_node(UkonenNode::new());
                let sfx = UkonenEdge{suffix_start: self.idx, suffix_stop: -1};
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
                let active_edge_identifier = cur_node.data.suffix_edge_ids[&self.seq[self.active_edge_seq_idx]];
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
                if self.seq[lookup_idx] == self.seq[self.idx] {
                    self.active_length += 1;
                    if self.previous_new_node != self.root_id {
                        self.graph.get_node_mut(&self.previous_new_node).data.link = self.active_node_id;
                        self.previous_new_node = self.active_node_id;
                    }
                    break;
                }

                // Resolve all partial suffixes have been accumulated during transition along edge
                // Break the matching on the active edge and introduce new leaf node for new suffix
                let split_node_id = Ukonen::split_suffix_edge(
                    &mut self.graph, 
                    active_edge_identifier, 
                    self.active_length, 
                    self.idx, 
                    &self.seq
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