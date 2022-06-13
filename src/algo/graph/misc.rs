use std::collections::{HashSet, HashMap};
use std::collections::VecDeque;
use std::fmt;

use crate::ds::tile::Tile;
use crate::ds::sequence::Sequence;
use crate::ds::graph::Graph;
use crate::ds::graph::GraphProperties;

use super::dfs::Dfs;


pub fn overlap_graph(sequences: &Tile, k: usize) -> Graph::<Sequence, u8> {

    // Instantiate empty graph
    // let gp = ;
    let mut g = Graph::<Sequence, u8>::new(GraphProperties{directed: true});
    
    // Add all nodes to  graph
    let mut node_ids : Vec<u64> = vec![];
    for seq in sequences {
        node_ids.push(g.add_node(seq.clone()));
    }

    // Connect overlap graph
    for (i, seq) in sequences.into_iter().enumerate() {
        let last = seq.into_iter().rev().take(k).rev(); 
        for (j, seq2) in sequences.into_iter().enumerate() {
            if j != i {
                let first = seq2.into_iter().take(k);
                // Check if suffix of `seq` is equal to prefix of `seq2`
                if first.zip(last.clone()).filter(|&(a, b)| a != b).count() == 0 {
                    g.add_edge(&node_ids[i], &node_ids[j], None).unwrap();
                }
            }
        }
    }
    g
}

// Connected components based on Kosaraju's Algorithm O(V+E)
pub fn connected_components<N, E>(g: &mut Graph<N, E>) -> Vec<HashSet::<u64>> 
            where N: fmt::Display + Clone , E: fmt::Display + Clone
{

    // Define dfs with backtrack, where nodes are added to stack after all their children are processed
    fn dfs_backtrack<N, E>(g: &Graph<N, E>, nid: u64, stack: &mut Vec::<u64>, visited: &mut HashSet::<u64>) 
        where N: fmt::Display  + Clone, E: fmt::Display + Clone
    {

        if visited.contains(&nid) {
            return;
        }

        visited.insert(nid);
        
        for node in g.out_neighbors(nid) {
            dfs_backtrack(g, *node, stack, visited);
        }

        stack.push(nid);

    }

    // Run the backtrack DFS
    let mut visited = HashSet::<u64>::new();
    let mut stack = Vec::<u64>::new();

    for n in g.nodes() {
        if !visited.contains(n) {
            dfs_backtrack(g, *n, &mut stack, &mut visited);
        }
    }

    // Reverse graph
    g.reverse();

    // Stack based DFS, no backtracking needed
    let mut dfs = Dfs::new(g);
    let mut cc = HashSet::<u64>::new();
    let mut components : Vec<HashSet::<u64>> = vec![];
    let mut processed = HashSet::<u64>::new();

    while !stack.is_empty() {
        let n = stack.pop().unwrap();
        if !processed.contains(&n) {
            dfs.init(&n);
            while let Ok(id) = dfs.process_node() {
                processed.insert(id);
                cc.insert(id);
            }
            components.push(cc.clone());
            cc.clear();
        }
    }

    // Reverse graph back to its initial form
    g.reverse();

    components
}

// pub fn cycles<N: fmt::Display , E: fmt::Display + Clone >(g: &Graph<N, E>) -> Vec<Vec<u64>> {

//     let mut processing = HashSet::<u64>::new();    
//     let mut visited = Vec::<u64>::new();
//     let mut dfs_stack = Vec::<u64>::new();
//     // Cycle metadata
//     let mut cycle_trace = HashMap::<u64, u64>::new();
//     let mut cycle_origins = Vec::<u64>::new();
//     let mut cycles = Vec::<Vec<u64>>::new();

//     let mut cur_node_id : u64;

//     for node in g.nodes() {
//         dfs_stack.push(*node);

//         while !dfs_stack.is_empty() {
//             cur_node_id = dfs_stack.pop().unwrap();

//             if !visited.contains(&cur_node_id) { 
//                 if processing.contains(&cur_node_id) { 
//                     cycle_origins.push(cur_node_id); 
//                 }
//                 else {
//                     processing.insert(cur_node_id);
//                     for s in g.out_neighbors(cur_node_id) {
//                         cycle_trace.insert(*s, cur_node_id);
//                         dfs_stack.push(*s);
//                     }
//                 }
//             }

//         }

//         visited.extend(processing.iter());
//         processing.clear();

//         for n in cycle_origins.iter() {
//             let mut cycle = Vec::<u64>::new();
            
//             let mut predecessor = cycle_trace[n];
//             while predecessor != *n {
//                 cycle.push(predecessor);
//                 predecessor = cycle_trace[&predecessor];
//             }

//             cycle.push(*n);
//             cycle.reverse();
//             cycles.push(cycle);
//         } 
        
//         cycle_trace.clear();
//         cycle_origins.clear();
//     }

//     cycles
// }