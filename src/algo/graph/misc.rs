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

fn unblock(nid: u64, blocked_set: &mut HashSet::<u64>, blocked_map: &mut HashMap::<u64, HashSet<u64>>) {
    blocked_set.remove(&nid);
    
    if !blocked_map.contains_key(&nid) {
        return;
    }

    let blocked_mapings = blocked_map.get(&nid).unwrap().clone();
    for n in blocked_mapings.iter() {
        blocked_map.get_mut(&nid).unwrap().remove(n);
        if blocked_set.contains(n) {
            unblock(*n, blocked_set, blocked_map);
        }
    }
}

// Find cycles within a graph using Johnson's algorithm
pub fn cycles<N, E>(g: &mut Graph<N, E>, test: u64) -> Vec<Vec::<u64>> 
            where N: fmt::Display + Clone , E: fmt::Display + Clone
{



    fn dfs_backtrack<N, E>(g: &Graph<N, E>, nid: u64, stack: &mut Vec::<u64>, blocked_set: &mut HashSet::<u64>, 
                            blocked_map: &mut HashMap::<u64, HashSet<u64>>, cycles: &mut Vec<Vec::<u64>>) -> bool
        where N: fmt::Display  + Clone, E: fmt::Display + Clone
    {
        println!("{}\n", nid+1);
        let mut f : bool = false;
        stack.push(nid);
        blocked_set.insert(nid);
        println!("ASD");
        for w in g.out_neighbors(nid) {
            if w == stack.first().unwrap() {
                cycles.push(stack.clone());
                f = true;
            } 
            else if !blocked_set.contains(w) {
                if dfs_backtrack(g, *w, stack, blocked_set, blocked_map, cycles) {
                    f = true;
                }
            }
        }

        if f {
            unblock(nid, blocked_set, blocked_map);
        } 
        else {
            for w in g.out_neighbors(nid) {
                //if (v∉B(w)) put v on B(w);
                if !blocked_map.contains_key(w) {
                    blocked_map.insert(*w, HashSet::new());
                }
                if !blocked_map.get(w).unwrap().contains(w) {
                    blocked_map.get_mut(w).unwrap().insert(nid);
                }
            }
        }
        // v = stack.pop();
        println!("{:?}", blocked_set);
        let v = stack.pop().unwrap();
        println!("POP: {}", v+1);

        return f;
    }


    let mut cycles : Vec<Vec::<u64>> = vec![];

    // while g.node_count() > 1 {

    let nid = g.nodes().next().expect("Can't get start node");
    let mut stack = Vec::<u64>::new();
    let mut blocked_set = HashSet::<u64>::new();
    let mut blocked_map = HashMap::<u64, HashSet::<u64>>::new();
            
    dfs_backtrack(g,test, &mut stack, &mut blocked_set, &mut blocked_map, &mut cycles);
    // }

    cycles
}
