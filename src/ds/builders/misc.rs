use std::collections::{HashSet, HashMap};
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


// Find cycles within a graph using Johnson's algorithm
pub fn cycles<N, E>(graph: &mut Graph<N, E>) -> Vec<Vec::<u64>> 
            where N: fmt::Display + Clone , E: fmt::Display + Clone
{

    fn unblock(cur_node: u64, blocked_set: &mut HashSet<u64>, blocked_map: &mut HashMap<u64, HashSet<u64>>) {
        let mut stack : Vec<u64> = vec![cur_node];
        while !stack.is_empty() {
            let node = stack.pop().unwrap();
            if blocked_set.contains(&node) {
                stack.extend(blocked_map.get(&node).unwrap_or(&HashSet::<u64>::new()).iter().cloned());
                blocked_set.remove(&node);
                blocked_map.remove(&node);
            }
        }
    }

    let mut cycles : Vec<Vec::<u64>> = vec![];

    let mut start_node : u64;
    let mut path : Vec<u64>;
    let mut blocked_set = HashSet::<u64>::new();
    let mut closed_set = HashSet::<u64>::new();
    let mut blocked_map = HashMap::<u64, HashSet::<u64>>::new();

    let mut stack : Vec<(u64, Vec<u64>)> = vec![];

    while graph.node_count() > 1 {

        start_node = *graph.nodes().next().unwrap();
        path = vec![start_node];
        
        blocked_set.clear();
        closed_set.clear();
        blocked_map.clear();
        stack.clear();

        blocked_set.insert(start_node);
        stack.push((start_node, graph.out_neighbors(start_node).cloned().collect()));

        while !stack.is_empty() {
            let (cur_node, neighbours) = stack.last_mut().unwrap();
            if !neighbours.is_empty() {
                let next_node = neighbours.pop().unwrap();
                if next_node == start_node {
                    cycles.push(path.clone());
                    closed_set.extend(path.clone());
                }
                else if !blocked_set.contains(&next_node) {
                    path.push(next_node);
                    stack.push((next_node, graph.out_neighbors(next_node).cloned().collect()));
                    closed_set.remove(&next_node);
                    blocked_set.insert(next_node);
                    continue
                }
            }
                
            if neighbours.is_empty() {
                if closed_set.contains(cur_node) {
                    unblock(*cur_node, &mut  blocked_set, &mut blocked_map);
                }
                else {
                    for nbr in graph.out_neighbors(*cur_node) {
                        if !blocked_map.contains_key(nbr) {
                            blocked_map.insert(*nbr, HashSet::<u64>::new());
                        }
                        if !blocked_map[nbr].contains(cur_node) {
                            blocked_map.get_mut(nbr).unwrap().insert(*cur_node);
                        }
                    }
                }
                stack.pop();
                path.pop();
            } 
        }

        graph.remove_node(start_node);
    }
    cycles
}


pub fn eulerian_circuit<N, E>(graph: &mut Graph<N, E>, source: u64) -> Vec<Vec<u64> >
            where N: fmt::Display + Clone , E: fmt::Display + Clone
{
    let mut stack : Vec<(u64, Graph<N, E>, Vec<u64>)> = vec![(source, graph.clone(), vec![])];
    let mut paths : Vec<Vec<u64>> = vec![vec![]];
    
    let edge_count = graph.edge_count();

    while !stack.is_empty() {
        let (cur_node, graph, mut cur_path) = stack.pop().unwrap();
        cur_path.push(cur_node);
        
        for neighbour in graph.out_neighbors(cur_node) {
            let mut g = graph.clone();
            let edge_id = g.has_edge(&cur_node, neighbour).unwrap();
            g.remove_edge(&edge_id);
            stack.push((*neighbour, g, cur_path.clone()))
        }
            

        if graph.node_degree(&cur_node).unwrap() == 0 {
            // Remove last entry from path -> redundant
            // Tail will be equivalent to the head of the circle
            cur_path.pop();
            if cur_path.len() == edge_count {
                paths.push(cur_path.clone());
            }
        }
    }
        
    // # Remove invalid paths
    paths
}

