use std::collections::{HashSet, HashMap};
use std::fmt;

use crate::ds::graph::Graph;

// Find cycles within a graph using Johnson's algorithm
// TODO: Move cycle detection as member of graph
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

// TODO: Move eulerian circuit detection as member of graph
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

