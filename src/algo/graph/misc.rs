use std::collections::HashSet;

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

pub fn connected_components(g: &Graph<u64, u8>) -> (u32, Vec<Vec<u64>>) {

    let mut dfs = Dfs::new(g);

    let mut processed = HashSet::<u64>::new();
    let mut ctr : u32 = 0;
    let mut cc : Vec<u64> = vec![];
    let mut components : Vec<Vec<u64>> = vec![];
    
    for n in g.nodes() {
        if !processed.contains(n) {
            dfs.init(n);
            cc.clear();
            while let Ok(id) = dfs.process_node() {
                processed.insert(id);
                cc.push(id);
            }
            components.push(cc.clone());
            ctr += 1;
        }
    }
    (ctr, components)
}

