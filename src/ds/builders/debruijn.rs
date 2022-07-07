use std::collections::HashMap;

use crate::ds::tile::Tile;
use crate::ds::sequence::Sequence;
use crate::ds::graph::Graph;
use crate::ds::graph::GraphProperties;

/// Build DeBruijn graph
pub struct DeBruijnBuilder {}

impl DeBruijnBuilder {

    pub fn new() -> DeBruijnBuilder {
        DeBruijnBuilder {}
    }

    pub fn build(&self, reads: &Tile) -> Graph<Sequence, usize>  {

        let mut graph = Graph::<Sequence, usize>::new(GraphProperties{directed: true});
        let mut present_nodes = HashMap::<Sequence, u64>::new();
        
        // Preallocate kmerid holders
        let mut kmer1_id = 0_u64;
        let mut kmer2_id = 0_u64;
    
        // Pre-calculate reverse complements of input reads
        let rev_reads = Tile::new();
        // reads.into_iter().for_each(|kmer|{ rev_reads.push(complement_dna(kmer.clone()))});
    
        // Construct deBruijn graph, iterate over all kmers and their complements 
        for kmer in reads.into_iter().chain(rev_reads.into_iter().as_ref()) {
    
            let kmer1 = Sequence::from(&kmer[0..kmer.len()-1]);
            let kmer2 = Sequence::from(&kmer[1..kmer.len()]);
    
            // Process Kmer1
            if present_nodes.contains_key(&kmer1) {
                kmer1_id = present_nodes[&kmer1];
            }
            else {
                kmer1_id = graph.add_node(kmer1.clone());
                present_nodes.insert(kmer1, kmer1_id);
            }
    
            // Process Kmer2
            if present_nodes.contains_key(&kmer2) {
                kmer2_id = present_nodes[&kmer2];
            }
            else {
                kmer2_id = graph.add_node(kmer2.clone());
                present_nodes.insert(kmer2, kmer2_id);
            }
    
            // Add edge connecting kmer1 with kmer2
            graph.add_edge(&kmer1_id, &kmer2_id, None).unwrap();
        }

        graph
    }
}

impl Default for DeBruijnBuilder {
    fn default() -> Self {
        Self::new()
    }
}
