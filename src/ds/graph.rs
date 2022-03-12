use hashbrown::{HashMap, HashSet};
use rand::{thread_rng, Rng};
use std::ops::{Index, IndexMut};
use std::hash::Hash;
use std::hash::Hasher;
use std::fs::File;
use std::io::Write;
use std::fmt;

#[derive(Clone, Debug, PartialEq)]
/// Graph operation error
pub enum GraphErr {
    /// There is no vertex with the given id in the graph
    NoSuchVertex,

    /// There is no such edge in the graph
    NoSuchEdge,

    /// Could not add an edge to the graph
    CannotAddEdge,

    /// The given weight is invalid
    InvalidWeight,

    /// The operation cannot be performed as it will
    /// create a cycle in the graph.
    CycleError,

    // #[cfg(feature = "dot")]
    // /// Could not render .dot file
    // CouldNotRender,

    // #[cfg(feature = "dot")]
    // /// The name of the graph is invalid. Check [this](https://docs.rs/dot/0.1.1/dot/struct.Id.html#method.new)
    // /// out for more information.
    // InvalidGraphName,
}

#[derive(Clone, Debug)]
pub struct GraphProperties { pub directed: bool}

#[derive(Clone, Debug)]
/// Edge internal struct
pub struct Edge<T: fmt::Display> {
    start: u64,
    end: u64,
    pub data: Option<T>,
}

#[derive(Clone, Debug)]
pub struct Node<T: fmt::Display> {
    id: u64,
    incoming: Vec<u64>,
    outgoing: Vec<u64>,
    pub data: T,
}

#[derive(Clone, Debug)]
/// Graph data-structure
pub struct Graph<N: fmt::Display, E: fmt::Display> {
    /// Mapping of vertex ids and vertex values
    nodes: HashMap<u64, Node<N>>,
    /// Mapping between edges and weights
    edges: HashMap<u64, Edge<E>>,
    /// Properties of the represented graph
    properties: GraphProperties,
}

impl<N: fmt::Display, E: fmt::Display + Clone> Graph<N, E> {

    pub fn new(props: GraphProperties) -> Graph<N, E> {
        Graph {
            nodes: HashMap::new(),
            edges: HashMap::new(),
            properties: props,
        }
    }

    pub fn add_node(&mut self, data: N) -> u64 {
        let mut rng = thread_rng();
        let nid = rng.gen_range(0..1000000000);
        self.nodes.insert(nid, Node{ id: nid, data: data, incoming: vec![], outgoing: vec![],});
        nid
    }   

    pub fn add_edge(&mut self, a: &u64, b: &u64, data: Option<E>) -> Result<(), GraphErr> {
        if self.has_edge(a, b) {
            return Ok(());
        }

        let id1 = if self.nodes.get(a).is_some() {
            *a
        } else {
            return Err(GraphErr::NoSuchVertex);
        };

        let id2 = if self.nodes.get(b).is_some() {
            *b
        } else {
            return Err(GraphErr::NoSuchVertex);
        };

        let mut rng = thread_rng();
        
        let mut eid = rng.gen_range(0..1000000000);
        self.edges.insert(eid, Edge{start: id1, end: id2, data: data.clone()});
        self.nodes.get_mut(&id1).unwrap().outgoing.push(eid);
        self.nodes.get_mut(&id2).unwrap().incoming.push(eid);

        if !self.properties.directed {
            eid = rng.gen_range(0..1000000000);
            self.edges.insert(eid, Edge{start: id2, end: id1, data: data.clone()});
            self.nodes.get_mut(&id2).unwrap().outgoing.push(eid);
            self.nodes.get_mut(&id1).unwrap().incoming.push(eid);
        }

        Ok(())
    }

    pub fn get_node(&self, id: &u64) -> &Node<N> {
        self.nodes.get(id).unwrap()
    }

    pub fn get_edge(&mut self, id: &u64) -> &mut Edge<E> {
        self.edges.get_mut(&id).unwrap()
    }    

    pub fn has_edge(&self, a: &u64, b: &u64) -> bool {
        match self.nodes.get(a) {
            Some(node) => node.outgoing.contains(b),
            None => false,
        }
    }    

    pub fn edges(&self) -> impl Iterator<Item=&u64> {
        self.edges.keys()
    }

    pub fn nodes(&self) -> impl Iterator<Item=&u64> {
        self.nodes.keys()
    }

    // pub fn nodes(&self, id: u64) -> &N {
    //     self.nodes.get(id).unwrap()
    // }

    // Obtain Iterator over neigbouring nodes or adjacent edges
    pub fn out_neighbors<'a>(&'a self, id: u64) -> Box<dyn Iterator<Item=&u64>+'a> {
        match self.nodes.get(&id) {
            Some(x) => Box::new(x.outgoing.iter().map(|idx| &self.edges.get(idx).unwrap().end)),
            _ => Box::new(::std::iter::empty())
        }
    }

    pub fn out_edges<'a>(&'a self, id: u64) -> Box<dyn Iterator<Item=&u64>+'a> {
        match self.nodes.get(&id) {
            Some(x) => Box::new(x.outgoing.iter()),
            _ => Box::new(::std::iter::empty())
        }
    }

    pub fn in_neighbors<'a>(&'a self, id: u64) -> Box<dyn Iterator<Item=&u64>+'a> {
        match self.nodes.get(&id) {
            Some(x) => Box::new(x.incoming.iter().map(|idx| &self.edges.get(idx).unwrap().end)),
            _ => Box::new(::std::iter::empty())
        }
    }

    pub fn in_edges<'a>(&'a self, id: u64) -> Box<dyn Iterator<Item=&u64>+'a> {
        match self.nodes.get(&id) {
            Some(x) => Box::new(x.incoming.iter()),
            _ => Box::new(::std::iter::empty())
        }
    }    

    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
    
    pub fn write_dot(&self, f: &str) {
        let mut w = File::create(f).unwrap();
        if self.properties.directed {
            writeln!(&mut w, "digraph sample {{").unwrap();
        }
        else {
            writeln!(&mut w, "strict graph sample {{").unwrap();
        }
        for (id, edge) in &self.edges {
            match (self.properties.directed, edge.data.as_ref()) {
                (true, Some(x)) =>  writeln!(&mut w, "  {} -> {} [ label = \"{}\" ];", edge.start, edge.end, x).unwrap(),
                (true, None) => writeln!(&mut w, "  {} -> {};", edge.start, edge.end).unwrap(),
                (false, Some(x)) => writeln!(&mut w, "  {} -- {} [ label = \"{}\" ];", edge.start, edge.end, x).unwrap(),
                (false, None) => writeln!(&mut w, "  {} -- {};", edge.start, edge.end).unwrap(),
            };  
        }
        for (id, node) in &self.nodes {
            writeln!(&mut w, "  {} [ label = \"{}\" ];", id, node.data).unwrap();
        }

        writeln!(&mut w, "}}").unwrap();
    }
}


/// DFS

pub struct Dfs<'a, N: fmt::Display, E: fmt::Display + Clone> {
    current: u64,
    stack: Vec<u64>,
    visited: HashSet<u64>,
    graph: &'a Graph<N,E>
}

impl<'a, N: fmt::Display, E: fmt::Display + Clone> Dfs<'a, N, E> {

    pub fn new(graph: &'a Graph<N, E>) -> Dfs<'a, N, E> {
        let mut stack = Vec::<u64>::new();
        let current = *graph.nodes().last().unwrap();
        stack.push(current);
        Dfs {
            current: current,
            stack: stack,
            visited: HashSet::new(),
            graph: graph
        }
    }

    pub fn init(&mut self, id: &u64) {
        self.current = *id;
        self.stack.clear();
        self.stack.push(*id);
    }

    pub fn process_node(&mut self) -> Result<u64, GraphErr> {
        // Obtain node from stack. Search is done, when stack empty.
        self.current = match self.stack.pop() {
            Some(x) => x,
            _ => return Err(GraphErr::NoSuchVertex)
        };
        for node in self.graph.out_neighbors(self.current) {
            if !self.visited.contains(node) {
                self.stack.push(*node);
            }
        }
        self.visited.insert(self.current);
        Ok(self.current)
    }
}


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
    suffix_start: usize,
    suffix_stop: i64
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
pub struct Ukonen<T:fmt::Display + Index<usize, Output=u8>> {
    // Sequence to be processed
    seq: T, 
    idx: usize,
    // State indicators
    remaining : usize,
    active_node_id : u64, 
    active_edge_id : i64,
    active_length : usize,
    end: u8,
    // Graph 
    root_id : u64,
    graph: Graph<UkonenNode, UkonenEdge>,

}

impl<T: fmt::Display + Clone + Index<usize, Output=u8>> Ukonen<T> {

    pub fn new(sequence: T) -> Ukonen<T> {

        let mut gr = Graph::<UkonenNode, UkonenEdge>::new(GraphProperties{directed: true});
        let root_id = gr.add_node(UkonenNode::new());

        Ukonen {
            seq : sequence, 
            idx : 0,
            remaining : 0,
            active_node_id : root_id, 
            active_edge_id : -1,
            active_length : 0,
            end : 0,
            root_id : root_id,
            graph :  gr
        }
    }    


    pub fn process(&mut self) -> &Graph<UkonenNode, UkonenEdge>{
        
        let value = self.seq[self.idx];
        let current_node  = self.graph.get_node(&self.active_node_id);

        // Increment remaining suffixes to be processed
        self.remaining += 1;
        // Set the current `end` character to new value
        self.end = value;

        if self.active_edge_id == -1 {
            if !current_node.data.suffix_edge_ids.contains_key(&value) {
                let nid = self.graph.add_node(UkonenNode::new());
                let suffix = UkonenEdge{suffix_start: self.idx, suffix_stop: -1};
                let eid = self.graph.add_edge(&self.active_node_id, &nid, Some(suffix));
                self.remaining -= 1;
            }
            else {
                self.active_edge_id = current_node.data.suffix_edge_ids[&value] as i64;
                self.active_length += 1;
            }
        }
        else {
            let active_edge = self.graph.get_edge(&(self.active_edge_id as u64));
            let lookup_idx = self.active_length + active_edge.data.as_ref().unwrap().suffix_start;
            
            if self.seq[lookup_idx] == value {
                self.active_length += 1;
            }
            else {
                
            }

            // if 
        }

        self.idx += 1;

        &self.graph

        // self.active_edge_id = current_node.data.suffix_edge_ids[&value];
        // let mut edge = self.graph.get_edge(&self.active_edge_id);
        // edge.data.as_mut().unwrap().start_idx = self.idx;
        // edge.data.as_mut().unwrap().stop_idx = -1;

    }

}