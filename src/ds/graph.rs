use hashbrown::{HashMap, HashSet};
use rand::{thread_rng, Rng};
use std::hash::Hasher;
use std::hash::Hash;
use std::io::Write;
use std::fs::File;
use std::fmt;
use std::io;

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
    pub start: u64,
    pub end: u64,
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
    /// Id of root node
    root: Option<u64>,
}

impl<N: fmt::Display, E: fmt::Display + Clone> Graph<N, E> {

    pub fn new(props: GraphProperties) -> Graph<N, E> {
        Graph {
            nodes: HashMap::new(),
            edges: HashMap::new(),
            properties: props,
            root: None,
        }
    }

    pub fn add_node(&mut self, data: N) -> u64 {
        let mut rng = thread_rng();
        // TODO: Make sure that indizes are unique
        let nid = rng.gen_range(0..10000000000); 
        self.nodes.insert(nid, Node{ id: nid, data: data, incoming: vec![], outgoing: vec![],});
        nid
    }   

    pub fn add_edge(&mut self, a: &u64, b: &u64, data: Option<E>) -> Result<u64, GraphErr> {
        
        let e = self.has_edge(a, b);
        if  e.is_ok() {
            return e;
        };

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
        
        // TODO: Make sure that indizes are unique
        let mut rng = thread_rng();
        let mut eid = rng.gen_range(0..10000000000);
        self.edges.insert(eid, Edge{start: id1, end: id2, data: data.clone()});
        self.nodes.get_mut(&id1).unwrap().outgoing.push(eid);
        self.nodes.get_mut(&id2).unwrap().incoming.push(eid);

        Ok(eid)
    }

    pub fn get_node(&self, id: &u64) -> & Node<N> {
        self.nodes.get(id).unwrap()
    }

    pub fn get_node_mut(&mut self, id: &u64) -> &mut Node<N> {
        self.nodes.get_mut(id).unwrap()
    }

    pub fn get_root(&self) -> Option<u64> {
        self.root
    }

    pub fn set_root(&mut self, id: u64) {
        self.root = Some(id);
    }

    pub fn get_edge(&self, id: &u64) -> & Edge<E> {
        // println!("ID: id");
        self.edges.get(id).unwrap()
    }    

    pub fn get_edge_mut(&mut self, id: &u64) -> &mut Edge<E> {
        self.edges.get_mut(id).unwrap()
    }

    pub fn has_edge(&self, a: &u64, b: &u64) -> Result<u64, GraphErr> {
        let node = self.nodes.get(a);
        let edge_id = match self.nodes.get(a) {
            Some(x) => x.outgoing.iter().find(|idx| self.edges.get(idx).unwrap().end == *b),
            _ => None
        };
        let result = match edge_id {
            Some(id) => Ok(*id),
            _ => Err(GraphErr::NoSuchEdge)
        };
        result
    }    

    pub fn remove_edge(&mut self, id: &u64) -> Option<Edge<E>> {

        // Check if edge exists in graph
        if !self.edges.contains_key(id) {
            return None;
        }

        // Remove outgoing adjacency in corresponding node
        let out_node_id = self.get_edge(id).start;
        let out_node_neighbours = &mut self.get_node_mut(&out_node_id).outgoing;
        let index_out = out_node_neighbours.iter().position(|x| x == id).unwrap();
        out_node_neighbours.remove(index_out);  

        // Remove incoming adjacency in corresponding node
        let in_node_id = self.get_edge(id).end;
        let in_node_neighbours = &mut self.get_node_mut(&in_node_id).incoming;
        let index_in = in_node_neighbours.iter().position(|x| x == id).unwrap();
        in_node_neighbours.remove(index_in);        

        // Remove edge from storage container
        self.edges.remove(id)
    }    

    pub fn edges(&self) -> impl Iterator<Item=&u64> {
        self.edges.keys()
    }

    pub fn nodes(&self) -> impl Iterator<Item=&u64> {
        self.nodes.keys()
    }
    
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
                (true, Some(data)) =>  writeln!(&mut w, "  {} -> {} [ label = \"{}\" ];", edge.start, edge.end, data).unwrap(),
                (true, None) => writeln!(&mut w, "  {} -> {};", edge.start, edge.end).unwrap(),
                (false, Some(data)) => writeln!(&mut w, "  {} -- {} [ label = \"{}\" ];", edge.start, edge.end, data).unwrap(),
                (false, None) => writeln!(&mut w, "  {} -- {};", edge.start, edge.end).unwrap(),
            };  
        }
        for (id, node) in &self.nodes {
            writeln!(&mut w, "  {} [ label = \"{}\" ];", id, node.data).unwrap();
        }

        writeln!(&mut w, "}}").unwrap();
    }
}
