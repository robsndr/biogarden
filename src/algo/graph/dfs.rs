use std::fmt;
use std::collections::HashMap;
use std::collections::HashSet;

use crate::ds::graph::Graph;
use crate::ds::graph::GraphErr;

pub struct Dfs<'a, N: fmt::Display + Clone, E: fmt::Display + Clone> {
    current: u64,
    stack: Vec<u64>,
    visited: HashSet<u64>,
    graph: &'a Graph<N,E>
}

impl<'a, N: fmt::Display + Clone, E: fmt::Display + Clone> Dfs<'a, N, E> {

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
