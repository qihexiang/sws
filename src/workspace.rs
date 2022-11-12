pub mod encode;
pub mod decode;
pub mod editor;
pub mod accessor;

use petgraph::{dot::Dot, stable_graph::{StableGraph, NodeIndex}, Directed};

use crate::definitions::{bond::Bond, atom::Atom};


/// a Workspace is a graph space that can deal with structures
#[derive(Debug)]
pub struct Workspace {
    pub graph: StableGraph<Atom, Bond, Directed>,
}

impl Workspace {
    /// create a workspace
    pub fn new() -> Self {
        Self {
            graph: StableGraph::new(),
        }
    }
    pub fn dot_representation(&self) -> Dot<&StableGraph<Atom, Bond, Directed>> {
        Dot::new(&self.graph)
    }
    /// Check if some atoms are in the same structure(has same structure root)
    pub fn in_same_structure(&self, atoms: &[NodeIndex]) -> bool {
        let mut roots = atoms
            .iter()
            .map(|node_index| self.find_root_of(*node_index))
            .collect::<Vec<_>>();
        if let Some(root_index) = roots.pop() {
            if let Some(_) = roots
                .into_iter()
                .find(|another_root| another_root != &root_index)
            {
                false
            } else {
                true
            }
        } else {
            false
        }
    }
}
