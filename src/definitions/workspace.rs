use petgraph::{
    dot::Dot, graph::NodeIndex, stable_graph::StableGraph, Directed, Direction::Incoming,
};

use super::{
    bond::{Bond, BondType},
    element::Element,
    atom::Atom,
};

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

    /// find the first node of each structures in the workspace
    pub fn find_strcutres_roots(&self) -> Vec<NodeIndex> {
        self.graph
            .node_indices()
            .filter(|index| {
                // 找到该原子对应的所有进入方向的键，但是排除成环键
                let node_incoming_amount = self
                    .graph
                    .edges_directed(*index, Incoming)
                    .filter(|edge| {
                        let bond = edge.weight();
                        !bond.is_ring_bond()
                    })
                    .collect::<Vec<_>>()
                    .len();
                if node_incoming_amount == 0 {
                    true
                } else {
                    false
                }
            })
            .collect()
    }

    /// Find a chemical bond and reverse it direction in graph.
    /// This shall not influcence the true structure
    /// but the support of @ chirality is still WIP.
    /// If it's a ring-bond, this function won't reverse it and return `Some(false)`
    fn reverse_bond(&mut self, from: NodeIndex, to: NodeIndex) -> Option<bool> {
        let current_edge = self.graph.find_edge(from, to)?;
        let bond = self.graph.remove_edge(current_edge)?;
        if bond.is_ring_bond() {
            self.graph.add_edge(from, to, bond);
            Some(false)
        } else {
            self.graph.add_edge(to, from, bond.reverse());
            Some(true)
        }
    }

    /// Reverse the direction from the given node recurvisely.
    /// Recurse will stop if there is no incoming edge or found a ring-key.
    fn reverse_recursive(&mut self, current: NodeIndex, last: Option<NodeIndex>) -> Option<()> {
        let incoming_edges = self
            .graph
            .neighbors_directed(current, Incoming)
            .collect::<Vec<_>>();
        for neighbor in incoming_edges {
            if !(Some(neighbor) == last) {
                // 当reverse_bond函数返回true时，这个键没有形成环，则向下递归
                if self.reverse_bond(neighbor, current)? {
                    self.reverse_recursive(neighbor, Some(current))?;
                }
            }
        }
        Some(())
    }

    /// reset the root of a structure.
    pub fn reset_root(&mut self, root_index: NodeIndex) -> Option<()> {
        self.reverse_recursive(root_index, None)
    }

    /// Find the root of given node
    pub fn find_root_of(&self, node_index: NodeIndex) -> Option<NodeIndex> {
        let incoming_edges = self
            .graph
            .edges_directed(node_index, Incoming)
            .filter(|edge| !edge.weight().is_ring_bond())
            .collect::<Vec<_>>()
            .len();
        if incoming_edges == 0 {
            Some(node_index)
        } else {
            let neighbors = self
                .graph
                .neighbors_directed(node_index, Incoming)
                .filter(|neighbor_index| {
                    self.graph
                        .find_edge(*neighbor_index, node_index)
                        .and_then(|edge_index| self.graph.edge_weight(edge_index))
                        .and_then(|bond| Some(!bond.is_ring_bond()))
                        .unwrap_or(false)
                })
                .collect::<Vec<_>>();
            let incoming = neighbors.get(0)?;
            self.find_root_of(*incoming)
        }
    }

    /// Find nodes by given conditions.
    pub fn filter_nodes(&self, find_fn: &dyn Fn(&Atom) -> bool) -> Vec<NodeIndex> {
        self.graph
            .node_indices()
            .filter(|index| find_fn(&self.graph[*index]))
            .collect()
    }

    /// Find first node match given condition
    pub fn find_node(&self, find_fn: &dyn Fn(&Atom) -> bool) -> Option<NodeIndex> {
        self.filter_nodes(find_fn).get(0).copied()
    }

    pub fn get_atom(&self, index: NodeIndex) -> Option<&Atom> {
        self.graph.node_weight(index)
    }

    pub fn get_atom_mut(&mut self, index: NodeIndex) -> Option<&mut Atom> {
        self.graph.node_weight_mut(index)
    }

    pub fn connect_new_atom(&mut self, atom: Atom, connect_to: NodeIndex, bond_type: Bond) {
        let new_node = self.graph.add_node(atom);
        self.graph.add_edge(connect_to, new_node, bond_type);
    }

    pub fn add_hydrogens(&mut self) {
        self.filter_nodes(&|_| true).iter().for_each(|index| {
            let node = self.get_atom(*index).unwrap();
            let hydrogens_to_add = if node.explicit_hydrogen != 0 {
                node.explicit_hydrogen
            } else if node.element.default_hydrogen() == 0 {
                0
            } else {
                let default_hydrogens = node.element.default_hydrogen();
                let aromatic = if node.aromatic { 1 } else { 0 };
                let charge = node.charge;
                let neighbors = self
                    .graph
                    .neighbors_undirected(*index)
                    .collect::<Vec<_>>()
                    .len() as isize;
                let need_to_add = (default_hydrogens as isize) + charge - neighbors - aromatic;
                if need_to_add >= 0 {
                    need_to_add as usize
                } else {
                    0
                }
            };
            let mut counter = 0usize;
            loop {
                if counter == hydrogens_to_add {
                    break;
                }
                counter += 1;
                self.connect_new_atom(
                    Atom::new("[H]", None).unwrap(),
                    *index,
                    Bond::new(BondType::Single, false),
                )
            }
        });
    }

    pub fn remove_hydrogens(&mut self) {
        self.filter_nodes(&|node| node.element == Element::H)
            .iter()
            .for_each(|index| {
                self.graph.remove_node(*index);
            })
    }

    pub fn dot_representation(&self) -> Dot<&StableGraph<Atom, Bond, Directed>> {
        Dot::new(&self.graph)
    }
}

