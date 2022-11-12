use petgraph::{
    stable_graph::{EdgeIndex, NodeIndex},
    Direction::{self, Incoming},
};

use crate::definitions::{atom::Atom, bond::Bond};
use super::Workspace;

/// Here implements function that can get indexes
impl Workspace {
    /// Get neighbor NodeIndex from outgoing edges.
    pub fn get_next_neighbors(&self, node: NodeIndex) -> Vec<NodeIndex> {
        self.graph.neighbors(node).collect()
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

    /// Find the root of structure by given node
    pub fn find_root_of(&self, structure_node: NodeIndex) -> Option<NodeIndex> {
        let incoming_edges = self
            .graph
            .edges_directed(structure_node, Incoming)
            .filter(|edge| !edge.weight().is_ring_bond())
            .collect::<Vec<_>>()
            .len();
        if incoming_edges == 0 {
            Some(structure_node)
        } else {
            let neighbors = self
                .graph
                .neighbors_directed(structure_node, Incoming)
                .filter(|neighbor_index| {
                    self.graph
                        .find_edge(*neighbor_index, structure_node)
                        .and_then(|edge_index| self.graph.edge_weight(edge_index))
                        .and_then(|bond| Some(!bond.is_ring_bond()))
                        .unwrap_or(false)
                })
                .collect::<Vec<_>>();
            let incoming = neighbors.get(0)?;
            self.find_root_of(*incoming)
        }
    }

    /// Get all atoms' NodeIndex into one Vector.
    pub fn get_atoms_of_structure(&self, structure_node: NodeIndex) -> Option<Vec<NodeIndex>> {
        let root = self.find_root_of(structure_node)?;
        Some([vec![root], self.search_all_atoms(root)].concat())
    }

    /// Find nodes by given conditions.
    pub fn filter_nodes_in_structure<F>(
        &self,
        structre_root: NodeIndex,
        find_fn: F,
    ) -> Option<Vec<NodeIndex>>
    where
        F: Fn(&Atom) -> bool,
    {
        Some(
            self.get_atoms_of_structure(structre_root)?
                .into_iter()
                .filter(|index| find_fn(&self.graph[*index]))
                .collect(),
        )
    }

    /// Find first node match given condition
    pub fn find_node_in_structure<F>(
        &self,
        structure_root: NodeIndex,
        find_fn: F,
    ) -> Option<NodeIndex>
    where
        F: Fn(&Atom) -> bool,
    {
        self.filter_nodes_in_structure(structure_root, find_fn)?
            .get(0)
            .copied()
    }

    /// Find nodes by given conditions.
    pub fn filter_nodes<F>(&self, find_fn: F) -> Vec<NodeIndex>
    where
        F: Fn(&Atom) -> bool,
    {
        self.graph
            .node_indices()
            .filter(|index| find_fn(&self.graph[*index]))
            .collect()
    }

    /// Find first node match given condition
    pub fn find_node<F>(&self, find_fn: F) -> Option<NodeIndex>
    where
        F: Fn(&Atom) -> bool,
    {
        self.filter_nodes(find_fn).get(0).copied()
    }
}

/// Here implements function to get instances like Atom and bond
impl Workspace {
    pub fn get_atom(&self, index: NodeIndex) -> Option<&Atom> {
        self.graph.node_weight(index)
    }

    pub fn get_atom_mut(&mut self, index: NodeIndex) -> Option<&mut Atom> {
        self.graph.node_weight_mut(index)
    }

    pub fn get_edge_undirected(
        &self,
        atom_a: NodeIndex,
        atom_to: NodeIndex,
    ) -> Option<(&Bond, EdgeIndex, Direction)> {
        self.graph
            .find_edge_undirected(atom_a, atom_to)
            .and_then(|(edge, direction)| {
                self.graph
                    .edge_weight(edge)
                    .and_then(|bond| Some((bond, edge, direction)))
            })
    }

    pub fn get_edge(&self, atom_from: NodeIndex, atom_to: NodeIndex) -> Option<(&Bond, EdgeIndex)> {
        self.graph.find_edge(atom_from, atom_to).and_then(|edge| {
            self.graph
                .edge_weight(edge)
                .and_then(|bond| Some((bond, edge)))
        })
    }

    pub fn get_edge_mut(
        &mut self,
        atom_from: NodeIndex,
        atom_to: NodeIndex,
    ) -> Option<(&mut Bond, EdgeIndex)> {
        self.graph.find_edge(atom_from, atom_to).and_then(|edge| {
            self.graph
                .edge_weight_mut(edge)
                .and_then(|bond| Some((bond, edge)))
        })
    }
}

/// Private functions used upon
impl Workspace {
    /// Search all atoms connect by bond in same structure
    fn search_all_atoms(&self, node: NodeIndex) -> Vec<NodeIndex> {
        let neighbors = self
            .get_next_neighbors(node)
            .iter()
            .copied()
            .filter(|neighbor| {
                let edge_index = self.graph.find_edge(node, *neighbor);
                if let Some(edge_index) = edge_index {
                    let edge = self.graph.edge_weight(edge_index);
                    if let Some(bond) = edge {
                        !bond.is_ring_bond()
                    } else {
                        false
                    }
                } else {
                    false
                }
            })
            .collect::<Vec<_>>();
        let recursive_neighbors = neighbors
            .iter()
            .map(|neighbor| self.search_all_atoms(*neighbor))
            .collect::<Vec<_>>()
            .concat();
        [neighbors, recursive_neighbors].concat()
    }
}
