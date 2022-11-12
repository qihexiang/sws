use petgraph::{
    stable_graph::{EdgeIndex, NodeIndex},
    Direction::Incoming,
};

use crate::definitions::{
    atom::Atom,
    bond::{Bond, BondType},
};

use super::Workspace;

impl Workspace {
    /// reset the root of a structure.
    pub fn reset_root(&mut self, root_index: NodeIndex) -> Option<()> {
        self.reverse_recursive(root_index, None)
    }

    pub fn add_hydrogen_to_structure(&mut self, node: NodeIndex) -> Option<()> {
        let atoms = self.get_atoms_of_structure(node)?;
        for atom in atoms {
            self.add_hydrogen_to_atom(atom);
        }
        Some(())
    }

    /// Connect two existed atoms.
    pub fn connect(
        &mut self,
        outgoing_from: NodeIndex,
        incoming_to: NodeIndex,
        bond_type: BondType,
    ) -> EdgeIndex {
        self.graph.add_edge(
            outgoing_from,
            incoming_to,
            Bond::new(
                bond_type,
                self.in_same_structure(&[outgoing_from, incoming_to]),
            ),
        )
    }
}

/// Implement private functions used upon
impl Workspace {
    /// Find a chemical bond and reverse it direction in graph.
    /// This shall not influcence the true structure
    /// but the support of @ chirality is still WIP.
    /// If it's a ring-bond, this function won't reverse it and return `Some(false)`
    fn reverse_bond(&mut self, from: NodeIndex, to: NodeIndex) -> Option<bool> {
        let current_edge = self.graph.find_edge(from, to)?;
        let bond = self.graph.remove_edge(current_edge)?;
        self.graph.add_edge(to, from, bond.reverse());
        if bond.is_ring_bond() {
            Some(false)
        } else {
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

    /// Add a new atom and connect it to an existed atom.
    fn connect_new_atom(
        &mut self,
        atom: Atom,
        connect_to: NodeIndex,
        bond_type: BondType,
    ) -> NodeIndex {
        let new_node = self.graph.add_node(atom);
        self.connect(connect_to, new_node, bond_type);
        new_node
    }

    /// Add hydrogen atoms to an existed atom.
    fn add_hydrogen_to_atom(&mut self, atom: NodeIndex) -> Option<Vec<NodeIndex>> {
        let node = self.get_atom_mut(atom)?;
        let explicit_hydrogen = node.explicit_hydrogen;
        node.explicit_hydrogen = 0;
        let hydrogens_to_add = if explicit_hydrogen != 0 {
            explicit_hydrogen
        } else if node.element.default_hydrogen() == 0 {
            0
        } else {
            let default_hydrogens = node.element.default_hydrogen();
            let aromatic = if node.aromatic { 1 } else { 0 };
            let charge = node.charge;
            let neighbors = self
                .graph
                .neighbors_undirected(atom)
                .collect::<Vec<_>>()
                .len() as isize;
            let need_to_add = (default_hydrogens as isize) + charge - neighbors - aromatic;
            if need_to_add >= 0 {
                need_to_add as usize
            } else {
                0
            }
        };
        let mut added_hydrogens = vec![];
        while added_hydrogens.len() != hydrogens_to_add {
            added_hydrogens.push(self.connect_new_atom(
                Atom::new("[H]").unwrap(),
                atom,
                BondType::Single,
            ))
        }
        Some(added_hydrogens)
    }
}
