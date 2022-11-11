use petgraph::{graph::NodeIndex, stable_graph::EdgeIndex};

use super::workspace::Workspace;

struct SmilesGenerator<'a> {
    workspace: &'a Workspace,
    previous_node: Option<NodeIndex>,
    current_node: Option<NodeIndex>,
    branch_stack: Vec<(NodeIndex, Vec<NodeIndex>)>,
    ring_bonds: Vec<EdgeIndex>,
}

impl<'a> Iterator for SmilesGenerator<'a> {
    type Item = String;

    fn next(&mut self) -> Option<Self::Item> {
        let current_node = self.current_node?;
        let mut fragment = String::new();

        let atom = self
            .workspace
            .get_atom(current_node)
            .expect("should always find an atom");

        if let Some(previous_node) = self.previous_node {
            let (bond, _) = self.workspace.get_edge(previous_node, current_node).expect(
                format!(
                    "No connect between {:?} {:?}",
                    previous_node, self.current_node
                )
                .as_str(),
            );
            if !(self.workspace.get_atom(previous_node).unwrap().aromatic
                && atom.aromatic
                && bond.is_aromatic())
                && !(bond.is_normal_single())
            {
                fragment.push_str(bond.as_str())
            }
            let atom_token = atom.to_token();
            fragment.push_str(&atom_token);
        } else {
            fragment.push_str(&atom.to_token());
        }

        let all_nexts = self.workspace.get_next_neighbors(current_node);

        let mut atom_ring_bonds = self
            .workspace
            .graph
            .neighbors_undirected(current_node)
            .filter(|node_index| {
                let (bond, _, _) = self
                    .workspace
                    .get_edge_undirected(current_node, *node_index)
                    .unwrap();
                bond.is_ring_bond()
            })
            .collect::<Vec<_>>();

        while let Some(ring_bond_neighbor) = atom_ring_bonds.pop() {
            let (bond, edge, _) = self
                .workspace
                .get_edge_undirected(current_node, ring_bond_neighbor)
                .unwrap();

            let ring_id = if let Some(position) = self
                .ring_bonds
                .iter()
                .position(|edge_index| edge == *edge_index)
            {
                position + 1
            } else {
                self.ring_bonds.push(edge);
                self.ring_bonds.len()
            };

            let bond_token = if !(self
                .workspace
                .get_atom(ring_bond_neighbor)
                .unwrap()
                .aromatic
                && atom.aromatic
                && bond.is_aromatic())
                && !(bond.is_normal_single())
                && !bond.is_normal_single()
                || self.ring_bonds.contains(&edge)
            {
                ""
            } else {
                bond.as_str()
            };

            let ring_id_token = if ring_id >= 10 {
                format!("%{}", ring_id)
            } else {
                ring_id.to_string()
            };

            let full_token = if atom_ring_bonds.len() >= 1 && bond_token != "" {
                format!("({}{})", bond_token, ring_id_token)
            } else {
                format!("{}{}", bond_token, ring_id_token)
            };

            fragment.push_str(&full_token);
        }

        let mut nexts = all_nexts
            .iter()
            .copied()
            .filter(|neighbor_index| {
                let (bond, _) = self
                    .workspace
                    .get_edge(current_node, *neighbor_index)
                    .unwrap();
                !bond.is_ring_bond()
            })
            .collect::<Vec<_>>();

        self.current_node = if nexts.len() > 0 {
            self.previous_node = self.current_node;
            if nexts.len() > 1 {
                fragment.push_str("(");
                let next = nexts.pop().unwrap();
                self.branch_stack.push((self.previous_node.unwrap(), nexts));
                Some(next)
            } else {
                nexts.pop()
            }
        } else {
            if let Some((prev, mut branches)) = self.branch_stack.pop() {
                self.previous_node = Some(prev);
                fragment.push_str(")");
                if branches.len() > 1 {
                    fragment.push_str("(");
                }
                let next = branches.pop().unwrap();
                if branches.len() >= 1 {
                    self.branch_stack.push((prev, branches));
                }
                Some(next)
            } else {
                None
            }
        };
        Some(fragment)
    }
}

impl<'a> SmilesGenerator<'a> {
    fn new(workspace: &'a Workspace, structure_node: NodeIndex) -> Option<Self> {
        Some(Self {
            workspace,
            previous_node: None,
            current_node: workspace.find_root_of(structure_node),
            branch_stack: vec![],
            ring_bonds: vec![],
        })
    }
}

impl Workspace {
    pub fn to_smiles(&self, node: NodeIndex) -> Option<String> {
        let generator = SmilesGenerator::new(self, node)?;
        let mut smiles = String::new();
        for fragment in generator {
            smiles.push_str(&fragment);
        }
        Some(smiles)
    }
}
