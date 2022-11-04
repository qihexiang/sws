use petgraph::{dot::Dot, graph::NodeIndex, Directed, stable_graph::StableGraph};

use crate::tokenizer::{smiles_tokenize, BRANCH_RE, NOTHING_RE, RING_BOND_RE};

use super::{
    bond::{Bond, BondType},
    element::Element,
    smiles_node::SmilesNode,
};

#[derive(Debug)]
pub struct Molecule(StableGraph<SmilesNode, Bond, Directed>);

impl Molecule {
    pub fn from_smiles(smiles: &str) -> Result<Self, String> {
        let mut graph: StableGraph<SmilesNode, Bond> = StableGraph::new();
        let mut construct_status = Status::new();
        let mut ring_status = RingStatus::new();
        let mut bond_to_connect: Option<BondType> = None;
        let tokens = smiles_tokenize(smiles);
        for token in tokens.into_iter() {
            if construct_status.is_none() {
                if let Some(node) = SmilesNode::new(token, None) {
                    construct_status.next(graph.add_node(node))
                } else {
                    return Err(format!(
                        "First token must be a SMILES atom token, but got {}",
                        token
                    ));
                }
            } else {
                let current_index = construct_status.get_index().map_err(String::from)?;
                if let Some(node) = SmilesNode::new(token, Some(current_index)) {
                    let node_index = graph.add_node(node);
                    graph.add_edge(
                        current_index,
                        node_index,
                        if let Some(bond_type) = bond_to_connect {
                            bond_to_connect = None;
                            Bond::new(bond_type, false)
                        } else {
                            let current_node = graph
                                .node_weight(current_index)
                                .expect("Current node shall be inserted into graph before now.");
                            let next_node = graph
                                .node_weight(node_index)
                                .expect("Next node is inserted into graph just now.");
                            if current_node.aromatic && next_node.aromatic {
                                Bond::new(BondType::Aromatic, false)
                            } else {
                                Bond::new(BondType::Single, false)
                            }
                        },
                    );
                    construct_status.next(node_index);
                } else if let Some(bond) = BondType::new(token) {
                    bond_to_connect = Some(bond)
                } else if NOTHING_RE.is_match(token) {
                    bond_to_connect = Some(BondType::NoBond)
                } else if let Some((bond_type, id)) = RingStatus::identify_ring(token) {
                    let ring = ring_status.ring(current_index, bond_type, id);
                    if let Some((another_index, bond_type)) = ring {
                        if let Some(bond_type) = bond_type {
                            graph.add_edge(current_index, another_index, bond_type);
                        } else {
                            graph.add_edge(current_index, another_index, {
                                let current_node = graph.node_weight(current_index).expect(
                                    "Current node shall be inserted into graph before now.",
                                );
                                let next_node = graph
                                    .node_weight(another_index)
                                    .expect("Next node is inserted into graph just now.");
                                if current_node.aromatic && next_node.aromatic {
                                    Bond::new(BondType::Aromatic, false)
                                } else {
                                    Bond::new(BondType::Single, false)
                                }
                            });
                        }
                    }
                } else if BRANCH_RE.is_match(token) {
                    match token {
                        "(" => {
                            construct_status.enter_branch();
                        }
                        ")" => {
                            construct_status.quit_branch();
                        }
                        _ => {
                            panic!("Unknown operator catched {}", token)
                        }
                    }
                }
            }
        }

        graph.retain_edges(|graph, edge| {
            let bond_type = graph.edge_weight(edge).unwrap();
            !bond_type.is_no_bond()
        });

        if construct_status.branch.len() != 0 {
            Err(format!(
                "Uncleaned brnach stack. Some opened stack not cloused: {:?}",
                construct_status
                    .branch
                    .iter()
                    .map(|index| graph.node_weight(*index).unwrap())
                    .collect::<Vec<&SmilesNode>>()
            ))
        } else if ring_status.waiting_to_connect.len() != 0 {
            Err(format!(
                "Some rings not closed: {:?}",
                ring_status.waiting_to_connect
            ))
        } else {
            Ok(Molecule(graph))
        }
    }

    pub fn filter_nodes<F>(&self, find_fn: F) -> Vec<NodeIndex>
    where
        F: Fn(&SmilesNode) -> bool,
    {
        self.0
            .node_indices()
            .filter(|index| find_fn(&self.0[*index]))
            .collect()
    }

    pub fn find_node<F>(&self, find_fn: F) -> Option<NodeIndex>
    where
        F: Fn(&SmilesNode) -> bool,
    {
        self.0
            .node_indices()
            .filter(|index| find_fn(&self.0[*index]))
            .collect::<Vec<NodeIndex>>()
            .get(0)
            .and_then(|index| Some(*index))
    }

    pub fn get_atom(&self, index: NodeIndex) -> Option<&SmilesNode> {
        self.0.node_weight(index)
    }

    pub fn get_atom_mut(&mut self, index: NodeIndex) -> Option<&mut SmilesNode> {
        self.0.node_weight_mut(index)
    }

    pub fn connect_new_atom(&mut self, atom: SmilesNode, connect_to: NodeIndex, bond_type: Bond) {
        let new_node = self.0.add_node(atom);
        self.0.add_edge(connect_to, new_node, bond_type);
    }

    pub fn add_hydrogens(&mut self) {
        self.filter_nodes(|_| true).iter().for_each(|index| {
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
                    .0
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
                    SmilesNode::new("[H]", None).unwrap(),
                    *index,
                    Bond::new(BondType::Single, false),
                )
            }
        });
    }

    pub fn remove_hydrogens(&mut self) {
        let hydrogens = self.filter_nodes(|node| node.element == Element::H);
        println!("{:?}", hydrogens);
        for atom in hydrogens {
            match self.0.remove_node(atom) {
                Some(node) => {println!("{} removed", node)}
                None => {println!("Failed to remove {:?}", atom)}
            }
        }
    }

    pub fn dot_representation(&self) -> Dot<&StableGraph<SmilesNode, Bond, Directed>> {
        Dot::new(&self.0)
    }
}

struct Status {
    branch: Vec<NodeIndex>,
    current: Option<NodeIndex>,
}

impl Status {
    fn new() -> Self {
        Self {
            branch: vec![],
            current: None,
        }
    }

    fn next(&mut self, index: NodeIndex) {
        self.current = Some(index);
    }

    fn enter_branch(&mut self) {
        self.branch.push(
            self.current
                .expect("At least one atom existed when enter a branch"),
        );
    }

    fn quit_branch(&mut self) {
        if self.branch.len() > 0 {
            self.current = self.branch.pop();
        } else {
            panic!("Quit a branch when there in fact on branches")
        }
    }

    fn get_index(&self) -> Result<NodeIndex, &'static str> {
        if let Some(index) = self.current {
            Ok(index)
        } else {
            Err("Can't get index from stop status.")
        }
    }

    fn is_none(&self) -> bool {
        if let Some(_) = self.current {
            false
        } else {
            true
        }
    }
}

struct RingStatus {
    waiting_to_connect: Vec<(NodeIndex, Option<Bond>, u8)>,
}

impl RingStatus {
    fn identify_ring(token: &str) -> Option<(Option<Bond>, u8)> {
        if let Some(captured) = RING_BOND_RE.captures(token) {
            let id = captured
                .name("ring_id")
                .and_then(|m| Some(m.as_str()))
                .and_then(|s| {
                    Some(
                        s.parse::<u8>()
                            .expect("Shall not return error because regex matched"),
                    )
                })
                .unwrap();
            let bond_type = captured
                .name("bond_type")
                .and_then(|m| Some(m.as_str()))
                .and_then(|s| BondType::new(s))
                .and_then(|bond_type| Some(Bond::new(bond_type, true)));
            Some((bond_type, id))
        } else {
            None
        }
    }

    fn new() -> Self {
        Self {
            waiting_to_connect: vec![],
        }
    }

    fn ring(
        &mut self,
        node_index: NodeIndex,
        bond_type: Option<Bond>,
        id: u8,
    ) -> Option<(NodeIndex, Option<Bond>)> {
        if let Some(target) = self.waiting_to_connect.iter().position(|item| item.2 == id) {
            let target_ring = self
                .waiting_to_connect
                .get(target)
                .expect("Never get None as target index is found by position method.");
            if let Some(given_bond_type) = bond_type {
                if let Some(target_bond_type) = target_ring.1 {
                    if given_bond_type == target_bond_type {
                        let node_index = target_ring.0;
                        self.waiting_to_connect.remove(target);
                        Some((node_index, Some(target_bond_type)))
                    } else {
                        None
                    }
                } else {
                    let node_index = target_ring.0;
                    self.waiting_to_connect.remove(target);
                    Some((node_index, Some(given_bond_type)))
                }
            } else {
                if let Some(target_bond_type) = target_ring.1 {
                    let node_index = target_ring.0;
                    self.waiting_to_connect.remove(target);
                    Some((node_index, Some(target_bond_type)))
                } else {
                    let node_index = target_ring.0;
                    self.waiting_to_connect.remove(target);
                    Some((node_index, None))
                }
            }
        } else {
            self.waiting_to_connect.push((node_index, bond_type, id));
            None
        }
    }
}
