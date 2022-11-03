use petgraph::{graph::NodeIndex, Graph, Undirected};

use crate::tokenizer::{smiles_tokenize, NOTHING_RE, BRANCH_RE, RING_BOND_RE};

use super::{bond::BondType, smiles_node::SmilesNode};

#[derive(Debug)]
pub struct Molecule(Graph<SmilesNode, BondType, Undirected>);

impl Molecule {
    pub fn from_smiles(smiles: &str) -> Result<Self, String> {
        let mut graph = Graph::new_undirected();
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
                        if let Some(bond) = bond_to_connect {
                            bond_to_connect = None;
                            bond
                        } else {
                            let current_node = graph
                                .node_weight(current_index)
                                .expect("Current node shall be inserted into graph before now.");
                            let next_node = graph
                                .node_weight(node_index)
                                .expect("Next node is inserted into graph just now.");
                            if current_node.aromatic && next_node.aromatic {
                                BondType::Aromatic
                            } else {
                                BondType::simple()
                            }
                        },
                    );
                    construct_status.next(node_index);
                } else if let Some(bond) = BondType::from_str(token) {
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
                                    BondType::Aromatic
                                } else {
                                    BondType::simple()
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
    waiting_to_connect: Vec<(NodeIndex, Option<BondType>, u8)>,
}

impl RingStatus {
    fn identify_ring(token: &str) -> Option<(Option<BondType>, u8)> {
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
                .and_then(BondType::from_str);
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
        bond_type: Option<BondType>,
        id: u8,
    ) -> Option<(NodeIndex, Option<BondType>)> {
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