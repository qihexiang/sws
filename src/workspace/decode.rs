use crate::tokenizer::{sws_tokenize, BRANCH_RE, NOTHING_RE, RING_BOND_RE};
use petgraph::graph::NodeIndex;

use crate::definitions::{
    atom::Atom,
    bond::{Bond, BondType},
};

use super::Workspace;

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

    fn next(&mut self, index: NodeIndex) -> NodeIndex {
        self.current = Some(index);
        index
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
                .and_then(|s| Some(s.strip_prefix("%").unwrap_or(s)))
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

impl Workspace {
    /// add a SMILES into workspace as a structure.
    pub fn add_structure(&mut self, smiles: &str) -> Result<NodeIndex, String> {
        let mut construct_status = Status::new();
        let mut ring_status = RingStatus::new();
        let mut bond_to_connect: Option<BondType> = None;
        let tokens = sws_tokenize(smiles);
        let result = if let Some(node) = Atom::new(tokens[0]) {
            Ok(construct_status.next(self.graph.add_node(node)))
        } else {
            Err(format!(
                "First token must be a SMILES atom token, but got {}",
                tokens[0]
            ))
        };
        for token in tokens[1..].into_iter() {
            let current_index = construct_status.get_index().map_err(String::from)?;
            if let Some(node) = Atom::new(token) {
                let node_index = self.graph.add_node(node);
                self.graph.add_edge(
                    current_index,
                    node_index,
                    if let Some(bond_type) = bond_to_connect {
                        bond_to_connect = None;
                        Bond::new(bond_type, false)
                    } else {
                        let current_node = self
                            .graph
                            .node_weight(current_index)
                            .expect("Current node shall be inserted into graph before now.");
                        let next_node = self
                            .graph
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
                if let Some((previous_index, bond)) = ring {
                    if let Some(bond) = bond {
                        self.graph.add_edge(previous_index, current_index, bond);
                    } else {
                        self.graph.add_edge(previous_index, current_index, {
                            let current_node = self
                                .graph
                                .node_weight(current_index)
                                .expect("Current node shall be inserted into graph before now.");
                            let next_node = self
                                .graph
                                .node_weight(previous_index)
                                .expect("Next node is inserted into graph just now.");
                            if current_node.aromatic && next_node.aromatic {
                                Bond::new(BondType::Aromatic, true)
                            } else {
                                Bond::new(BondType::Single, true)
                            }
                        });
                    }
                }
            } else if BRANCH_RE.is_match(token) {
                match *token {
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

        if construct_status.branch.len() != 0 {
            Err(format!(
                "Uncleaned brnach stack. Some opened stack not cloused: {:?}",
                construct_status
                    .branch
                    .iter()
                    .map(|index| self.graph.node_weight(*index).unwrap())
                    .collect::<Vec<&Atom>>()
            ))
        } else if ring_status.waiting_to_connect.len() != 0 {
            Err(format!(
                "Some rings not closed: {:?}",
                ring_status.waiting_to_connect
            ))
        } else {
            result
        }
    }
}
