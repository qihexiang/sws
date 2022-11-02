use petgraph::{graph::NodeIndex, Graph, Undirected};
use std::{str::FromStr};
use strum_macros::EnumString;

use crate::tokenizer::{
    smiles_tokenize, AROMATIC_ORGANIC_RE, BOND_RE, BRANCH_RE, NAGETIVE_RE, NOTHING_RE,
    ORGANIC_SUBSET_RE, RING_BOND_RE, STANDARD_NODE_RE,
};

#[derive(EnumString, Debug)]
enum Element {
    H = 1,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F,
    Ne,
    Na,
    Mg,
    Al,
    Si,
    P,
    S,
    Cl,
    Ar,
    K,
    Ca,
    Sc,
    Ti,
    V,
    Cr,
    Mn,
    Fe,
    Co,
    Ni,
    Cu,
    Zn,
    Ga,
    Ge,
    As,
    Se,
    Br,
    Kr,
}

#[derive(Debug)]
enum Chirality {
    Clockwise,
    Counter,
}

impl Chirality {
    pub fn new(s: &str) -> Option<Self> {
        match s {
            "@" => Some(Self::Counter),
            "@@" => Some(Self::Clockwise),
            _ => None,
        }
    }
}

#[derive(Debug)]
struct SmilesNode {
    element: Element,
    isotope: Option<u16>,
    charge: i8,
    chirality: Option<Chirality>,
    explicit_hydrogen: u8,
    selector: Option<String>,
    aromatic: bool,
}

impl SmilesNode {
    fn new(token: &str) -> Option<Self> {
        if let Some(captured) = ORGANIC_SUBSET_RE.captures(token) {
            let aromatic = AROMATIC_ORGANIC_RE.is_match(token);
            let element = captured
                .name("element")
                .and_then(|m| Some(m.as_str()))
                .and_then(|s| {
                    Some({
                        let mut first_capitalized = s[0..1].to_uppercase();
                        first_capitalized.push_str(&s[1..]);
                        first_capitalized
                    })
                })
                .and_then(|s| Some(Element::from_str(&s).unwrap()))
                .unwrap();
            let chirality = captured
                .name("chirality")
                .and_then(|m| Some(m.as_str()))
                .and_then(|s| Chirality::new(s));
            Some(SmilesNode {
                element,
                isotope: None,
                charge: 0,
                chirality,
                explicit_hydrogen: 0,
                selector: None,
                aromatic,
            })
        } else if let Some(captured) = STANDARD_NODE_RE.captures(token) {

            let element = captured
                .name("element")
                .and_then(|m| Some(m.as_str()))
                .and_then(|s| {
                    Some({
                        let mut first_capitalized = s[0..1].to_uppercase();
                        first_capitalized.push_str(&s[1..]);
                        first_capitalized
                    })
                })
                .and_then(|s| Some(Element::from_str(&s).unwrap()))
                .unwrap();
            let aromatic = AROMATIC_ORGANIC_RE.is_match(captured.name("element").unwrap().as_str());
            let chirality = captured
                .name("chirality")
                .and_then(|m| Some(m.as_str()))
                .and_then(|s| Chirality::new(s));
            let isotope: Option<u16> =
                captured
                    .name("isotope")
                    .and_then(|m| match m.as_str().parse() {
                        Ok(isotope) => Some(isotope),
                        Err(_) => None,
                    });
            let explicit_hydrogen: u8 = captured
                .name("explicit_hydrogen")
                .and_then(|m| match m.as_str().parse() {
                    Ok(isotope) => Some(isotope),
                    Err(_) => None,
                })
                .map_or(0, |v| v);
            let charge: i8 = {
                if let Some(charge) = captured
                    .name("charge_num")
                    .and_then(|m| Some(m.as_str().parse().unwrap()))
                {
                    charge
                } else if let Some(charge_str) =
                    captured.name("charge").and_then(|m| Some(m.as_str()))
                {
                    charge_str.len() as i8 * {
                        if NAGETIVE_RE.is_match(charge_str) {
                            -1
                        } else {
                            1
                        }
                    }
                } else {
                    0
                }
            };
            let selector = captured
                .name("selector")
                .and_then(|m| {
                    Some(
                        m.as_str()
                            .strip_prefix("{")
                            .unwrap()
                            .strip_suffix("}")
                            .unwrap(),
                    )
                })
                .map(String::from);
            Some(SmilesNode {
                element,
                isotope,
                charge,
                chirality,
                explicit_hydrogen,
                selector,
                aromatic,
            })
        } else {
            None
        }
    }
}

#[derive(PartialEq, Clone, Copy, Debug)]
enum BondType {
    Single(SingleBondType),
    Double,
    Triple,
    Quad,
    Aromatic,
    NoBond,
}

#[derive(PartialEq, Clone, Copy, Debug)]
enum SingleBondType {
    Normal,
    LeftUp,
    RightUp,
}

impl BondType {
    fn simple() -> Self {
        Self::Single(SingleBondType::Normal)
    }
    fn from_str(s: &str) -> Option<Self> {
        if BOND_RE.is_match(s) {
            Some(match s {
                "-" => Self::Single(SingleBondType::Normal),
                "/" => Self::Single(SingleBondType::RightUp),
                "\\" => Self::Single(SingleBondType::LeftUp),
                "=" => Self::Double,
                "#" => Self::Triple,
                "$" => Self::Quad,
                ":" => Self::Aromatic,
                _ => panic!("This shall never happend"),
            })
        } else {
            None
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

    fn next(&mut self, index: NodeIndex<u32>) {
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
            let id: u8 = captured.name("ring_id").unwrap().as_str().parse().unwrap();
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
            let target_ring = self.waiting_to_connect.get(target).unwrap();
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

#[derive(Debug)]
struct Molecule(Graph<SmilesNode, BondType, Undirected>);

impl Molecule {
    fn from_smiles(smiles: &str) -> Result<Self, String> {
        let mut graph = Graph::new_undirected();
        let mut construct_status = Status::new();
        let mut ring_status = RingStatus::new();
        let mut bond_to_connect: Option<BondType> = None;
        let tokens = smiles_tokenize(smiles);
        for token in tokens.into_iter() {
            if construct_status.is_none() {
                if let Some(node) = SmilesNode::new(token) {
                    construct_status.next(graph.add_node(node))
                } else {
                    return Err(format!(
                        "First token must be a SMILES atom token, but got {}",
                        token
                    ));
                }
            } else {
                let current_index = construct_status.get_index().unwrap();
                if let Some(node) = SmilesNode::new(token) {
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
        Ok(Molecule(graph))
    }
}

#[test]
fn generate_node() {
    let smiles_strs = vec![
        "c1cc[13c]cc1",
        "c1ccccc1C@@(N)(P)S",
        "c1ccccc1/C=C/(N)P",
        "c1ccccc1/C=[C{selected}]/(N)P",
        "N/C=[C{selected}]\\(N(C)[O-])P",
        "c1c([OH2+])cccc1[P{selected}]6(c2ccccc2)Cc3c(cccc3)C[PH2{selected}](c4ccccc4)(c5ccccc5)[Fe+2]6",
        "C=1CCC1",
        "C12C3C4C1C5C2C3C45"
    ];

    for smiles in smiles_strs {
        println!("{}", smiles);
        let mole = Molecule::from_smiles(smiles);
        if let Ok(mole) = mole {
            println!("{:?}", mole);
        }
    }
}
