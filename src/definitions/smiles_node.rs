use std::{str::FromStr, fmt::Display};

use super::{chirality::ChiralityType, element::Element};
use crate::tokenizer::{AROMATIC_ORGANIC_RE, NAGETIVE_RE, ORGANIC_SUBSET_RE, STANDARD_NODE_RE};
use petgraph::graph::NodeIndex;

#[derive(Debug)]
pub struct SmilesNode {
    pub element: Element,
    pub isotope: Option<u16>,
    pub charge: isize,
    pub chirality: (Option<ChiralityType>, Option<NodeIndex>),
    pub explicit_hydrogen: usize,
    pub selector: Option<String>,
    pub aromatic: bool,
    pub react_id: Option<usize>,
}

impl SmilesNode {
    pub fn new(token: &str, current_node_index: Option<NodeIndex>) -> Option<Self> {
        if let Some(captured) = ORGANIC_SUBSET_RE.captures(token) {
            let (element, aromatic, chirality_type) = Self::minimal_node_info(&captured)?;
            Some(SmilesNode {
                element,
                isotope: None,
                charge: 0,
                chirality: (chirality_type, current_node_index),
                explicit_hydrogen: 0,
                selector: None,
                aromatic,
                react_id: None,
            })
        } else if let Some(captured) = STANDARD_NODE_RE.captures(token) {
            let (element, aromatic, chirality_type) = Self::minimal_node_info(&captured)?;
            let isotope: Option<u16> =
                captured
                    .name("isotope")
                    .and_then(|m| match m.as_str().parse() {
                        Ok(isotope) => Some(isotope),
                        Err(_) => None,
                    });
            let explicit_hydrogen = captured
                .name("explicit_hydrogen")
                .and_then(|m| match m.as_str().parse() {
                    Ok(isotope) => Some(isotope),
                    Err(_) => None,
                })
                .map_or(0, |v| v);
            let charge = {
                if let Some(charge) = captured
                    .name("charge_num")
                    .and_then(|m| Some(m.as_str().parse::<isize>().unwrap()))
                {
                    charge
                } else if let Some(charge_str) =
                    captured.name("charge").and_then(|m| Some(m.as_str()))
                {
                    charge_str.len() as isize * {
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
            let react_id = captured
                .name("react_id")
                .and_then(|m| Some(m.as_str()))
                .and_then(|s| Some(s.parse::<usize>().unwrap()));
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
                chirality: (chirality_type, current_node_index),
                explicit_hydrogen,
                selector,
                aromatic,
                react_id,
            })
        } else {
            None
        }
    }

    fn minimal_node_info(
        captured: &regex::Captures,
    ) -> Option<(Element, bool, Option<ChiralityType>)> {
        let element = captured.name("element").and_then(|m| Some(m.as_str()))?;
        let aromatic = AROMATIC_ORGANIC_RE.is_match(element);
        let chirality_type = captured
            .name("chirality")
            .and_then(|m| Some(m.as_str()))
            .and_then(ChiralityType::new);
        let mut capitalized = element[0..1].to_uppercase();
        capitalized.push_str(&element[1..]);
        let element = Element::from_str(&capitalized).expect("Invalid given element.");
        Some((element, aromatic, chirality_type))
    }
}

impl Display for SmilesNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.element.as_ref())
    }
}
