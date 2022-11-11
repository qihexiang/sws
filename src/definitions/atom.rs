use std::{fmt::Display, str::FromStr};

use super::{chirality::ChiralityType, element::Element};
use crate::tokenizer::{AROMATIC_ORGANIC_RE, NAGETIVE_RE, ORGANIC_SUBSET_RE, STANDARD_NODE_RE};

#[derive(Debug)]
pub struct Atom {
    pub element: Element,
    pub isotope: Option<u16>,
    pub charge: isize,
    pub chirality_type: Option<ChiralityType>,
    pub explicit_hydrogen: usize,
    pub selector: Option<String>,
    pub aromatic: bool,
    pub react_id: Option<usize>,
}

impl Atom {
    pub fn to_token(&self) -> String {
        let mut token = String::new();
        if self.element.is_organic_subset()
            && self.isotope == None
            && self.charge == 0
            && self.explicit_hydrogen == 0
            && self.selector == None
            && self.react_id == None
        {
            token.push_str(&self.core_token());
        } else {
            token.push_str("[");
            if let Some(isotope) = self.isotope {
                token.push_str(&isotope.to_string());
            }
            token.push_str(&self.core_token());
            if self.explicit_hydrogen != 0 {
                token.push_str("H");
                if self.explicit_hydrogen > 1 {
                    token.push_str(&self.explicit_hydrogen.to_string());
                }
            }
            token.push_str(&self.charge_token());
            if let Some(react_id) = self.react_id {
                token.push_str(":");
                token.push_str(&react_id.to_string());
            }
            if let Some(selector) = &self.selector {
                token.push_str("{");
                token.push_str(selector);
                token.push_str("}");
            }
            token.push_str("]");
        }
        token
    }

    fn core_token(&self) -> String {
        let mut core_token = String::new();
        let element_token = self.element.as_ref();
        if self.aromatic && ORGANIC_SUBSET_RE.is_match(element_token) {
            core_token.push_str(&element_token.to_lowercase());
        } else {
            core_token.push_str(element_token);
        }
        if let Some(chirality) = &self.chirality_type {
            core_token.push_str(chirality.as_str());
        }
        core_token
    }

    fn charge_token(&self) -> String {
        let mut charge_token = String::new();
        if self.charge < 0 {
            charge_token.push_str("-");
        } else if self.charge > 0 {
            charge_token.push_str("+");
        }
        if self.charge.abs() > 1 {
            charge_token.push_str(&self.charge.abs().to_string())
        }
        charge_token
    }

    pub fn new(token: &str) -> Option<Self> {
        if let Some(captured) = ORGANIC_SUBSET_RE.captures(token) {
            let (element, aromatic, chirality_type) = Self::minimal_node_info(&captured)?;
            Some(Atom {
                element,
                isotope: None,
                charge: 0,
                chirality_type,
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
                .and_then(|_| match captured.name("explicit_hydrogen_num") {
                    Some(num) => Some(num.as_str().parse::<usize>().unwrap_or(1)),
                    None => Some(1),
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
            Some(Atom {
                element,
                isotope,
                charge,
                chirality_type,
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

impl Display for Atom {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.element.as_ref())
    }
}
