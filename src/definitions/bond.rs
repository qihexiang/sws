use std::fmt::Display;

use crate::tokenizer::BOND_RE;

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum BondType {
    Single,
    UpSingle,
    DownSingle,
    Double,
    Triple,
    Quad,
    Aromatic,
    NoBond,
}

impl BondType {
    pub fn new(bond_str: &str) -> Option<Self>{
        if BOND_RE.is_match(bond_str) {
            Some(match bond_str {
                "-" => BondType::Single,
                "/" => BondType::UpSingle,
                "\\" => BondType::DownSingle,
                "=" => BondType::Double,
                "#" => BondType::Triple,
                "$" => BondType::Quad,
                ":" => BondType::Aromatic,
                _ => panic!("This shall never happend"),
            })
        } else {
            None
        }
    }
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub struct Bond {
    bond_type: BondType,
    ring: bool
}

impl Bond {
    pub fn new(bond_type: BondType, ring: bool) -> Self {
        Bond {
            bond_type, ring
        }
    }

    pub fn is_no_bond(&self) -> bool {
        if let BondType::NoBond = self.bond_type {
            true
        } else {
            false
        }
    }

    pub fn is_ring_bond(&self) -> bool {
        self.ring
    } 

    pub fn reverse(self) -> Self {
        match self.bond_type {
            BondType::UpSingle => Bond::new(BondType::DownSingle, self.ring),
            BondType::DownSingle => Bond::new(BondType::UpSingle, self.ring),
            _ => self
        }
    } 
}

impl Display for Bond {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}",
            match self.bond_type {
                BondType::Double => "=",
                BondType::Triple => "#",
                BondType::Quad => "$",
                BondType::Aromatic => ":",
                BondType::Single => "-",
                BondType::UpSingle => "\\",
                BondType::DownSingle => "/",
                BondType::NoBond => "."
            },
            if self.ring {
                " ring key"
            } else {
                ""
            }
        )
    }
}
