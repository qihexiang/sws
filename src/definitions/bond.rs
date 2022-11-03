use crate::tokenizer::BOND_RE;

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum BondType {
    Single(SingleBondType),
    Double,
    Triple,
    Quad,
    Aromatic,
    NoBond,
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum SingleBondType {
    Normal,
    LeftUp,
    RightUp,
}

impl BondType {
    pub fn simple() -> Self {
        Self::Single(SingleBondType::Normal)
    }
    
    pub fn from_str(s: &str) -> Option<Self> {
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