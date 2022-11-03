pub mod definitions;
mod tokenizer;

use definitions::molecule::Molecule;

pub fn create_mole(smiles: &str) -> Molecule {
    Molecule::from_smiles(smiles).unwrap()
}
