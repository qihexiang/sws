pub mod definitions;
mod tokenizer;

use definitions::{molecule::Molecules, selectors::REPLACE_SELECTOR};

// pub fn create_mole(smiles: &str) -> Molecule {
//     Molecule::add_smiles(smiles).unwrap()
// }

// #[test]
// fn create_mole_example() {
//     let mut mole = create_mole("[P{Replacer(-,2)}]Cc3c(cccc3)CP(c4ccccc4)(c5ccccc5)");
//     mole.add_hydrogens();
//     println!("{}", mole.dot_representation());
//     mole.remove_hydrogens();
//     println!("{}", mole.dot_representation());
// }
