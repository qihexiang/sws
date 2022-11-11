mod element;
mod bond;
mod chirality;
mod atom;
pub mod selectors;
pub mod workspace;
pub mod decode;
pub mod encode;

#[test]
fn generate_node() {
    use workspace::Workspace;
    use element::Element;
    let smiles_strs = vec![
        "c1cc[13c]cc1",
        "c1ccccc1C@@(N)(P)S",
        "c1ccccc1/C=C/(N)P",
        "c1ccccc1/C=[C{selected}]/(N)\\P",
        "N/C=[C{selected}]\\(N(C)[O-])/P",
        "c1c([OH2+])cccc1[P]6(c2ccccc2)Cc3c(cccc3)C[PH2](c4ccccc4)(c5ccccc5)[Fe+2]6",
        "C=1CCC1",
        "C1%12C3C4C1C5C%12C3C45",
        "[I-].[Na+].C=CCBr",
        "[Fe++]",
        "N[C@@](F)(C)C(=O)O",
        "O1CCCC[C@@H]1C",
        "[CH2:1]=[CH:2][CH2:1][CH2:3][C:4](C)[CH2:3]"
    ];

    let mut ws = Workspace::new();
    for smiles in smiles_strs.iter() {
        ws.add_smiles(smiles).unwrap();
        // println!("{}", smiles);
    }

    // let atoms = ws.filter_nodes(&|atom| atom.element == Element::C);
    // ws.reset_root(atoms[4]);

    for structure_root in ws.find_strcutres_roots() {
        ws.add_hydrogen_to_structure(structure_root);
        let smiles = ws.to_smiles(structure_root);
        if let Some(smiles) = smiles {
            println!("{}", smiles);
        } else {
            println!("Failed to generate SMILES")
        }
    }
}
