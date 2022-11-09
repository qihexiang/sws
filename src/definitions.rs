mod element;
mod bond;
mod chirality;
mod atom;
pub mod selectors;
pub mod workspace;
pub mod smiles;

#[test]
fn generate_node() {
    use workspace::Workspace;
    use element::Element;
    let smiles_strs = vec![
        // "c1cc[13c]cc1",
        // "c1ccccc1C@@(N)(P)S",
        // "c1ccccc1/C=C/(N)P",
        // "c1ccccc1/C=[C{selected}]/(N)P",
        // "N/C=[C{selected}]\\(N(C)[O-])P",
        "c1c([OH2+])cccc1[P{selected}]6(c2ccccc2)Cc3c(cccc3)C[PH2{selected}](c4ccccc4)(c5ccccc5)[Fe+2]6",
        // "C=1CCC1",
        // "C1%12C3C4C1C5C%12C3C45",
        // "[I-].[Na+].C=CCBr",
        // "[Fe++]",
        // "N[C@@](F)(C)C(=O)O",
        // "O1CCCC[C@@H]1C",
        // "[CH2:1]=[CH:2][CH2:1][CH2:3][C:4](C)[CH2:3]"
    ];

    let mut mole = Workspace::new();
    for smiles in smiles_strs.iter() {
        mole.add_smiles(smiles).unwrap();
        println!("{}", smiles);
    }

    println!("{}", mole.dot_representation());
    let p_atoms = mole.filter_nodes(&|atom| atom.element == Element::P);
    mole.reset_root(p_atoms[0]);
    assert_eq!(p_atoms.get(0).copied(), mole.find_root_of(p_atoms[1]));
    println!("{}", mole.dot_representation());
}
