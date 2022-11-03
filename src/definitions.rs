mod element;
mod bond;
mod chirality;
mod smiles_node;
pub mod molecule;

#[test]
fn generate_node() {
    use molecule::Molecule;
    let smiles_strs = vec![
        "c1cc[13c]cc1",
        "c1ccccc1C@@(N)(P)S",
        "c1ccccc1/C=C/(N)P",
        "c1ccccc1/C=[C{selected}]/(N)P",
        "N/C=[C{selected}]\\(N(C)[O-])P",
        "c1c([OH2+])cccc1[P{selected}]6(c2ccccc2)Cc3c(cccc3)C[PH2{selected}](c4ccccc4)(c5ccccc5)[Fe+2]6",
        "C=1CCC1",
        "C1%12C3C4C1C5C%12C3C45",
        "[I-].[Na+].C=CCBr>>[Na+].[Br-].C=CCI",
        "[Fe++]",
        "N[C@@](F)(C)C(=O)O",
        "O1CCCC[C@@H]1C",
        "[CH2:1]=[CH:2][CH2:1][CH2:3][C:4](C)[CH2:3]"
    ];

    for smiles in smiles_strs {
        println!("{}", smiles);
        let mole = Molecule::from_smiles(smiles);
        if let Ok(mole) = mole {
            println!("{:?}", mole);
        }
    }
}
