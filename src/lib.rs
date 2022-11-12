pub mod definitions;
mod tokenizer;
pub mod workspace;
use definitions::bond::BondType;
use petgraph::graph::NodeIndex;
use rand::prelude::*;
use workspace::Workspace;

/// 连接规则：Out-In，R-RIn
pub fn random_generate_structures(
    start: &str,
    end: &str,
    duals: Vec<&str>,
    singles: Vec<&str>,
    duals_amount: usize,
) -> Option<String> {
    let mut ws = Workspace::new();
    let start_point = ws.add_smiles(start).unwrap();
    let find_with_selector = |ws: &Workspace, root: NodeIndex, target_selector: &str| {
        ws.find_node_in_structure(root, |atom| {
            if let Some(selector) = &atom.selector {
                let selectors = selector.split(";").collect::<Vec<&str>>();
                selectors.contains(&target_selector)
            } else {
                false
            }
        })
    };

    let remove_selector = |ws: &mut Workspace, target: NodeIndex, target_selector: &str| {
        ws.get_atom_mut(target).and_then(|atom| {
            if let Some(selectors) = &atom.selector {
                let selectors = selectors.split(";").collect::<Vec<_>>();
                let outgoing_index = selectors
                    .iter()
                    .position(|selector| selector == &target_selector)
                    .unwrap();
                let before = &selectors[0..outgoing_index].join(";");
                let after = &selectors[outgoing_index + 1..].join(";");
                let mut updated_selectors = String::new();
                updated_selectors.push_str(before);
                updated_selectors.push_str(";");
                updated_selectors.push_str(after);
                if updated_selectors == ";" {
                    atom.selector = None
                } else {
                    atom.selector = Some(updated_selectors);
                }
                Some(&atom.selector)
            } else {
                None
            }
        });
    };

    for _ in 0..duals_amount {
        let dual_root = ws.add_smiles(*random_take_one(&duals)).unwrap();
        let outgoing = find_with_selector(&ws, start_point, "Out")
            .expect("No Out selector found, but still have something to add.");
        let incoming = find_with_selector(&ws, dual_root, "In")
            .expect("Give fragment doesn't have In selector");
        ws.reset_root(incoming);
        ws.connect(outgoing, incoming, BondType::Single);
        remove_selector(&mut ws, outgoing, "Out");
        remove_selector(&mut ws, incoming, "In");
    }
    let outgoing = find_with_selector(&ws, start_point, "Out")
        .expect("No In selector found, but still have something to add.");
    let end_root = ws.add_smiles(end).unwrap();
    let end_incoming =
        find_with_selector(&ws, end_root, "In").expect("Give fragment doesn't have Out selector");
    ws.reset_root(end_incoming);
    ws.connect(outgoing, end_incoming, BondType::Single);
    remove_selector(&mut ws, outgoing, "Out");
    remove_selector(&mut ws, end_incoming, "In");
    while let Some(r_outgoing) = find_with_selector(&ws, start_point, "R") {
        let r_root = ws.add_smiles(*random_take_one(&singles)).unwrap();
        let r_incoming = find_with_selector(&ws, r_root, "RIn")
            .expect("Single replacer must have at least one RIn selector");
        ws.connect(r_outgoing, r_incoming, BondType::Single);
        remove_selector(&mut ws, r_outgoing, "R");
        remove_selector(&mut ws, r_incoming, "RIn");
    }
    ws.to_smiles(start_point)
}

fn random_take_one<'a, E>(collection: &'a Vec<E>) -> &'a E {
    let size = collection.len();
    let mut rng = thread_rng();
    let rand_value = rng.gen_range(0..size);
    &collection[rand_value]
}

#[test]
fn generate_structures() {
    for _ in 0..100 {
        let result = random_generate_structures(
            "[P{R;R;Out}]",
            "[P{R;R;In}]",
            vec![
                "[c{In}]1cccc[c{Out}]1",
                "[c{In}]1ccc[c{Out}]c1",
                "[c{In}]1cc[c{Out}]cc1",
                "[C{R;In;Out}]",
            ],
            vec![
                "[H{RIn}]",
                "[C{RIn}]",
                "[OH{RIn}]",
                "[NH2{RIn}]",
                "[c{RIn}]1ccccc1",
            ],
            2,
        )
        .unwrap();
        println!("{}", result);
    }
}
