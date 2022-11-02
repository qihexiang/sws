use lazy_static::lazy_static;
use regex::Regex;

lazy_static! {
    pub static ref SSMILES_RE: Regex = Regex::new(r"\[([1-9][0-9]*)?((br?|cl?|n|o|p|s|f|i)|([A-Z][a-z]?))(@{0,2})(H([1-9][0-9]*)?)?(((\+|\-)([1-9][0-9]*))|(\+*)|(\-*))(\{.+?\})?\]|(((br?|cl?|n|o|p|s|f|i)|(Br?|Cl?|N|O|P|S|F|I))(@{0,2}))|\(|\)|\.|((\-|=|#|$|:)?([1-9]|(%[1-9][0-9]+)))|\-|=|#|$|:|/|\\").unwrap();
    pub static ref BOND_RE: Regex = Regex::new(r"^(\-|=|#|$|:|/|\\)$").unwrap();
    pub static ref NOTHING_RE: Regex = Regex::new(r"^(\.)$").unwrap();
    pub static ref RING_BOND_RE: Regex = Regex::new(r"^((?P<bond_type>\-|=|#|$|:|\|//)?(?P<ring_id>[1-9]|(%[1-9][0-9]+)))$").unwrap();
    pub static ref ORGANIC_SUBSET_RE: Regex = Regex::new("^((?P<element>(br?|cl?|n|o|p|s|f|i)|(Br?|Cl?|N|O|P|S|F|I))(?P<chirality>@{0,2}))$").unwrap();
    pub static ref AROMATIC_ORGANIC_RE: Regex = Regex::new("^(br?|cl?|n|o|p|s|f|i)$").unwrap();
    pub static ref STANDARD_NODE_RE: Regex = Regex::new(r"^(\[(?P<isotope>[1-9][0-9]*)?(?P<element>(br?|cl?|n|o|p|s|f|i)|([A-Z][a-z]?))(?P<chirality>@{0,2})(H(?P<explicit_hydrogen>[1-9][0-9]*)?)?(?P<charge>((?P<charge_num>(\+|\-)([1-9][0-9]*))|(\+*)|(\-*)))(?P<selector>\{.+?\})?\])$").unwrap();
    pub static ref BRANCH_RE: Regex = Regex::new(r"^(\(|\))$").unwrap();
    pub static ref NAGETIVE_RE: Regex = Regex::new(r"^(\-+)$").unwrap();
    pub static ref POSITIVE_RE: Regex = Regex::new(r"^(\++)$").unwrap();
}

pub fn smiles_tokenize(smiles: &str) -> Vec<&str> {
    SSMILES_RE.find_iter(smiles).map(|r| r.as_str()).collect::<Vec::<&str>>()
}

#[test]
fn re_test() {
    let smiles_strs = vec![
        "c1cc[13c]cc1",
        "c1ccccc1C@@(N)(P)S",
        "c1ccccc1/C=C/(N)P",
        "c1ccccc1/C=[C{selected}]/(N)P",
        "N/C=[C{selected}]\\(N(C)[O-])P",
        "c1ccccc1[P{selected}](c2ccccc2)Cc3c(cccc3)C[P{selected}](c4ccccc4)c5ccccc5",
        "C=1CCC1"
    ];
    for smiles in smiles_strs.iter() {
        // for caps in SSMILES_RE.captures_iter(smiles) {
        //     let tokens: Vec<&str> = caps.iter().map(|c| c.unwrap().as_str()).collect();
        //     println!("{:?}", tokens);
        // }
        println!("{} {:?}", smiles, SSMILES_RE.find_iter(smiles).map(|r| r.as_str()).collect::<Vec::<&str>>())
    }
}
