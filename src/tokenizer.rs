use lazy_static::lazy_static;
use regex::Regex;

lazy_static! {
    pub static ref SSMILES_RE: Regex = Regex::new(r"\[([1-9][0-9]*)?((br?|cl?|n|o|p|s|f|i)|([A-Z][a-z]?))(@{0,2})(H([1-9][0-9]*)?)?(((\+|\-)([1-9][0-9]*))|(\+*)|(\-*))(:([0-9]+))?(\{.+?\})?\]|(((br?|cl?|n|o|p|s|f|i)|(Br?|Cl?|N|O|P|S|F|I))(@{0,2}))|\(|\)|\.|((\-|=|#|$|:)?([1-9]|(%[1-9][0-9]+)))|\-|=|#|$|:|/|\\").unwrap();
    pub static ref BOND_RE: Regex = Regex::new(r"^(\-|=|#|$|:|/|\\)$").unwrap();
    pub static ref NOTHING_RE: Regex = Regex::new(r"^(\.)$").unwrap();
    pub static ref RING_BOND_RE: Regex = Regex::new(r"^((?P<bond_type>\-|=|#|$|:|\|//)?(?P<ring_id>([1-9])|(%[1-9][0-9]+)))$").unwrap();
    pub static ref ORGANIC_SUBSET_RE: Regex = Regex::new("^((?P<element>(br?|cl?|n|o|p|s|f|i)|(Br?|Cl?|N|O|P|S|F|I))(?P<chirality>@{0,2}))$").unwrap();
    pub static ref AROMATIC_ORGANIC_RE: Regex = Regex::new("^(br?|cl?|n|o|p|s|f|i)$").unwrap();
    pub static ref STANDARD_NODE_RE: Regex = Regex::new(r"^(\[(?P<isotope>[1-9][0-9]*)?(?P<element>(br?|cl?|n|o|p|s|f|i)|([A-Z][a-z]?))(?P<chirality>@{0,2})(?P<explicit_hydrogen>H(?P<explicit_hydrogen_num>[1-9][0-9]*)?)?(?P<charge>((?P<charge_num>(\+|\-)([1-9][0-9]*))|(\+*)|(\-*)))(:(?P<react_id>[0-9]+))?(?P<selector>\{.+?\})?\])$").unwrap();
    pub static ref BRANCH_RE: Regex = Regex::new(r"^(\(|\))$").unwrap();
    pub static ref NAGETIVE_RE: Regex = Regex::new(r"^(\-+)$").unwrap();
    pub static ref POSITIVE_RE: Regex = Regex::new(r"^(\++)$").unwrap();
}

pub fn smiles_tokenize(smiles: &str) -> Vec<&str> {
    SSMILES_RE.find_iter(smiles).map(|r| r.as_str()).collect::<Vec::<&str>>()
}
