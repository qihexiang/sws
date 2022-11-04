use lazy_static::lazy_static;
use regex::Regex;

lazy_static! {
    pub static ref REPLACE_SELECTOR: Regex = Regex::new(r"Replacer\((?P<keytype>[\-=#$])((\s)*,(\s)*(?P<amount>[1-9][0-9]*))?\)").unwrap();
}

#[test]
fn match_replace_selector() {
    let result = vec![
        ("=", 1), ("-", 2), ("=", 1), ("-",1)
    ];
    vec![
        "Replacer(=)", "Replacer(-, 2)", "Replacer(=,1)", "Replacer(-)"
    ].iter()
    .map(|s| REPLACE_SELECTOR.captures(s))
    .map(|cap| {
        cap.and_then(|cap| Some((cap.name("keytype").and_then(|m| Some(m.as_str())), cap.name("amount").and_then(|m| Some(m.as_str().parse::<u8>().unwrap_or(1))).unwrap_or(1))))
    })
    .enumerate()
    .for_each(|(idx, res)| {
        if let Some((Some(keytype), amount)) = res {
            let (target_type, target_amount) = result[idx];
            assert_eq!(keytype, target_type);
            assert_eq!(amount, target_amount);
        } else {
            debug_assert!(false, "{:?}", res)
        }
    })
}