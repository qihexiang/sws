#[derive(Debug)]
pub enum ChiralityType {
    Clockwise,
    Counter,
}

impl ChiralityType {
    pub fn new(s: &str) -> Option<Self> {
        match s {
            "@" => Some(Self::Counter),
            "@@" => Some(Self::Clockwise),
            _ => None,
        }
    }
}