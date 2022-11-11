#[derive(Debug, PartialEq)]
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

    pub fn as_str(&self) -> &str {
        match self {
            Self::Clockwise => "@@",
            Self::Counter => "@"
        }
    }
}