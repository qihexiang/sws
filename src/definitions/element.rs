use strum_macros::{AsRefStr, EnumString};

#[derive(EnumString, AsRefStr, Debug, PartialEq)]
pub enum Element {
    H = 1,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F,
    Ne,
    Na,
    Mg,
    Al,
    Si,
    P,
    S,
    Cl,
    Ar,
    K,
    Ca,
    Sc,
    Ti,
    V,
    Cr,
    Mn,
    Fe,
    Co,
    Ni,
    Cu,
    Zn,
    Ga,
    Ge,
    As,
    Se,
    Br,
    Kr,
    Rb,
    Sr,
    Y,
    Zr,
    Nb,
    Mo,
    Tc,
    Ru,
    Rh,
    Rd,
    Ag,
    Cd,
    In,
    Sn,
    Sb,
    Te,
    I,
    Xe,
    Cs,
    Ba,
    La,
    Ce,
    Pr,
    Nd,
    Pm,
    Sm,
    Eu,
    Gd,
    Tb,
    Dy,
    Ho,
    Er,
    Tm,
    Yb,
    Lu,
    Hf,
    Ta,
    W,
    Re,
    Os,
    Ir,
    Pt,
    Au,
    Hg,
    Tl,
    Pb,
    Bi,
    Po,
    At,
    Rn,
}

impl Element {
    pub fn default_hydrogen(&self) -> usize {
        match self {
            Self::F | Self::Cl | Self::Br | Self::I => 1,
            Self::O | Self::S => 2,
            Self::B | Self::N | Self::P => 3,
            Self::C | Self::Si => 4,
            _ => 0,
        }
    }

    pub fn is_organic_subset(&self) -> bool {
        match self {
            Self::B
            | Self::C
            | Self::N
            | Self::O
            | Self::F
            | Self::P
            | Self::S
            | Self::Cl
            | Self::Br 
            | Self::I => true,
            _ => false,
        }
    }
}
