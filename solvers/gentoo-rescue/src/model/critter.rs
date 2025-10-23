use point::Point;

use crate::model::color::Color;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum CritterKind {
    Penguin,
    Seal,
}

impl TryFrom<&str> for CritterKind {
    type Error = ();

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "penguin" => Ok(CritterKind::Penguin),
            "seal" => Ok(CritterKind::Seal),
            _ => Err(()),
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct Critter {
    pub(crate) kind: CritterKind,
    pub(crate) color: Color,
    pub(crate) location: Point,
}

impl Critter {
    pub(crate) fn index_char(index: usize) -> char {
        if index < 10 {
            (b'0' + (index as u8)) as char
        } else if index < 10 + 26 {
            (b'a' + ((index - 10) as u8)) as char
        } else if index < 10 + 26 + 26 {
            (b'A' + ((index - 10 - 26) as u8)) as char
        } else {
            unimplemented!("Too many critters!")
        }
    }
}

impl std::fmt::Display for Critter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!(
            "{row} {col} {color:?} {kind:?}",
            row = self.location.y + 1,
            col = self.location.x + 1,
            color = self.color,
            kind = self.kind,
        ))
    }
}
