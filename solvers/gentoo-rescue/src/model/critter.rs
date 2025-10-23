use point::Point;

use crate::model::{color::Color, thing::ThingKind};

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
    pub(crate) carrying: Option<ThingKind>,
}

impl std::fmt::Display for Critter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!(
            "{row} {col} {color:?} {kind:?}{carrying}",
            row = self.location.y + 1,
            col = self.location.x + 1,
            color = self.color,
            kind = self.kind,
            carrying = if let Some(thing) = self.carrying {
                format!(" w/{thing:?}")
            } else {
                String::new()
            }
        ))
    }
}
