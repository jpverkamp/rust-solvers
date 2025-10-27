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
    kind: CritterKind,
    color: Color,
    location: Point,
    carrying: Option<ThingKind>,
    escaped: bool,
}

impl Critter {
    pub fn new(kind: CritterKind, color: Color, location: Point) -> Self {
        Critter {
            kind,
            color,
            location,
            carrying: None,
            escaped: false,
        }
    }

    pub fn kind(&self) -> CritterKind {
        self.kind
    }

    pub fn color(&self) -> Color {
        self.color
    }

    pub fn location(&self) -> Point {
        self.location
    }

    pub fn move_to(&mut self, p: Point) {
        self.location = p;
    }

    pub fn carrying(&self) -> Option<ThingKind> {
        self.carrying
    }

    pub fn pick_up(&mut self, t: ThingKind) {
        self.carrying = Some(t);
    }

    pub fn escaped(&self) -> bool {
        self.escaped
    }

    pub fn escape(&mut self) {
        self.escaped = true;
        self.location = Point { x: -10, y: -10 };
    }
}

impl std::fmt::Display for Critter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.pad(&format!(
            "{row} {col} {color:?} {kind:?}{carrying}{escaped}",
            row = self.location.y + 1,
            col = self.location.x + 1,
            color = self.color,
            kind = self.kind,
            carrying = if let Some(thing) = self.carrying {
                format!(" w/{thing:?}")
            } else {
                String::new()
            },
            escaped = if self.escaped {
                " {ESC}".to_string()
            } else {
                String::new()
            }
        ))
    }
}
