use point::Point;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum ThingKind {
    Spring,
    Hammer,
    Crutch,
}

impl TryFrom<&str> for ThingKind {
    type Error = ();

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "spring" => Ok(ThingKind::Spring),
            "hammer" => Ok(ThingKind::Hammer),
            "crutch" => Ok(ThingKind::Crutch),
            _ => Err(()),
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Thing {
    pub(crate) kind: ThingKind,
    pub(crate) location: Point,
}
