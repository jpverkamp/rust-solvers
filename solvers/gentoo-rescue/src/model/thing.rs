use point::Point;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum ThingKind {
    Spring,
    Hammer,
}

impl TryFrom<&str> for ThingKind {
    type Error = ();

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "spring" => Ok(ThingKind::Spring),
            "hammer" => Ok(ThingKind::Hammer),
            _ => Err(()),
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Thing {
    pub(crate) kind: ThingKind,
    pub(crate) location: Point,
}
