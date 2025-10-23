#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum WallKind {
    #[default]
    Empty,
    Solid,
    Cracked,
}

impl From<char> for WallKind {
    fn from(value: char) -> Self {
        match value {
            '.' => WallKind::Empty,
            '|' | '-' => WallKind::Solid,
            ':' | '~' => WallKind::Cracked,
            _ => unimplemented!("unknown wall kind {value}"),
        }
    }
}

impl WallKind {
    pub(crate) fn as_vertical_char(self) -> char {
        match self {
            WallKind::Empty => ' ',
            WallKind::Solid => '|',
            WallKind::Cracked => ':',
        }
    }

    pub(crate) fn as_horizontal_char(self) -> char {
        match self {
            WallKind::Empty => ' ',
            WallKind::Solid => '-',
            WallKind::Cracked => '╌',
        }
    }
}
