use crate::model::color::Color;

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum WallKind {
    #[default]
    Empty,
    Solid,
    Cracked,
    Color(Color),
}

impl From<char> for WallKind {
    fn from(value: char) -> Self {
        match value {
            '.' => WallKind::Empty,
            '|' | '-' => WallKind::Solid,
            ':' | '~' => WallKind::Cracked,
            _ => {
                if let Ok(c) = Color::try_from(value) {
                    WallKind::Color(c)
                } else {
                    unimplemented!("unknown wall kind {value}")
                }
            }
        }
    }
}

impl WallKind {
    pub(crate) fn as_vertical_char(self) -> char {
        match self {
            WallKind::Empty => ' ',
            WallKind::Solid => '|',
            WallKind::Cracked => ':',
            WallKind::Color(_) => '⁞', // TODO: Print colors?
        }
    }

    pub(crate) fn as_horizontal_char(self) -> char {
        match self {
            WallKind::Empty => ' ',
            WallKind::Solid => '-',
            WallKind::Cracked => '╌',
            WallKind::Color(_) => '┈', // TODO: Print colors?
        }
    }
}
