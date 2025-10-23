use crate::model::color::Color;

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum Tile {
    #[default]
    Water,
    Floor,
    CrackedFloor,
    Nest(Color),
}

impl From<char> for Tile {
    fn from(value: char) -> Self {
        match value {
            '~' => Tile::Water,
            '.' => Tile::Floor,
            'x' => Tile::CrackedFloor,
            _ => unimplemented!("unknown tile {value}"),
        }
    }
}

impl TryFrom<&str> for Tile {
    type Error = ();

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "nest" => Ok(Tile::Nest(Color::default())),
            _ => Err(()),
        }
    }
}

impl From<Tile> for char {
    fn from(val: Tile) -> Self {
        match val {
            Tile::Water => '~',
            Tile::Floor => '.',
            Tile::CrackedFloor => 'x',
            Tile::Nest(_color) => 'o', // TODO: Support this somehow?
        }
    }
}
