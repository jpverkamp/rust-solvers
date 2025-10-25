use point::Point;

use crate::model::color::Color;

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum Tile {
    #[default]
    Water,
    Floor,
    CrackedFloor,
    Nest {
        color: Color,
        dusty: bool,
    },
    Wall,
    Dust,
    Teleport(Point),
}

impl From<char> for Tile {
    fn from(value: char) -> Self {
        match value {
            '~' => Tile::Water,
            '.' => Tile::Floor,
            'x' => Tile::CrackedFloor,
            '#' => Tile::Wall,
            '*' => Tile::Dust,
            _ => unimplemented!("unknown tile {value}"),
        }
    }
}

impl TryFrom<&str> for Tile {
    type Error = ();

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "nest" => Ok(Tile::Nest {
                color: Color::default(),
                dusty: false,
            }),
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
            Tile::Nest { .. } => 'o', // TODO: Support this somehow?
            Tile::Wall => '#',
            Tile::Dust => '*',
            Tile::Teleport(_) => '§',
        }
    }
}
