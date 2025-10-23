use direction::Direction;
use point::Point;

use crate::local::Local;

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum Color {
    #[default]
    Red,
    Yellow,
    Green,
    Blue,
}

impl From<&str> for Color {
    fn from(value: &str) -> Self {
        match value {
            "red" => Color::Red,
            "yellow" => Color::Yellow,
            "green" => Color::Green,
            "blue" => Color::Blue,
            _ => unimplemented!("Unknown color {value}"),
        }
    }
}

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
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

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
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

#[derive(Debug, Clone, Default)]
pub struct Global {
    pub width: usize,
    pub height: usize,
    tiles: Vec<Tile>,
    h_walls: Vec<WallKind>,
    v_walls: Vec<WallKind>,
    initial_critters: Vec<Critter>,
}

impl Global {
    pub(crate) fn tile_at(&self, p: Point) -> Tile {
        if p.x < 0 || p.x >= (self.width as isize) || p.y < 0 || p.y >= (self.height as isize) {
            return Tile::Water;
        }

        let index = (p.y * (self.width as isize) + p.x) as usize;
        self.tiles[index]
    }

    pub(crate) fn wall_at(&self, p: Point, d: Direction) -> WallKind {
        if p.x < 0 || p.x >= (self.width as isize) || p.y < 0 || p.y >= (self.height as isize) {
            return WallKind::Empty;
        }

        let x = p.x as usize;
        let y = p.y as usize;

        match d {
            Direction::Up => self.h_walls[x + y * self.width],
            Direction::Down => self.h_walls[x + (y + 1) * self.width],
            Direction::Left => self.v_walls[x + y * (self.width + 1)],
            Direction::Right => self.v_walls[(x + 1) + y * (self.width + 1)],
        }
    }

    pub(crate) fn make_local(&self) -> Local {
        Local {
            critters: self.initial_critters.clone(),
        }
    }

    pub(crate) fn critter(&self, i: usize) -> Critter {
        self.initial_critters[i]
    }
}

impl From<&str> for Global {
    fn from(input: &str) -> Self {
        let mut lines = input.lines();

        // First read the tiles
        let mut width = None;
        let mut tiles = vec![];

        while let Some(line) = lines.next()
            && !line.is_empty()
        {
            if width.is_none() {
                width = Some(line.len());
            }
            assert_eq!(
                width,
                Some(line.len()),
                "Unexpected line width while reading tiles: {line}"
            );

            for c in line.chars() {
                tiles.push(Tile::from(c));
            }
        }

        let width = width.unwrap();
        let height = tiles.len() / width;

        // Then the vertical walls
        let mut v_walls = vec![];
        let mut row = 0;

        while let Some(line) = lines.next()
            && !line.is_empty()
        {
            assert_eq!(
                width + 1,
                line.len(),
                "Unexpected line width while reading vertical walls: {line}"
            );

            for c in line.chars() {
                v_walls.push(c.into());
            }
            row += 1;
        }
        assert_eq!(row, height, "Not enough rows when reading vertical walls");

        // Then the horizontal walls
        let mut h_walls = vec![];
        let mut row = 0;

        while let Some(line) = lines.next()
            && !line.is_empty()
        {
            assert_eq!(
                width,
                line.len(),
                "Unexpected line width while reading horizontal walls: {line}"
            );

            for c in line.chars() {
                h_walls.push(c.into());
            }

            row += 1;
        }
        assert_eq!(
            row,
            height + 1,
            "Not enough rows when reading horizontal walls"
        );

        // Then the critters (plus setting nests in tiles)
        let mut critters = vec![];
        while let Some(line) = lines.next()
            && !line.is_empty()
        {
            let parts: Vec<_> = line.split_ascii_whitespace().collect();
            assert_eq!(parts.len(), 4, "Malformed critter at {line}");

            // TODO: Why did I make these 1 based...
            let row = parts[0]
                .parse::<usize>()
                .expect("Critter row should be a number")
                - 1;
            let col = parts[1]
                .parse::<usize>()
                .expect("Critter col should be a number")
                - 1;
            let color = Color::from(parts[2]);

            if let Ok(kind) = CritterKind::try_from(parts[3]) {
                critters.push(Critter {
                    kind,
                    color,
                    location: Point::from((col, row)),
                })
            } else if let Ok(mut tile) = Tile::try_from(parts[3]) {
                // TOOD: Can this be done better?
                if let Tile::Nest(nest_color) = &mut tile {
                    *nest_color = color;
                }

                let index = row * width + col;
                tiles[index] = tile;
            }
        }

        Global {
            width,
            height,
            tiles,
            h_walls,
            v_walls,
            initial_critters: critters,
        }
    }
}
