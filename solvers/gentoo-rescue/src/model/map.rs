use direction::Direction;
use point::Point;

use crate::model::{
    color::Color,
    critter::{Critter, CritterKind},
    tile::Tile,
    wall::WallKind,
};

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub(crate) struct Map {
    pub(crate) width: usize,
    pub(crate) height: usize,

    tiles: Vec<Tile>,

    h_walls: Vec<WallKind>,
    v_walls: Vec<WallKind>,

    pub(crate) critters: Vec<Critter>,
    pub(crate) active_critter: usize,
}

impl Map {
    pub(crate) fn tile_at(&self, p: Point) -> Tile {
        if p.x < 0 || p.x >= (self.width as isize) || p.y < 0 || p.y >= (self.height as isize) {
            return Tile::Water;
        }

        let index = (p.y * (self.width as isize) + p.x) as usize;
        self.tiles[index]
    }

    pub(crate) fn break_floor(&mut self, p: Point) {
        assert!(
            p.x >= 0 || p.x < (self.width as isize) || p.y >= 0 || p.y < (self.height as isize),
            "Tried to break a floor out of bounds at {p:?}"
        );

        let index = (p.y * (self.width as isize) + p.x) as usize;
        assert_eq!(
            self.tiles[index],
            Tile::CrackedFloor,
            "Tried to break a non-cracked floor at {p:?}"
        );
        self.tiles[index] = Tile::Water;
    }

    pub(crate) fn wall_index(&self, p: Point, d: Direction) -> Option<(bool, usize)> {
        if p.x < 0 || p.x >= (self.width as isize) || p.y < 0 || p.y >= (self.height as isize) {
            return None;
        }

        let x = p.x as usize;
        let y = p.y as usize;

        match d {
            Direction::Up => Some((true, x + y * self.width)),
            Direction::Down => Some((true, x + (y + 1) * self.width)),
            Direction::Left => Some((false, x + y * (self.width + 1))),
            Direction::Right => Some((false, (x + 1) + y * (self.width + 1))),
        }
    }

    pub(crate) fn wall_at(&self, p: Point, d: Direction) -> WallKind {
        match self.wall_index(p, d) {
            Some((true, index)) => self.h_walls[index],
            Some((false, index)) => self.v_walls[index],
            None => WallKind::Empty,
        }
    }

    pub(crate) fn break_wall(&mut self, p: Point, d: Direction) {
        if let Some(wall) = match self.wall_index(p, d) {
            Some((true, index)) => self.h_walls.get_mut(index),
            Some((false, index)) => self.v_walls.get_mut(index),
            None => None,
        } {
            assert_eq!(
                *wall,
                WallKind::Cracked,
                "Tried to break a n non-cracked wall at {p:?} {d:?}"
            );
            *wall = WallKind::Empty
        }
    }
}

impl From<&str> for Map {
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

        Map {
            width,
            height,
            tiles,
            h_walls,
            v_walls,
            critters,
            active_critter: 0,
        }
    }
}
