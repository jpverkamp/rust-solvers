use std::rc::Rc;

use direction::Direction;
use point::Point;

use crate::model::{
    color::Color,
    critter::{Critter, CritterKind},
    thing::{Thing, ThingKind},
    tile::Tile,
    toggle::ToggleRule,
    wall::WallKind,
};

#[derive(Debug, Clone)]
pub(crate) struct Map {
    // The size of the map
    pub(crate) width: usize,
    pub(crate) height: usize,

    // The tiles that make up the floors of the map, solid, cracked, water
    // May change (cracked -> water)
    tiles: Vec<Tile>,

    // The horizontal and vertical walls of the level
    // May change (cracked -> water)
    h_walls: Vec<WallKind>,
    v_walls: Vec<WallKind>,

    // The critters moving around the level and the one we're currently moving
    pub(crate) critters: Vec<Critter>,

    // Anything a critter could pick up and carry
    pub(crate) things: Vec<Thing>,

    // Toggles that change the map when walked over
    pub(crate) toggle_rules: Rc<Vec<ToggleRule>>,

    // Sublevels that we should try to enter
    pub(crate) sublevels: Rc<Vec<(Point, String)>>,

    // === State variables while solving ===

    // Current state of teleporters
    // Used to detect infinite loops and avoid double teleports
    pub(crate) used_teleports: Vec<(Direction, Point)>,
    pub(crate) teleport_cooldown: bool,
}

// Only check equality for things that will actually change and aren't state variables
impl PartialEq for Map {
    fn eq(&self, other: &Self) -> bool {
        self.width == other.width
            && self.height == other.height
            && self.tiles == other.tiles
            && self.h_walls == other.h_walls
            && self.v_walls == other.v_walls
            && self.critters == other.critters
            && self.things == other.things
    }
}

impl Eq for Map {}

// Same for hash: Ignore state variables and things that don't change
impl std::hash::Hash for Map {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.width.hash(state);
        self.height.hash(state);
        self.tiles.hash(state);
        self.h_walls.hash(state);
        self.v_walls.hash(state);
        self.critters.hash(state);
        self.things.hash(state);
    }
}

impl Map {
    // Return what tile is at a given location
    pub(crate) fn tile_at(&self, p: Point) -> Tile {
        if p.x < 0 || p.x >= (self.width as isize) || p.y < 0 || p.y >= (self.height as isize) {
            return Tile::Water;
        }

        let index = (p.y * (self.width as isize) + p.x) as usize;
        self.tiles[index]
    }

    // Break the floor at a given location
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

    // Helper to calculate the index and map used for a given wall
    fn wall_index(&self, p: Point, d: Direction) -> Option<(bool, usize)> {
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

    // Return what wall is at a given location in a given direction
    pub(crate) fn wall_at(&self, p: Point, d: Direction) -> WallKind {
        match self.wall_index(p, d) {
            Some((true, index)) => self.h_walls[index],
            Some((false, index)) => self.v_walls[index],
            None => WallKind::Empty,
        }
    }

    // Break the wall in a given location/direction
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

    // Toggle a floor
    pub(crate) fn toggle_floor(&mut self, p: Point) {
        assert!(
            p.x >= 0 || p.x < (self.width as isize) || p.y >= 0 || p.y < (self.height as isize),
            "Tried to toggle a floor out of bounds at {p:?}"
        );

        let index = (p.y * (self.width as isize) + p.x) as usize;
        self.tiles[index] = match self.tiles[index] {
            Tile::Floor => Tile::Water,
            Tile::Water => Tile::Floor,
            other => {
                unimplemented!("Can only toggle floor/water tiles, but tile at {p:?} is {other:?}")
            }
        }
    }

    // Toggle a wall
    pub(crate) fn toggle_wall(&mut self, p: Point, d: Direction) {
        if let Some(wall) = match self.wall_index(p, d) {
            Some((true, index)) => self.h_walls.get_mut(index),
            Some((false, index)) => self.v_walls.get_mut(index),
            None => None,
        } {
            *wall = match *wall {
                WallKind::Empty => WallKind::Solid,
                WallKind::Solid => WallKind::Empty,
                other => unimplemented!(
                    "Can only toggle empty/solid walls, but wall at {p:?} {d:?} is {other:?}"
                ),
            }
        }
    }
}

impl From<&str> for Map {
    fn from(input: &str) -> Self {
        let mut lines = input.lines().peekable();

        // Skip any leading # comment or empty lines
        while let Some(&line) = lines.peek() {
            if line.is_empty() || line.starts_with('#') {
                lines.next();
            } else {
                break;
            }
        }

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

        // Then the critters, things, and setting nests in tiles
        let mut critters = vec![];
        let mut things = vec![];
        let mut toggle_rules = vec![];
        let mut sublevels = vec![];

        for line in lines {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let parts: Vec<_> = line.split_ascii_whitespace().collect();
            assert!(
                parts.len() >= 3,
                "Need at least a row, col, and name in {line}"
            );

            // Toggles always start with "toggle"
            if parts[0] == "toggle" {
                toggle_rules.push(ToggleRule::from(line));
                continue;
            }

            // TODO: Why did I make these 1 based...
            let row = parts[0]
                .parse::<usize>()
                .expect("Critter row should be a number")
                - 1;
            let col = parts[1]
                .parse::<usize>()
                .expect("Critter col should be a number")
                - 1;

            if parts.len() == 4
                && let Ok(kind) = CritterKind::try_from(parts[3])
            {
                // Try to load critters: 1 1 red penguin
                let color = Color::from(parts[2]);
                critters.push(Critter::new(kind, color, (col, row).into()));
            } else if parts.len() == 4
                && let Ok(mut tile) = Tile::try_from(parts[3])
            {
                // Try to load a nest: 1 1 red nest
                let color = Color::from(parts[2]);
                let index = row * width + col;

                // TOOD: Can this be done better?
                if let Tile::Nest {
                    color: nest_color,
                    dusty,
                } = &mut tile
                {
                    *nest_color = color;
                    if tiles[index] == Tile::Dust {
                        *dusty = true;
                    }
                }

                tiles[index] = tile;
            } else if parts.len() == 3
                && let Ok(kind) = ThingKind::try_from(parts[2])
            {
                // Things are ... things(?) the critters can carry around
                things.push(Thing {
                    kind,
                    location: Point::from((col, row)),
                });
            } else if parts[2] == "teleport" {
                let dest_row = parts[3]
                    .parse::<usize>()
                    .expect("Teleport destination row should be a number")
                    - 1;
                let dest_col = parts[4]
                    .parse::<usize>()
                    .expect("Teleport destination col should be a number")
                    - 1;

                let index = row * width + col;
                tiles[index] = Tile::Teleport((dest_col, dest_row).into());
            } else if parts[2] == "level" {
                let name = parts[3];
                sublevels.push((Point::from((col, row)), name.to_string()));
            } else {
                // If we made it this far, it's a bad object (probably?)
                panic!("Malformed object: {line}");
            }
        }

        // Any things that start on a tile with a critter get picked up immediately
        for critter in &mut critters {
            if let Some(index) = things.iter().position(|t| t.location == critter.location()) {
                critter.pick_up(things[index].kind);
                things.remove(index);
            }
        }

        Map {
            width,
            height,
            tiles,
            h_walls,
            v_walls,
            critters,
            things,
            toggle_rules: Rc::new(toggle_rules),
            sublevels: Rc::new(sublevels),

            used_teleports: vec![],
            teleport_cooldown: false,
        }
    }
}
