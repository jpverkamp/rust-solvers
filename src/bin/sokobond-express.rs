use std::{io::{self, BufRead, BufReader}, ops::{Add, AddAssign, Sub}};
use fxhash::{FxHashMap, FxHashSet};
use solver::{Solver, State};
use uuid::Uuid;

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
struct Point { 
    x: i32,
    y: i32 
}

impl Add<Point> for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        Point {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl AddAssign<Point> for Point {
    fn add_assign(&mut self, other: Point) {
        *self = *self + other;
    }
}

impl Sub<Point> for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        Point {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Point {
    fn manhattan_distance(&self, other: &Point) -> i32 {
        (self.x - other.x).abs() + (self.y - other.y).abs()
    }
}


#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
enum Direction {
    Up,
    Down,
    Left,
    Right,
}

impl Into<Direction> for &str {
    fn into(self) -> Direction {
        match self {
            "Up" => Direction::Up,
            "Down" => Direction::Down,
            "Left" => Direction::Left,
            "Right" => Direction::Right,
            _ => panic!("Unknown direction: {}", self),
        }
    }
}

impl From<Direction> for Point {
    fn from(value: Direction) -> Self {
        match value {
            Direction::Up => Point { x: 0, y: -1 },
            Direction::Down => Point { x: 0, y: 1 },
            Direction::Left => Point { x: -1, y: 0 },
            Direction::Right => Point { x: 1, y: 0 },
        }
    }
}

// Global state, the map and final goal

#[derive(Clone, Debug)]
struct Map {
    // Size of the map
    width: u32,
    height: u32,

    // Unwalkable parts of the map 
    unwalkable: FxHashSet<Point>,

    // Map elements
    atoms: Vec<Atom>,
    tracks: Vec<Track>,
}

// Local state 

#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq)]
enum Element {
    Ion,
    Hydrogen,
    Oxygen,
    Nitrogen,
    Carbon,
}

impl Element {
    fn max_electrons(&self) -> u8 {
        match self {
            Element::Ion => 1,
            Element::Hydrogen => 1,
            Element::Oxygen => 2,
            Element::Nitrogen => 3,
            Element::Carbon => 4,
        }
    }
}

impl Into<Element> for &str {
    fn into(self) -> Element {
        match self {
            "Ion" => Element::Ion,
            "Hydrogen" => Element::Hydrogen,
            "Oxygen" => Element::Oxygen,
            "Nitrogen" => Element::Nitrogen,
            "Carbon" => Element::Carbon,
            _ => panic!("Unknown element kind: {}", self),
        }
    }
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
struct Atom {
    kind: Element,
    position: Point,
    max_electrons: u8,
    free_electrons: u8,
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
struct Track {
    from: (Point, Direction),
    to: (Point, Direction),
    via: Vec<Point>,
}

impl Track {
    fn contains(&self, point: &Point) -> bool {
        point == &self.from.0 || point == &self.to.0 || self.via.contains(point)
    }
}

impl From<&str> for Map {
    fn from(value: &str) -> Self {
        let mut line_reader = BufReader::new(value.as_bytes());

        #[derive(Copy, Clone, Debug)]
        enum ReadState {
            Default,
            Map,
            Track,
        }

        let mut state = ReadState::Default;
        let mut line_no = 0;

        let mut map = Map {
            width: 0,
            height: 0,
            unwalkable: FxHashSet::default(),
            atoms: vec![],
            tracks: vec![],
        };

        let mut current_track = None;

        loop {
            line_no += 1;

            let mut line = String::new();
            match line_reader.read_line(&mut line) {
                Ok(0) => break,
                Ok(_) => {},
                Err(e) => panic!("Error reading line: {}", e),
            }
            let line = line.trim();
            
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            match state {
                ReadState::Default => {
                    if line.starts_with("map") {

                        // A map block sets the width and height and can contain 'cutouts'
                        let parts = line.split_ascii_whitespace().collect::<Vec<_>>();
                        assert!(parts.len() == 3, "[{line_no}] Excepted 'map <width> <height>', got {line}");
                        map.width = parts[1].parse::<u32>().unwrap();
                        map.height = parts[2].parse::<u32>().unwrap();
                        state = ReadState::Map;

                    } else if line.starts_with("track") {

                        // A track contains a start point, end point, and a series of points along the way
                        current_track = Some(Track {
                            from: (Point { x: 0, y: 0 }, Direction::Up),
                            to: (Point { x: 0, y: 0 }, Direction::Up),
                            via: vec![],
                        });
                        state = ReadState::Track;

                    } else {
                        panic!("[{line_no}] Unexpected line: {}", line);
                    }
                }
                ReadState::Map => {
                    if line.starts_with("cut") {

                        // A cut is a point on the map tracks cannot go onto
                        let parts = line.split_ascii_whitespace().collect::<Vec<_>>();
                        assert!(parts.len() == 3, "[{line_no}] In map mode, expected 'cut <x> <y>', got {line}");
                        
                        let x = parts[1].parse::<i32>().unwrap();
                        let y = parts[2].parse::<i32>().unwrap();

                        map.unwalkable.insert(Point { x, y });

                    } else if line.starts_with("element") {

                        let mut parts = line.split_ascii_whitespace();
                        parts.next().unwrap(); // Skip the 'element' keyword

                        let kind: Element = parts.next().unwrap().into();
                        let x = parts.next().unwrap().parse::<i32>().unwrap();
                        let y = parts.next().unwrap().parse::<i32>().unwrap();

                        let max_electrons = kind.max_electrons();
                        let free_electrons = parts
                            .next()
                            .and_then(|e| Some(e.parse::<u8>().unwrap()))
                            .unwrap_or(max_electrons);

                        map.atoms.push(Atom {
                            kind,
                            position: Point { x, y },
                            max_electrons,
                            free_electrons,
                        });
                        map.unwalkable.insert(Point { x, y });

                    } else if line.starts_with("end") {

                        // End of the map block, go back to parsing other blocks
                        state = ReadState::Default;

                    } else {
                        panic!("[{line_no}] In map mode, unexpected line: {}", line);
                    }
                
                },
                ReadState::Track => {
                    if line.starts_with("end") {
                        map.tracks.push(current_track.unwrap());
                        current_track = None;
                        state = ReadState::Default;
                        continue;    
                    }
                
                    let mut parts = line.split_ascii_whitespace();
                    parts.next().unwrap(); // Skip keyword

                    let x = parts.next().unwrap().parse::<i32>().unwrap();
                    let y = parts.next().unwrap().parse::<i32>().unwrap();

                    if line.starts_with("from") {

                        let direction = parts.next().unwrap().into();
                        assert!(parts.next().is_none());
                        current_track.as_mut().unwrap().from = (Point { x, y }, direction);

                    } else if line.starts_with("to") {

                        let direction = parts.next().unwrap().into();
                        assert!(parts.next().is_none());
                        current_track.as_mut().unwrap().to = (Point { x, y }, direction);

                    } else if line.starts_with("via") {

                        assert!(parts.next().is_none());
                        current_track.as_mut().unwrap().via.push(Point { x, y });

                    } else {
                        panic!("[{line_no}] In track mode, unexpected line: {}", line);
                    }
                },
            }
        }

        map
    }
}

#[derive(Clone, Debug)]
struct Molecule {
    position: Point,
    atoms: FxHashMap<Uuid, Atom>,
    bonds: FxHashMap<(Uuid, Uuid), u8>,
}

impl Molecule {
    fn free_electrons(&self) -> u8 {
        self.atoms.iter().map(|(_, a)| a.free_electrons).sum()
    }

    fn intersect(&self, other: &Molecule) -> bool {
        for (_, atom) in &self.atoms {
            for (_, other_atom) in &other.atoms {
                if self.position + atom.position == other_atom.position + other.position {
                    return true;
                }
            }
        }

        false
    }
}

#[derive(Clone, Debug)]
struct Bond {
    from: Uuid,
    to: Uuid,
    count: u8,
}

struct Simulation {
    molecules: FxHashMap<Uuid, Molecule>,
    tracks: Vec<Track>,
}

impl Simulation {
    fn new(map: &Map, tracks: &Vec<Track>) -> Self {
        let mut molecules = FxHashMap::default();
        let tracks = tracks.clone();

        // Create a molecule from each initial (non-ion) atom
        for atom in map.atoms.iter() {
            if let Element::Ion = atom.kind {
                continue;
            }

            let id = Uuid::new_v4();
            
            let mut atoms = FxHashMap::default();
            atoms.insert(Uuid::new_v4(), atom.clone());

            let bonds = FxHashMap::default();

            molecules.insert(
                Uuid::new_v4(),
                Molecule { position: atom.position, atoms, bonds });
        }

        log::debug!("Initial molecules: {molecules:?}");

        Simulation { molecules, tracks }
    }

    fn run(&mut self) -> Result<Vec<Molecule>, ()> {
        let mut track_positions = vec![];
        for track in self.tracks.iter() {
            track_positions.push(track.from.0);
        }

        // Run the simulation until all tracks 
        let mut tick = 0;
        loop {
            tick += 1;
            log::debug!("[simulate] Tick {tick}, currently at {track_positions:?}");

            // If all tracks are done, stop
            if self.tracks.iter().all(|track| track.via.len() == tick) {
                break;
            }

            // Each track moves forward one
            for (i, track) in self.tracks.iter().enumerate() {
                let src = track_positions[i];

                let dst = if track.via.len() < tick {
                    continue; // this track has finished
                } else if self.tracks[i].via.len() == tick {
                    track.to.0
                } else {
                    track.via[tick - 1]
                };

                let offset = dst - src;
                log::debug!("- Moving on track {i} from {src:?} to {dst:?}");

                
                let molecule_index = self
                    .molecules
                    .iter()
                    .find_map(
                        |(id, m)| 
                        if m.position == src {
                            Some(id)
                        } else {
                            None
                        })
                    .unwrap();
                log::debug!("- Molecule is at index {molecule_index:?}");
                
                self.try_move(molecule_index, offset)?;
                self.bond();

                let next = *track.via.last().unwrap() + track.to.1.into();
                track_positions[i] = next;
            }

        }

        Ok(self.molecules.values().cloned().collect())
    }


    // Helper function to move a molecule to a new position
    fn try_move(&mut self, index: &Uuid, offset: Point) -> Result<(), ()> {
        log::debug!("Trying to move molecule {index:?} by {offset:?}");

        // A molecule cannot hit another molecule
        for (other_index, other_molecule) in self.molecules.iter() {
            if other_index == index {
                continue;
            }

            if self.molecules[index].intersect(other_molecule) {
                return Err(());
            }
        }

        // If we passed all of the tests, move it!
        self.molecules.get_mut(index).unwrap().position += offset;

        // All is well
        Ok(())
    }

    // Helper function to try to bond any molecules that are touching
    fn bond(&mut self) { 
        // For each pair of molecules
        'bonding: loop {
            // let candidate = None;
            
            // Find a candidate bond
            for (i, molecule) in self.molecules.iter().enumerate() {
                for (j, other_molecule) in self.molecules.iter().enumerate() {
                }
            }

            break;
        }
    }
}

impl State<Map, ()> for Vec<Track> {
    fn next_states(&self, global: &Map) -> Option<Vec<(i64, (), Self)>>
    where
        Self: Sized {

        let mut next_states = vec![];
        
        for (i, track) in self.iter().enumerate() {
            for direction in &[Direction::Up, Direction::Down, Direction::Left, Direction::Right] {
                let next = *track.via.last().unwrap() + (*direction).into();
                if global.unwalkable.contains(&next) {
                    continue;
                }

                let mut new_state = self.clone();
                new_state[i].via.push(next);
                next_states.push((1, (), new_state));
            }
        }

        if next_states.is_empty() {
            None
        } else {
            Some(next_states)
        }
    }

    fn heuristic(&self, _: &Map) -> i64 {
        return 0;
    }

    fn is_valid(&self, _: &Map) -> bool {
        return true;
    }

    fn is_solved(&self, map: &Map) -> bool {
        if !self
            .iter()
            .all(|track| *track.via.last().unwrap() + track.to.1.into() == track.to.0) {
            return false;
        }

        log::debug!("Checking if all atoms are bonded");

        let simulation = Simulation::new(map, self);
        
        let molecules = simulation.run();
        if molecules.is_err() {
            return false;
        }
        let molecules = molecules.unwrap();

        if molecules.iter().all(|m| m.free_electrons() == 0) {
            return true;
        }

        false
    }

    fn stringify(&self, _: &Map) -> String {
        todo!()
    }
}

fn main() {
    env_logger::init();

    let input = io::read_to_string(io::stdin()).unwrap();
    let map = Map::from(input.as_str());
    let tracks = map.tracks.clone();

    let mut solver = Solver::new(map.clone(), tracks);

    while let Some(state) = solver.next() {
        if solver.states_checked() % 100000 != 0 {
            continue;
        }
        log::info!("{solver}, state:\n{}", state.stringify(&map));
    }

    let solution = solver.get_solution();
    if solution.is_none() {
        log::info!("No solution found");
        return;
    }
    let solution = solution.unwrap();

    println!("{solution:?}");

}