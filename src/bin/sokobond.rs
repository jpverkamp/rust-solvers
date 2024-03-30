use std::io;
use std::ops::Add;
use std::ops::Sub;

use solver::{Solver, State};

// A point in 2D space
#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
struct Point(isize, isize);

impl Add<Point> for Point {
    type Output = Point;

    fn add(self, rhs: Point) -> Self::Output {
        Point(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl Sub<Point> for Point {
    type Output = Point;

    fn sub(self, rhs: Point) -> Self::Output {
        Point(self.0 - rhs.0, self.1 - rhs.1)
    }
}

impl Point {
    const ZERO: Point = Point(0, 0);

    fn manhattan_distance(&self, other: Point) -> isize {
        (self.0 - other.0).abs() + (self.1 - other.1).abs()
    }
}

// Global state: a set of walls
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct Map {
    width: usize,
    height: usize,
    walls: Vec<bool>,
    splitters: Vec<Point>,
}

impl Map {
    fn load(input: &str) -> (Map, LocalState) {
        let mut width = 0;
        let mut height = 0;
        let mut walls = Vec::new();
        let mut splitters = Vec::new();

        let mut molecules = Vec::new();

        for (y, line) in input.lines().enumerate() {
            width = width.max(line.len());
            height += 1;

            for (x, c) in line.chars().enumerate() {
                let pt = Point(x as isize, y as isize);

                match Element::try_from(c) {
                    // A new element, convert it to a molecule
                    // Primaries are uppercase and put at the start of the list
                    // The rest are added to the end
                    // Technically only one primary is supported, this will take the last if multiple
                    Ok(e) => {
                        let molecule = Molecule::new(pt, e);

                        if c.is_uppercase() {
                            molecules.insert(0, molecule);
                        } else {
                            molecules.push(molecule);
                        }

                        walls.push(false);
                    }
                    Err(_) => {
                        match c {
                            // Splitters are offset between the grid lines
                            '\\' => {
                                walls.push(false);
                                splitters.push(pt);
                            }
                            ' ' | '-' => walls.push(false),
                            'x' | 'X' | '#' => walls.push(true),
                            _ => panic!("unknown character: {}", c),
                        }
                    }
                }
            }
        }

        // Try to bond each original pair of molecules
        'settled: loop {
            for i in 0..molecules.len() {
                for j in (i+1)..molecules.len() {
                    let mut primary = molecules[i].clone();

                    if primary.try_bind(Point::ZERO, &molecules[j]) {
                        molecules[i] = primary;
                        molecules.remove(j);
                        continue 'settled;
                    }
                }
            }
            break;
        }

        (
            Map {
                width,
                height,
                walls,
                splitters,
            },
            LocalState { molecules },
        )
    }

    fn is_wall(&self, x: usize, y: usize) -> bool {
        if x >= self.width || y >= self.height {
            return true;
        }

        self.walls[y * self.width + x]
    }
}

// Represents single elements
#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
enum Element {
    Hydrogen,
    Helium,
    Nitrogen,
    Carbon,
    Oxygen,
}

impl TryFrom<char> for Element {
    type Error = ();

    fn try_from(value: char) -> Result<Self, Self::Error> {
        use Element::*;

        match value {
            'h' | 'H' => Ok(Hydrogen),
            'e' | 'E' => Ok(Helium),
            'n' | 'N' => Ok(Nitrogen),
            'c' | 'C' => Ok(Carbon),
            'o' | 'O' => Ok(Oxygen),
            _ => Err(()),
        }
    }
}

impl Into<char> for &Element {
    fn into(self) -> char {
        match self {
            Element::Hydrogen => 'H',
            Element::Helium => 'E',
            Element::Nitrogen => 'N',
            Element::Carbon => 'C',
            Element::Oxygen => 'O',
        }
    }
}

impl Element {
    fn free_electrons(&self) -> usize {
        use Element::*;

        match self {
            Hydrogen => 1,
            Helium => 0,
            Nitrogen => 3,
            Carbon => 4,
            Oxygen => 2,
        }
    }
}

// Represent elements joined together into a molecule
// Each element contains the Element, offset, and remaining free electrons
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct Molecule {
    offset: Point,
    elements: Vec<(Element, Point, usize)>,
}

impl Molecule {
    fn new(offset: Point, element: Element) -> Molecule {
        Molecule {
            offset,
            elements: vec![(element, Point::ZERO, element.free_electrons())],
        }
    }

    // Test for molecular helium
    fn is_helium(&self) -> bool {
        self.elements.len() == 1 && self.elements[0].0 == Element::Helium
    }

    // How many free electrons in the hole molecule
    fn free_electrons(&self) -> usize {
        self.elements
            .iter()
            .map(|(_, _, free)| *free)
            .sum::<usize>()
    }

    // If the given molecule at an offset would intersect with a wall
    fn intersects_wall(&self, offset: Point, map: &Map) -> bool {
        for (_, element_offset, _) in &self.elements {
            let target = self.offset + *element_offset + offset;
            if map.is_wall(target.0 as usize, target.1 as usize) {
                return true;
            }
        }

        false
    }

    // If the given molecule + an offset would intersect another molecule
    fn intersects(&self, offset: Point, other: &Molecule) -> bool {
        for (_, element_offset, _) in &self.elements {
            for (_, other_element_offset, _) in &other.elements {
                if self.offset + *element_offset + offset == other.offset + *other_element_offset {
                    return true;
                }
            }
        }

        false
    }

    // Try to bind two molecules together
    // Offset is between the centers of the molecules
    // Updates and returns true if successful
    fn try_bind(&mut self, offset: Point, other: &Molecule) -> bool {
        let mut bound = false;

        // Make local mutable copies
        let mut other = other.clone();

        // Go through each molecule pairwise
        for (_, src_offset, src_free) in self.elements.iter_mut() {
            // Skip our elements that are no longer free
            if *src_free == 0 {
                continue;
            }

            for (_, dst_offset, dst_free) in other.elements.iter_mut() {
                let src_pt = self.offset + offset + *src_offset;
                let dst_pt = other.offset + *dst_offset;

                // Not adjacent
                if src_pt.manhattan_distance(dst_pt) != 1 {
                    continue;
                }

                // Not enough free electrons
                if *dst_free == 0 {
                    continue;
                }

                // Bind the two elements
                *src_free -= 1;
                *dst_free -= 1;
                bound = true;
            }
        }

        // If we bound anything, add the other elements to our molecule
        if bound {
            for (element, element_offset, free_electrons) in other.elements {
                self.elements.push((
                    element,
                    other.offset - self.offset + element_offset,
                    free_electrons,
                ));
            }
            true
        } else {
            false
        }
    }
}

#[cfg(test)]
mod test_molecule {
    use super::*;

    #[test]
    fn test_wall_intersection_hit() {
        let a = Molecule {
            offset: Point::ZERO,
            elements: vec![(Element::Hydrogen, Point::ZERO, 1)],
        };

        let map = Map {
            width: 3,
            height: 3,
            walls: vec![
                true, true, true, true, false, true, true, true, true, // 3x3
            ],
            splitters: Vec::new(),
        };

        assert!(a.intersects_wall(Point(1, 0), &map));
        assert!(!a.intersects_wall(Point(1, 1), &map));
    }

    #[test]
    fn test_molecule_intersection() {
        let a = Molecule {
            offset: Point::ZERO,
            elements: vec![(Element::Hydrogen, Point::ZERO, 1)],
        };

        let b = Molecule {
            offset: Point(1, 0),
            elements: vec![(Element::Hydrogen, Point::ZERO, 1)],
        };

        assert!(a.intersects(Point(1, 0), &b));
        assert!(!a.intersects(Point(0, 0), &b));
        assert!(!a.intersects(Point(0, 1), &b));
    }

    #[test]
    fn test_bind() {
        let mut a = Molecule::new(Point::ZERO, Element::Hydrogen);
        let b = Molecule::new(Point(1, 0), Element::Hydrogen);

        let bound = a.try_bind(Point::ZERO, &b);

        assert!(bound);
        assert_eq!(a.elements.len(), 2);
        assert_eq!(a.elements[0].2, 0);
    }

    #[test]
    fn test_nobind_no_free() {
        let mut a = Molecule::new(Point::ZERO, Element::Hydrogen);
        let b = Molecule::new(Point(1, 0), Element::Helium);

        let bound = a.try_bind(Point(1, 0), &b);

        assert!(!bound);
        assert!(a.elements[0].2 == 1);
    }

    #[test]
    fn test_nobind_too_far() {
        let mut a = Molecule {
            offset: Point::ZERO,
            elements: vec![(Element::Hydrogen, Point::ZERO, 1)],
        };

        let b = Molecule {
            offset: Point(2, 0),
            elements: vec![(Element::Hydrogen, Point::ZERO, 1)],
        };

        let bound = a.try_bind(Point(2, 0), &b);

        assert!(!bound);
        assert!(a.elements[0].2 == 1);
    }

    #[test]
    fn test_single_bind_o2_not_double() {
        let mut a = Molecule {
            offset: Point::ZERO,
            elements: vec![(Element::Oxygen, Point::ZERO, 2)],
        };

        let b = Molecule {
            offset: Point(1, 0),
            elements: vec![(Element::Oxygen, Point::ZERO, 2)],
        };

        let bound = a.try_bind(Point::ZERO, &b);

        assert!(bound);
        assert_eq!(a.elements.len(), 2);
        assert_eq!(a.elements[0].2, 1);
        assert_eq!(a.elements[1].2, 1);
    }
}

// The local state will be a primary molecule and a list of other molecules
// The primary is the one that can move
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct LocalState {
    molecules: Vec<Molecule>,
}

impl LocalState {
    // Try to move the ith molecule by the given offset
    // Will also move all other touching molecules out of the way
    // Returns the new state if successful
    fn try_move(&mut self, map: &Map, index: usize, offset: Point) -> bool {
        if index > self.molecules.len() {
            return false;
        }

        // Collect all molecules that would move
        let mut would_move = vec![index];

        // Recursively add any molecules we'd hit
        // Each will try to move by the same offset (represent a chain of moving)
        'settling: loop {
            // For each currently moving molecule
            for i in 0..self.molecules.len() {
                if !would_move.contains(&i) {
                    continue;
                }

                // If moving that molecule into another would cause an intersection, add it to the list
                for j in 0..self.molecules.len() {
                    if i == j || would_move.contains(&j) {
                        continue;
                    }

                    if self.molecules[i].intersects(offset, &self.molecules[j]) {
                        would_move.push(j);
                        continue 'settling;
                    }
                }
            }

            break;
        }

        // Check each moving molecule to see if it would hit a wall
        for i in would_move.iter() {
            if self.molecules[*i].intersects_wall(offset, map) {
                return false;
            }
        }

        // Assume move is allowed, update all moving molecules
        for i in would_move.iter() {
            self.molecules[*i].offset = self.molecules[*i].offset + offset;
        }

        // Try to bind all pairs of molecules
        // TODO: Inefficient, we only really need to check moved (?); especially with all the cloning
        'settling: loop {
            for i in 0..self.molecules.len() {
                for j in 0..self.molecules.len() {
                    if i == j {
                        continue;
                    }

                    let mut primary = self.molecules[i].clone();

                    if primary.try_bind(Point::ZERO, &self.molecules[j]) {
                        self.molecules[i] = primary;
                        self.molecules.remove(j);
                        continue 'settling;
                    }
                }
            }

            break;
        }

        true
    }
}

#[cfg(test)]
mod test_localstate {
    #[test]
    fn test_move_basic() {
        use super::*;

        let (map, mut state) = Map::load("H-#");

        assert!(state.try_move(&map, 0, Point(1, 0)));
        assert_eq!(state.molecules[0].offset, Point(1, 0));
    }

    #[test]
    fn test_bump_into_wall() {
        use super::*;

        let (map, mut state) = Map::load("-H#");

        assert!(!state.try_move(&map, 0, Point(1, 0)));
    }

    #[test]
    fn test_move_push() {
        use super::*;

        let (map, mut state) = Map::load("Ee-#");

        assert!(state.try_move(&map, 0, Point(1, 0)));
        assert_eq!(state.molecules.len(), 2);
        assert_eq!(state.molecules[0].offset, Point(1, 0));
        assert_eq!(state.molecules[1].offset, Point(2, 0));
    }

    #[test]
    fn test_move_and_bind() {
        use super::*;

        let (map, mut state) = Map::load("H-h#");

        assert!(state.try_move(&map, 0, Point(1, 0)));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].offset, Point(1, 0));
        assert_eq!(state.molecules[0].elements.len(), 2);

        assert_eq!(state.to_string(&map), "-HH#\n".to_string())
    }

    #[test]
    fn test_push_and_bind() {
        use super::*;

        let (map, mut state) = Map::load("Hh-#");

        assert!(state.try_move(&map, 0, Point(1, 0)));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].offset, Point(1, 0));
        assert_eq!(state.molecules[0].elements.len(), 2);
    }

    #[test]
    fn test_no_move_push_into_wall() {
        use super::*;

        let (map, mut state) = Map::load("Ee#");

        assert!(!state.try_move(&map, 0, Point(1, 0)));
        assert_eq!(state.molecules[0].offset, Point(0, 0));
        assert_eq!(state.molecules[1].offset, Point(1, 0));
    }

    #[test]
    fn test_push_unbound() {
        use super::*;

        let (map, mut state) = Map::load("-H-\no--\n-h-\n---");

        assert!(state.try_move(&map, 0, Point(0, 1)));
        assert!(state.try_move(&map, 0, Point(0, 1)));

        println!("{}", state.to_string(&map));

        // Two molecules
        assert_eq!(state.molecules.len(), 2);

        // First is OH with a free still on the O
        assert_eq!(state.molecules[0].offset, Point(1, 2));
        assert_eq!(state.molecules[0].elements[0].2, 0);
        assert_eq!(state.molecules[0].elements[1].2, 1);

        // Second is h with a free still open
        assert_eq!(state.molecules[1].offset, Point(1, 3));
        assert_eq!(state.molecules[1].elements[0].2, 1);
    }
}

// The step the primary molecule takes each tick
#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
enum Step {
    North,
    South,
    East,
    West,
}

impl Into<Point> for Step {
    fn into(self) -> Point {
        match self {
            Step::North => Point(0, -1),
            Step::South => Point(0, 1),
            Step::East => Point(1, 0),
            Step::West => Point(-1, 0),
        }
    }
}

impl State<Map, Step> for LocalState {
    fn is_valid(&self, _global: &Map) -> bool {
        if self.is_solved(_global) {
            return true;
        }

        // The primary molecule must have free electrons
        if !(self.molecules[0].is_helium() || self.molecules[0].free_electrons() > 0) {
            return false;
        }

        // Any non-primary (other than helium) must have free electrons
        if self.molecules.iter().skip(1).any(|m| !m.is_helium() && m.free_electrons() == 0) {
            return false;
        }

        return true;
    }

    fn is_solved(&self, _global: &Map) -> bool {
        self.molecules.iter().all(|m| m.is_helium() || m.free_electrons() == 0)

    }

    fn next_states(&self, map: &Map) -> Option<Vec<(i64, Step, LocalState)>> {
        let mut next_states = Vec::new();

        // Try to move the primary each direction
        for step in [Step::North, Step::South, Step::East, Step::West].iter() {
            let mut next_state = self.clone();
            if next_state.try_move(map, 0, (*step).into()) {
                next_states.push((1, *step, next_state));
            }
        }

        if next_states.is_empty() {
            return None;
        } else {
            return Some(next_states);
        }
    }

    fn heuristic(&self, _global: &Map) -> i64 {
        self.molecules.iter().map(|m| m.free_electrons() as i64).sum()
    }

    fn to_string(&self, map: &Map) -> String {
        let mut grid = vec![vec!['-'; map.width]; map.height];

        // DEBUG
        // if cfg!(debug_assertions) {
        //     output.push_str(format!("{}x{}: {:?}", map.width, map.height, self).as_str());
        // }

        for (i, molecule) in self.molecules.iter().enumerate() {
            for (element, offset, _) in &molecule.elements {
                let offset = molecule.offset + *offset;
                let mut c: char = element.into();
                if i > 0 {
                    c = c.to_ascii_lowercase();
                }
                grid[offset.1 as usize][offset.0 as usize] = c;
            }
        }

        for y in 0..map.height {
            for x in 0..map.width {
                if map.walls[y * map.width + x] {
                    grid[y][x] = '#';
                }
            }
        }

        let mut output = String::new();
        for y in 0..map.height {
            for x in 0..map.width {
                output.push(grid[y][x]);
            }
            output.push('\n');
        }
        output
    }
}

fn solve(input: &str) -> Option<String> {
    let (map, molecules) = Map::load(input);

    log::info!("Initial state:\n{}", molecules.to_string(&map));

    let mut solver = Solver::new(map.clone(), molecules.clone());

    while let Some(state) = solver.next() {
        log::info!(
            "{} states checked, {} in queue, {} invalidated, {} seconds, heuristic: {}, state:\n{}\n",
            solver.states_checked(),
            solver.in_queue(),
            solver.states_invalidated(),
            solver.time_spent(),
            state.heuristic(&map),
            state.to_string(&map),
        );
    }

    let solution = solver.get_solution();
    if solution.is_none() {
        return None;
    }
    let solution = solution.unwrap();

    log::info!("Solved after {} states in {} seconds:\n{}",
        solver.states_checked(),
        solver.time_spent(),
        solution.to_string(&map),
    );

    let mut steps = String::new();
    for step in solver.path(&molecules, &solution).unwrap() {
        match step {
            Step::North => steps += "W",
            Step::South => steps += "S",
            Step::East => steps += "D",
            Step::West => steps += "A",
        }
    }
    Some(steps)
}

fn main() {
    env_logger::init();

    let input = io::read_to_string(io::stdin()).unwrap();
    let solution = solve(&input);

    if solution.is_none() {
        panic!("no solution found");
    }
    
    let solution = solution.unwrap();
    println!("{}", solution);
}

#[cfg(test)]
mod test_solutions {
    use super::*;

    macro_rules! test {
        ($name:ident, $folder:expr, $file:expr, $expected:expr) => {
            #[test]
            fn $name() {
                let input = include_str!(concat!("../../data/sokobond/", $folder, "/", $file));
                let solution = solve(input);
                assert_eq!(solution.expect("solution exists"), $expected);
            }            
        };
    }

    test!{test_01_01, "01 - Yellow", "01 - Let's Go.txt", "WWDDWWDD"}
    test!{test_01_02, "01 - Yellow", "02 - Cell.txt", "DWDDAA"}
    test!{test_01_03, "01 - Yellow", "03 - Loop.txt", "WWDDSSWWAAASSSSD"}
    test!{test_01_04, "01 - Yellow", "04 - Monty Hall.txt", "WSDDWSAAWWWDAAA"}
    test!{test_01_05, "01 - Yellow", "05 - Knot.txt", "WDSASAWDDDAAAA"}
    test!{test_01_06, "01 - Yellow", "06 - Push.txt", "DDWDSD"}
    test!{test_01_07, "01 - Yellow", "07 - Structure.txt", "DSSWAAA"}
    test!{test_01_08, "01 - Yellow", "08 - Lotus.txt", "SDSSWWAADSSADSSWAWA"}
    test!{test_01_09, "01 - Yellow", "09 - Coyote.txt", "WWAASASDDWWAASSSAWWSSDDDWWWAASSSAAWW"}
    test!{test_01_10, "01 - Yellow", "10 - Roadrunner.txt", "AWWDDSSWWDDSDAWAASSASDWWWDDSDSWAAASS"}

    test!{test_02_01, "02 - Orange", "01 - Suit.txt", "SSDDDWWSSAAASSWWWWDDDSSSSA"}
    test!{test_02_02, "02 - Orange", "02 - Chimney.txt", "WAWSSDDWAAADWSS"}
    test!{test_02_03, "02 - Orange", "03 - Heart.txt", "WWDWASAAWDAWWDSSASDWDDWW"}
    test!{test_02_04, "02 - Orange", "04 - Factory.txt", "AASSDDDDDWWWWSSSSAAAAAWWWWWDSWD"}
    test!{test_02_05, "02 - Orange", "05 - Scoop.txt", "DAASSDSDWWSAAAWWDDSSDDW"}
    test!{test_02_06, "02 - Orange", "06 - Block.txt", "ASSSDDWSDDWWWADSSSAAWAW"}
    test!{test_02_07, "02 - Orange", "07 - Chandelier.txt", "DSSWDSSAAWDWWASSWDDSSWWAASWDDSSA"}
    test!{test_02_08, "02 - Orange", "08 - Creature.txt", "WAAWWDDWDSWAASAASSDDSSWWDDSDW"}
    test!{test_02_09, "02 - Orange", "09 - Rosie.txt", "WDWWWSDSSAAASDDDWWWWAA"}
    test!{test_02_10, "02 - Orange", "10 - Kruskal.txt", "DWDSASAWAWDSSAWWDDSSAAWWWSSSDDWWWSSS"}

    test!{test_03_01, "03 - Gray", "01 - Helium.txt", "WDDDDDSAAWASSASSDWWWSSSS"}
    test!{test_03_02, "03 - Gray", "02 - Tee.txt", "WWDDSAWAASWDDSAWAAAASDWDDSWASWDDSSSSADWWWASDSSSAWDWW"}
    test!{test_03_03, "03 - Gray", "03 - Freedom.txt", "WWWSSSDWWSSAAWWSSAWSDDWWWADSSSAAWAWDWD"}
    // test!{test_03_04, "03 - Gray", "04 - Against the Wall.txt", ""}
}