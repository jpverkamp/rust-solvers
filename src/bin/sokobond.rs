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
}

impl Map {
    fn load(input: &str) -> (Map, LocalState) {
        let mut width = 0;
        let mut height = 0;
        let mut walls = Vec::new();

        let mut primary = None;
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
                        let molecule = Molecule {
                            offset: pt,
                            elements: vec![(e, Point::ZERO, e.free_electrons())],
                        };

                        if c.is_uppercase() {
                            if primary.is_some() {
                                panic!("Multiple primary molecules");
                            }

                            primary = Some(molecule);
                        } else {
                            molecules.push(molecule);
                        }

                        walls.push(false);
                    }
                    Err(_) => walls.push(c != '-'),
                }
            }
        }

        (
            Map {
                width,
                height,
                walls,
            },
            LocalState {
                primary: primary.expect("Must have a primary"),
                others: molecules,
            },
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
    Oxygen,
}

impl TryFrom<char> for Element {
    type Error = ();

    fn try_from(value: char) -> Result<Self, Self::Error> {
        use Element::*;

        match value {
            'h' | 'H' => Ok(Hydrogen),
            'o' | 'O' => Ok(Oxygen),
            _ => Err(()),
        }
    }
}

impl Into<char> for &Element {
    fn into(self) -> char {
        match self {
            Element::Hydrogen => 'H',
            Element::Oxygen => 'O',
        }
    }
}

impl Element {
    fn free_electrons(&self) -> usize {
        use Element::*;

        match self {
            Hydrogen => 1,
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
    // How many free electrons in the hole molecule
    fn free_electrons(&self) -> usize {
        self.elements.iter().map(|(_, _, free)| *free).sum::<usize>()
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
                // Not adjacent
                if src_offset.manhattan_distance(offset + *dst_offset) != 1 {
                    continue;
                }

                // Not enough free electrons
                if *dst_free == 0 {
                    continue;
                }

                // We'll automatically bind as many free electrons as we can
                let bind_electrons = *src_free.min(dst_free);

                // We're good! Bind it, using up that many electrons from each
                *src_free -= bind_electrons;
                *dst_free -= bind_electrons;
                bound = true;
            }
        }

        // If we bound anything, add the other elements to our molecule
        if bound {
            for (element, element_offset, free_electrons) in other.elements {
                self.elements.push((element, offset + element_offset, free_electrons));
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

        assert!(a.intersects(Point(0, 0), &b));
        assert!(!a.intersects(Point(0, 1), &b));
    }

    #[test]
    fn test_bind() {
        let mut a = Molecule {
            offset: Point::ZERO,
            elements: vec![(Element::Hydrogen, Point::ZERO, 1)],
        };

        let b = Molecule {
            offset: Point(1, 0),
            elements: vec![(Element::Hydrogen, Point::ZERO, 1)],
        };

        let bound = a.try_bind(Point(1, 0), &b);

        assert!(bound);
        assert_eq!(a.elements.len(), 2);
        assert_eq!(a.elements[0].2, 0);
    }

    #[test]
    fn test_nobind_no_free() {
        let mut a = Molecule {
            offset: Point::ZERO,
            elements: vec![(Element::Hydrogen, Point::ZERO, 1)],
        };

        let b = Molecule {
            offset: Point(1, 0),
            elements: vec![(Element::Hydrogen, Point::ZERO, 0)],
        };

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
    fn test_double_bind_o2() {
        let mut a = Molecule {
            offset: Point::ZERO,
            elements: vec![(Element::Oxygen, Point::ZERO, 2)],
        };

        let b = Molecule {
            offset: Point(1, 0),
            elements: vec![(Element::Oxygen, Point::ZERO, 2)],
        };

        let bound = a.try_bind(Point(1, 0), &b);

        assert!(bound);
        assert_eq!(a.elements.len(), 2);
        assert_eq!(a.elements[0].2, 0);
        assert_eq!(a.elements[1].2, 0);
    }
}

// The local state will be a primary molecule and a list of other molecules
// The primary is the one that can move
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct LocalState {
    primary: Molecule,
    others: Vec<Molecule>,
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
        // TODO: If we have no free electrons (and we're not solved), we're invalid

        return true;
    }

    fn is_solved(&self, _global: &Map) -> bool {
        // The puzzle is solved once there is no longer any free electrons
        self.primary.free_electrons() == 0 && self.others.iter().all(|m| m.free_electrons() == 0)
    }

    fn next_states(&self, map: &Map) -> Option<Vec<(i64, Step, LocalState)>> {
        let mut next_states = Vec::new();

        // Try to move the primary each direction
        'next_step: for step in [Step::North, Step::South, Step::East, Step::West].iter() {
            // Check for intersections

            // Would hit a wall
            if self.primary.intersects_wall((*step).into(), map) {
                log::debug!("intersected {:?} with offset {step:?}", self.primary);
                log::debug!("{step:?} failed, wall"); 
                continue 'next_step;
            }

            // Would hit another molecule
            for molecule in self.others.iter() {
                if self.primary.intersects((*step).into(), molecule) {
                    log::debug!("{step:?} failed, molecule"); 
                    continue 'next_step;
                }
            }

            // Move is allowed, update primary
            let mut new_state = self.clone();
            new_state.primary.offset = new_state.primary.offset + (*step).into();

            // Try to bind with each other molecule
            let mut bound_indexes = Vec::new();

            for (i, other) in new_state.others.iter().enumerate() {
                if new_state.primary.try_bind(other.offset - new_state.primary.offset, other) {
                    bound_indexes.push(i);
                }
            }

            for i in bound_indexes.iter().rev() {
                new_state.others.remove(*i);
            }

            // Valid state, queue it
            next_states.push((1, *step, new_state));
        }

        // if cfg!(debug_assertions) {} {
        //     log::debug!("--- next states generated:");
        //     for (score, step, state) in &next_states {
        //         log::debug!("{}: {:?}", score, step);
        //         if cfg!(debug_assertions) {
        //             println!("{}", state.to_string(map)));
        //         } 
        //     }
        // }

        if next_states.is_empty() {
            return None;
        } else {
            return Some(next_states);
        }
    }

    fn heuristic(&self, _global: &Map) -> i64 {
        return 0;
    }

    fn to_string(&self, map: &Map) -> String {
        let mut grid = vec![vec![' '; map.width]; map.height];

        // DEBUG
        // if cfg!(debug_assertions) {
        //     output.push_str(format!("{}x{}: {:?}", map.width, map.height, self).as_str());
        // }

        for (element, offset, _) in &self.primary.elements {
            let offset = self.primary.offset + *offset;
            grid[offset.1 as usize][offset.0 as usize] = element.into();
        }

        for molecule in self.others.iter() {
            for (element, offset, _) in &molecule.elements {
                let offset = molecule.offset + *offset;
                let c: char = element.into();
                grid[offset.1 as usize][offset.0 as usize] = c.to_ascii_lowercase();
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

fn main() {
    env_logger::init();

    let input = io::read_to_string(io::stdin()).unwrap();
    let (map, molecules) = Map::load(&input);

    println!("Initial state:");
    println!("{}", molecules.to_string(&map));

    let mut solver = Solver::new(map.clone(), molecules.clone());

    while let Some(state) = solver.next() {
        if solver.states_checked() % 100000 != 0 {
            continue;
        }

        println!("===== ===== ===== ===== =====");
        println!("{}", state.to_string(&map));
        println!(
            "{} states checked, {} in queue, {} invalidated, {} seconds, heuristic: {}",
            solver.states_checked(),
            solver.in_queue(),
            solver.states_invalidated(),
            solver.time_spent(),
            state.heuristic(&map),
        );

        if solver.states_checked() > 100 {
            break;
        }
    }

    let solution = solver.get_solution();
    if solution.is_none() {
        println!("no solution found");
        return;
    }
    let solution = solution.unwrap();

    println!("Solved state:");
    println!("{}", solution.to_string(&map));

    let mut steps = String::new();
    for step in solver.path(&molecules, &solution).unwrap() {
        match step {
            Step::North => steps += "W",
            Step::South => steps += "S",
            Step::East => steps += "D",
            Step::West => steps += "A",
        }
    }
    println!("path: {}", steps);

    println!(
        "{} states, {} seconds",
        solver.states_checked(),
        solver.time_spent()
    );
}

// TODO: Rebecca
