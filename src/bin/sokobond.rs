use std::io;
use std::ops::Add;
use std::ops::Sub;

use solver::{Solver, State};

const SINGLE_HORIZONTAL: char = '-';
const SINGLE_VERTICAL: char = '|';
const DOUBLE_HORIZONTAL: char = '=';
const DOUBLE_VERTICAL: char = '‖';
const TRIPLE_HORIZONTAL: char = '≡';
const TRIPLE_VERTICAL: char = '⦀';

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

                    if primary.try_bond(Point::ZERO, &molecules[j]) {
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
    bonds: Vec<(Point, Point, usize)>,
}

impl Molecule {
    fn new(offset: Point, element: Element) -> Molecule {
        Molecule {
            offset,
            elements: vec![(element, Point::ZERO, element.free_electrons())],
            bonds: Vec::new(),
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

    // Try to bond two molecules together
    // Offset is between the centers of the molecules
    // Updates and returns true if successful
    fn try_bond(&mut self, offset: Point, other: &Molecule) -> bool {
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

                // bond the two elements
                bound = true;

                self.bonds.push((
                    offset + *src_offset,
                    other.offset - self.offset + *dst_offset, 
                    1
                ));

                *src_free -= 1;
                *dst_free -= 1;
            }
        }

        // If we bound anything, add the other elements and bonds to our molecule
        if bound {
            for (element, element_offset, free_electrons) in other.elements {
                self.elements.push((
                    element,
                    other.offset - self.offset + element_offset,
                    free_electrons,
                ));
            }

            for (src, dst, count) in other.bonds {
                self.bonds.push((
                    other.offset - self.offset + src,
                    other.offset - self.offset + dst,
                    count,
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
        let a = Molecule::new(Point::ZERO, Element::Hydrogen);

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
        let a = Molecule::new(Point::ZERO, Element::Hydrogen);
        let b = Molecule::new(Point(1, 0), Element::Hydrogen);

        assert!(a.intersects(Point(1, 0), &b));
        assert!(!a.intersects(Point(0, 0), &b));
        assert!(!a.intersects(Point(0, 1), &b));
    }

    #[test]
    fn test_bond() {
        let mut a = Molecule::new(Point::ZERO, Element::Hydrogen);
        let b = Molecule::new(Point(1, 0), Element::Hydrogen);

        let bound = a.try_bond(Point::ZERO, &b);

        assert!(bound);
        assert_eq!(a.elements.len(), 2);
        assert_eq!(a.elements[0].2, 0);
    }

    #[test]
    fn test_nobond_no_free() {
        let mut a = Molecule::new(Point::ZERO, Element::Hydrogen);
        let b = Molecule::new(Point(1, 0), Element::Helium);

        let bound = a.try_bond(Point(1, 0), &b);

        assert!(!bound);
        assert!(a.elements[0].2 == 1);
    }

    #[test]
    fn test_nobond_too_far() {
        let mut a = Molecule::new(Point::ZERO, Element::Hydrogen);
        let b = Molecule::new(Point(2, 0), Element::Hydrogen);

        let bound = a.try_bond(Point(2, 0), &b);

        assert!(!bound);
        assert!(a.elements[0].2 == 1);
    }

    #[test]
    fn test_single_bond_o2_not_double() {
        let mut a = Molecule::new(Point::ZERO, Element::Oxygen);
        let b = Molecule::new(Point(1, 0), Element::Oxygen);

        let bound = a.try_bond(Point::ZERO, &b);

        assert!(bound);
        assert_eq!(a.elements.len(), 2);
        assert_eq!(a.elements[0].2, 1);
        assert_eq!(a.elements[1].2, 1);
    }
}

// The local state will be a primary molecule and a list of other molecules
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

        // Apply splitter before moving, since the split half doesn't move
        // TODO: Splitting multiple bonds
        let mut to_split = None;
        'splitting: for splitter in &map.splitters {
            for (i, (a, b, _count)) in self.molecules[index].bonds.iter().enumerate() {
                let a = self.molecules[index].offset + *a;
                let b = self.molecules[index].offset + *b;

                // Vertical bonds have the same x
                let is_vertical = a.0 == b.0;
                
                // We'll hit a vertical splitter if the offset is horizontal and we're moving across it
                // Ignore bonds that are moving the wrong way
                if is_vertical && offset.0 == 0 || !is_vertical && offset.1 == 0{
                    continue;
                }

                // Moving 'positive' is down or right
                let is_positive = offset.0 > 0 || offset.1 > 0;

                // Because either x or y is the same for all bonds, min is top/left and max is bottom/right
                // This will always match the splitter if we're moving across it right or down
                let pre_min = Point(a.0.min(b.0), a.1.min(b.1));

                // The post point is the one after we've moved
                let post_a = a + offset;
                let post_b = b + offset;
                let post_min = Point(post_a.0.min(post_b.0), post_a.1.min(post_b.1));
                
                // If we're moving positive, the min (top left) will equal the splitter
                if is_positive && splitter != &pre_min {
                    continue;
                }

                // If we're moving negative, then the *post* min will equal the splitter
                if !is_positive && splitter != &post_min {
                    continue;
                }

                // If we made it this far in the loop, this is the bond to break
                to_split = Some(i);
                break 'splitting;
            }
        };

        // If we found a bond to split, do that
        if let Some(i) = to_split {
            // Create a new copy of the molecule to modify
            let mut src = self.molecules[index].clone();

            // Add the count of the bond to the free electrons of both
            {
                let (bond_a, bond_b, count) = src.bonds[i];
                for (_, pt, free) in &mut src.elements {
                    if *pt == bond_a || *pt == bond_b {
                        *free += count;
                    }
                }
            }

            // Remove the bond
            src.bonds.remove(i);

            // Create the new molecule (we'll remove the halves from each)
            let mut dst = src.clone();

            

            // Now starting at the origin in src, remove anything we can get to in dst
            let mut connected_elements = Vec::new();
            let mut connected_bonds = Vec::new();

            let mut todo = vec![Point::ZERO];
            let mut done = vec![];
            
            while let Some(pt) = todo.pop() {
                done.push(pt);

                for (i, (_, offset, _)) in src.elements.iter().enumerate() {
                    if *offset == pt {
                        connected_elements.push(i);
                    }
                }

                for (i, (src, dst, _)) in src.bonds.iter().enumerate() {
                    if src == &pt {
                        connected_bonds.push(i);
                    }

                    if todo.contains(dst) || done.contains(dst) {
                        continue;
                    }

                    todo.push(*dst);
                }
            }

            // Now remove the connected elements and bonds from the dst
            for i in connected_elements.iter().rev() {
                dst.elements.remove(*i);
            }
            for i in connected_bonds.iter().rev() {
                dst.bonds.remove(*i);
            }

            // Update dst's offset so it has a 0,0
            let new_zero = dst.elements[0].1;
            dst.offset = dst.offset + new_zero;
            for (_, element_offset, _) in &mut dst.elements {
                *element_offset = *element_offset - new_zero;
            }

            // And remove everything else from src
            for i in (0..src.elements.len()).rev() {
                if !connected_elements.contains(&i) {
                    src.elements.remove(i);
                }
            }
            for i in (0..src.bonds.len()).rev() {
                if !connected_bonds.contains(&i) {
                    src.bonds.remove(i);
                }
            }

            // Now, update the original molecule and add the new one
            self.molecules[index] = src;
            self.molecules.push(dst);
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

        // Try to bond all pairs of molecules
        // TODO: Inefficient, we only really need to check moved (?); especially with all the cloning
        'settling: loop {
            for i in 0..self.molecules.len() {
                for j in 0..self.molecules.len() {
                    if i == j {
                        continue;
                    }

                    let mut primary = self.molecules[i].clone();

                    if primary.try_bond(Point::ZERO, &self.molecules[j]) {
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
    fn test_move_and_bond() {
        use super::*;

        let (map, mut state) = Map::load("H-h#");

        assert!(state.try_move(&map, 0, Point(1, 0)));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].offset, Point(1, 0));
        assert_eq!(state.molecules[0].elements.len(), 2);
    }

    #[test]
    fn test_push_and_bond() {
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

    #[test]
    fn test_move_across_splitter() {
        use super::*;

        let (map, mut state) = Map::load("H\\--\nh---");

        // Map will look like:
        // H - - -
        // |  +
        // H - - -
        
        // Move once, still together
        assert!(state.try_move(&map, 0, Point(1, 0)));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 2);
        
        // Move again, now we're split
        assert!(state.try_move(&map, 0, Point(1, 0)));
        assert_eq!(state.molecules.len(), 2);
        assert_eq!(state.molecules[0].elements.len(), 1);
        assert_eq!(state.molecules[1].elements.len(), 1);

        // Molecule 0 is the one that moved
        assert_eq!(state.molecules[0].offset, Point(2, 0));
        assert_eq!(state.molecules[1].offset, Point(1, 1));

        // Move once more, doesn't rejoin or anything
        assert!(state.try_move(&map, 0, Point(1, 0)));
        assert_eq!(state.molecules.len(), 2);
        assert_eq!(state.molecules[0].offset, Point(3, 0));
        assert_eq!(state.molecules[1].offset, Point(1, 1));
    }

    #[test]
    fn test_partial_split() {
        // Map will look like:
        // O-H - -
        // |  +
        // H - - -

        use super::*;

        let (map, mut state) = Map::load("Oh\\--\nh----");

        // Move twice, still together
        assert!(state.try_move(&map, 0, Point(1, 0)));
        assert!(state.try_move(&map, 0, Point(1, 0)));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 3);

        // Move again, now we're split
        assert!(state.try_move(&map, 0, Point(1, 0)));
        assert_eq!(state.molecules.len(), 2);
        assert_eq!(state.molecules[0].elements.len(), 2);
        assert_eq!(state.molecules[1].elements.len(), 1);

        // Molecule 0 is the one that moved
        assert_eq!(state.molecules[0].offset, Point(3, 0));
        assert_eq!(state.molecules[1].offset, Point(2, 1));

        // Move down, rejoin
        assert!(state.try_move(&map, 0, Point(0, 1)));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 3);
        assert_eq!(state.molecules[0].offset, Point(3, 1));


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
    fn is_valid(&self, map: &Map) -> bool {
        if self.is_solved(map) {
            return true;
        }

        // If we have any splitters, anything goes :)
        if !map.splitters.is_empty() {
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
        let mut grid = vec![vec![' '; map.width * 2]; map.height * 2];

        for y in 0..map.height {
            for x in 0..map.width {
                if map.walls[y * map.width + x] {
                    grid[y * 2][x * 2] = '#';
                } else {
                    grid[y * 2][x * 2] = '-';
                }
            }
        }

        for (i, molecule) in self.molecules.iter().enumerate() {
            // Add the elements
            for (element, offset, _) in &molecule.elements {
                let offset = molecule.offset + *offset;
                let mut c: char = element.into();
                if i > 0 {
                    c = c.to_ascii_lowercase();
                }
                grid[2 * offset.1 as usize][2 * offset.0 as usize] = c;
            }

            // Add the bonds
            for (src, dst, count) in &molecule.bonds {
                let src = molecule.offset + *src;
                let dst = molecule.offset + *dst;
                assert!(src.manhattan_distance(dst) == 1);

                // Offset from src to dst
                let dx = dst.0 - src.0;
                let dy = dst.1 - src.1;

                // Use abs to choose horizontal or vertical; count to choose single, double, or triple
                let c = match (dx.abs(), dy.abs(), count) {
                    (1, 0, 1) => SINGLE_HORIZONTAL,
                    (0, 1, 1) => SINGLE_VERTICAL,
                    (1, 0, 2) => DOUBLE_HORIZONTAL,
                    (0, 1, 2) => DOUBLE_VERTICAL,
                    (1, 0, 3) => TRIPLE_HORIZONTAL,
                    (0, 1, 3) => TRIPLE_VERTICAL,
                    _ => panic!("invalid bond in {molecule:?}: {src:?} {dst:?}"),
                };

                // Place the bond on the offset grid
                grid[(dy + 2 * src.1) as usize][(dx + 2 * src.0) as usize] = c;
            }
        }

        // Add splitters
        for splitter in &map.splitters {
            grid[1 + 2 * splitter.1 as usize][1 + 2 * splitter.0 as usize] = '+';
        }

        let mut output = String::new();
        for y in 0..(2*map.height) {
            for x in 0..(2*map.width) {
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
        // if solver.states_checked() % 10000 != 0 {
        //     continue;
        // }

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
    // test!{test_03_05, "03 - Gray", "05 - Pathways.txt", ""}

    test!{test_04_01, "04 - Red", "01 - Split.txt", "DDWWDSSDDWWWAASSDWDSSAAWDSDDAAWWWASSWWDSAWDDS"}
}