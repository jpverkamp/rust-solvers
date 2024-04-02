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
#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
struct Point {
    x: isize,
    y: isize,
}

impl Add<Point> for Point {
    type Output = Point;

    fn add(self, rhs: Point) -> Self::Output {
        Point {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Sub<Point> for Point {
    type Output = Point;

    fn sub(self, rhs: Point) -> Self::Output {
        Point {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl Point {
    const ZERO: Point = Point { x: 0, y: 0 };

    fn manhattan_distance(&self, other: Point) -> isize {
        (self.x - other.x).abs() + (self.y - other.y).abs()
    }
}

// Possible kinds of map modifiers
#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
enum ModifierKind {
    Weaken,
    Strengthen,
    Rotate,
}

impl From<char> for ModifierKind {
    fn from(value: char) -> Self {
        match value {
            '/' => ModifierKind::Weaken,
            '+' => ModifierKind::Strengthen,
            '@' => ModifierKind::Rotate,
            _ => panic!("unknown modifier: {}", value),
        }
    }
}

impl Into<char> for ModifierKind {
    fn into(self) -> char {
        match self {
            ModifierKind::Weaken => '/',
            ModifierKind::Strengthen => '+',
            ModifierKind::Rotate => '@',
        }
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
struct Modifier {
    kind: ModifierKind,
    location: Point,
}

// Global state: a set of walls
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct Map {
    width: usize,
    height: usize,
    walls: Vec<bool>,
    modifiers: Vec<Modifier>,
}

impl Map {
    fn load(input: &str) -> (Map, LocalState) {
        let mut width = 0;
        let mut height = 0;

        let mut walls = Vec::new();
        let mut modifiers = Vec::new();

        let mut molecules = Vec::new();

        // Version 1 files are just the grid of elements

        // Version 2 files start with v2 on the first line
        // After that, grid elements are spaced out with extra items in between, like:
        // H - -
        //    /
        // H - -

        let mut lines = input.lines().peekable();

        // In v2, build a new input array as alternating lines
        let grid = if lines.peek().is_some_and(|l| l.starts_with("v2")) {
            lines.next(); // Skip the v2 line
            let mut new_input = String::new();

            let mut y = 0;
            while let Some(line) = lines.next() {
                y += 1;
                if y % 2 == 1 {
                    for (x, c) in line.chars().enumerate() {
                        if x % 2 == 0 {
                            new_input.push(c);
                        } else if c != ' ' {
                            panic!("unexpected character in v2 grid spacing: {}", c);
                        }
                    }
                    new_input.push('\n');
                } else {
                    for (x, c) in line.chars().enumerate() {
                        if c == ' ' {
                            continue;
                        }

                        let x = x / 2;
                        let y = y / 2 - 1;
                        modifiers.push(Modifier {
                            location: Point {
                                x: x as isize,
                                y: y as isize,
                            },
                            kind: ModifierKind::from(c),
                        });
                    }
                }
            }

            new_input
        } else {
            String::from(input)
        };

        // Now grid contains just the elements, walls, and empty space
        for (y, line) in grid.lines().enumerate() {
            width = width.max(line.len());
            height += 1;

            for (x, c) in line.chars().enumerate() {
                let pt = Point {
                    x: x as isize,
                    y: y as isize,
                };

                match ElementKind::try_from(c) {
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
                    Err(_) => match c {
                        ' ' | '-' => walls.push(false),
                        'x' | 'X' | '#' => walls.push(true),
                        _ => panic!("unknown character: {}", c),
                    },
                }
            }
        }

        // Try to bond each original pair of molecules
        'settled: loop {
            for i in 0..molecules.len() {
                for j in (i + 1)..molecules.len() {
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

        // Modifiers have to be applied in a specific order: weaken, strengthen, rotate(?)
        modifiers.sort_by_key(|m| m.kind);

        (
            Map {
                width,
                height,
                walls,
                modifiers,
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

#[cfg(test)]
mod test_map {
    #[test]
    fn test_load() {
        use super::*;

        let (map, state) = Map::load("H-#");

        assert_eq!(map.width, 3);
        assert_eq!(map.height, 1);
        assert!(map.is_wall(2, 0));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 1);
    }

    #[test]
    fn test_load_prebond() {
        use super::*;

        let (_, state) = Map::load("Hh-#");

        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 2);
        assert_eq!(state.molecules[0].bonds.len(), 1);
        assert_eq!(state.molecules[0].bonds[0].count, 1);
    }
}

// Represents single elements
#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
enum ElementKind {
    Hydrogen,
    Helium,
    Nitrogen,
    Carbon,
    Oxygen,
}

impl TryFrom<char> for ElementKind {
    type Error = ();

    fn try_from(value: char) -> Result<Self, Self::Error> {
        use ElementKind::*;

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

impl Into<char> for ElementKind {
    fn into(self) -> char {
        match self {
            ElementKind::Hydrogen => 'H',
            ElementKind::Helium => 'E',
            ElementKind::Nitrogen => 'N',
            ElementKind::Carbon => 'C',
            ElementKind::Oxygen => 'O',
        }
    }
}

impl ElementKind {
    fn free_electrons(&self) -> usize {
        use ElementKind::*;

        match self {
            Hydrogen => 1,
            Helium => 0,
            Nitrogen => 3,
            Carbon => 4,
            Oxygen => 2,
        }
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
struct Element {
    kind: ElementKind,
    offset: Point,
    free_electrons: usize,
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
struct Bond {
    a: Point,
    b: Point,
    count: usize,
}

// Represent elements joined together into a molecule
// Each element contains the Element, offset, and remaining free electrons
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct Molecule {
    offset: Point,
    elements: Vec<Element>,
    bonds: Vec<Bond>,
}

impl Molecule {
    fn new(offset: Point, element: ElementKind) -> Molecule {
        Molecule {
            offset,
            elements: vec![Element {
                kind: element,
                offset: Point::ZERO,
                free_electrons: element.free_electrons(),
            }],
            bonds: Vec::new(),
        }
    }

    // Test for molecular helium
    fn is_helium(&self) -> bool {
        self.elements.len() == 1 && self.elements[0].kind == ElementKind::Helium
    }

    // How many free electrons in the hole molecule
    fn free_electrons(&self) -> usize {
        self.elements
            .iter()
            .map(|e| e.free_electrons)
            .sum::<usize>()
    }

    // If the given molecule at an offset would intersect with a wall
    fn intersects_wall(&self, offset: Point, map: &Map) -> bool {
        for element in &self.elements {
            let target = self.offset + element.offset + offset;
            if map.is_wall(target.x as usize, target.y as usize) {
                return true;
            }
        }

        false
    }

    // If the given molecule + an offset would intersect another molecule
    fn intersects(&self, offset: Point, other: &Molecule) -> bool {
        for element in &self.elements {
            for other_element in &other.elements {
                if self.offset + element.offset + offset == other.offset + other_element.offset {
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
        for a in self.elements.iter_mut() {
            // Skip our elements that are no longer free
            if a.free_electrons == 0 {
                continue;
            }

            for b in other.elements.iter_mut() {
                let real_a = self.offset + offset + a.offset;
                let real_b = other.offset + b.offset;

                // Not adjacent
                if real_a.manhattan_distance(real_b) != 1 {
                    continue;
                }

                // Not enough free electrons
                if b.free_electrons == 0 {
                    continue;
                }

                // Bond the two elements
                bound = true;

                self.bonds.push(Bond {
                    a: offset + a.offset,
                    b: other.offset - self.offset + b.offset,
                    count: 1,
                });

                a.free_electrons -= 1;
                b.free_electrons -= 1;
            }
        }

        // If we bound anything, add the other elements and bonds to our molecule
        if bound {
            for element in other.elements {
                self.elements.push(Element {
                    kind: element.kind,
                    offset: other.offset - self.offset + element.offset,
                    free_electrons: element.free_electrons,
                });
            }

            for bond in other.bonds {
                self.bonds.push(Bond {
                    a: other.offset - self.offset + bond.a,
                    b: other.offset - self.offset + bond.b,
                    count: bond.count,
                });
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
        let a = Molecule::new(Point::ZERO, ElementKind::Hydrogen);

        let map = Map {
            width: 3,
            height: 3,
            walls: vec![
                true, true, true, true, false, true, true, true, true, // 3x3
            ],
            modifiers: vec![],
        };

        assert!(a.intersects_wall(Point { x: 1, y: 0 }, &map));
        assert!(!a.intersects_wall(Point { x: 1, y: 1 }, &map));
    }

    #[test]
    fn test_molecule_intersection() {
        let a = Molecule::new(Point::ZERO, ElementKind::Hydrogen);
        let b = Molecule::new(Point { x: 1, y: 0 }, ElementKind::Hydrogen);

        assert!(a.intersects(Point { x: 1, y: 0 }, &b));
        assert!(!a.intersects(Point { x: 0, y: 0 }, &b));
        assert!(!a.intersects(Point { x: 0, y: 1 }, &b));
    }

    #[test]
    fn test_bond() {
        let mut a = Molecule::new(Point::ZERO, ElementKind::Hydrogen);
        let b = Molecule::new(Point { x: 1, y: 0 }, ElementKind::Hydrogen);

        let bound = a.try_bond(Point::ZERO, &b);

        assert!(bound);
        assert_eq!(a.elements.len(), 2);
        assert_eq!(a.elements[0].free_electrons, 0);
    }

    #[test]
    fn test_nobond_no_free() {
        let mut a = Molecule::new(Point::ZERO, ElementKind::Hydrogen);
        let b = Molecule::new(Point { x: 1, y: 0 }, ElementKind::Helium);

        let bound = a.try_bond(Point { x: 1, y: 0 }, &b);

        assert!(!bound);
        assert!(a.elements[0].free_electrons == 1);
    }

    #[test]
    fn test_nobond_too_far() {
        let mut a = Molecule::new(Point::ZERO, ElementKind::Hydrogen);
        let b = Molecule::new(Point { x: 2, y: 0 }, ElementKind::Hydrogen);

        let bound = a.try_bond(Point { x: 2, y: 0 }, &b);

        assert!(!bound);
        assert!(a.elements[0].free_electrons == 1);
    }

    #[test]
    fn test_single_bond_o2_not_double() {
        let mut a = Molecule::new(Point::ZERO, ElementKind::Oxygen);
        let b = Molecule::new(Point { x: 1, y: 0 }, ElementKind::Oxygen);

        let bound = a.try_bond(Point::ZERO, &b);

        assert!(bound);
        assert_eq!(a.elements.len(), 2);
        assert_eq!(a.elements[0].free_electrons, 1);
        assert_eq!(a.elements[1].free_electrons, 1);
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

        // Collect all map modifiers that we are trying to cross (this may take multiple passes)
        let mut modifiers_applied = Vec::new();
        loop {
            let mut bond_to_modify = None;
            
            'find_bond: for modifier in map.modifiers.iter() {
                if modifiers_applied.contains(&modifier) {
                    continue;
                }

                for (bond_index, bond) in self.molecules[index].bonds.iter().enumerate() {
                    let real_a = bond.a + self.molecules[index].offset;
                    let real_b = bond.b + self.molecules[index].offset;

                    // Vertical bonds have the same x
                    let is_vertical = bond.a.x == bond.b.x;

                    // We'll hit a vertical splitter if the offset is horizontal and we're moving across it
                    // Ignore bonds that are moving the wrong way
                    if is_vertical && offset.x == 0 || !is_vertical && offset.y == 0 {
                        continue;
                    }

                    // Moving 'positive' is down or right
                    let is_positive = offset.x > 0 || offset.y > 0;

                    // Because either x or y is the same for all bonds, min is top/left and max is bottom/right
                    // This will always match the splitter if we're moving across it right or down
                    let pre_min = Point {
                        x: real_a.x.min(real_b.x),
                        y: real_a.y.min(real_b.y),
                    };

                    // The post point is the one after we've moved
                    let post_a = real_a + offset;
                    let post_b = real_b + offset;
                    let post_min = Point {
                        x: post_a.x.min(post_b.x),
                        y: post_a.y.min(post_b.y),
                    };

                    // If we're moving positive, the min (top left) will equal the splitter
                    if is_positive && modifier.location != pre_min {
                        continue;
                    }

                    // If we're moving negative, then the *post* min will equal the splitter
                    if !is_positive && modifier.location != post_min {
                        continue;
                    }

                    // We have a bond to try to modify
                    bond_to_modify = Some((bond_index, modifier));
                    break 'find_bond;
                }
            }

            // We found no more bonds to modify
            if bond_to_modify.is_none() {
                break;
            }

            // Note it, so we don't apply the same modifier more than once per move
            let (bond_index, modifier) = bond_to_modify.unwrap();
            modifiers_applied.push(modifier);

            // Figure out which elements we're dealing with
            let bond = self.molecules[index].bonds[bond_index];

            let el_a_index = self.molecules[index]
                .elements
                .iter()
                .position(|el| el.offset == bond.a)
                .unwrap();

            let el_b_index = self.molecules[index]
                .elements
                .iter()
                .position(|el| el.offset == bond.b)
                .unwrap();

            // Handle different modifier types
            match modifier.kind {
                ModifierKind::Weaken => {
                    // Reduce the bond and give back electrons
                    self.molecules[index].elements[el_a_index].free_electrons += 1;
                    self.molecules[index].elements[el_b_index].free_electrons += 1;
                    self.molecules[index].bonds[bond_index].count -= 1;

                    // If it was more than a single bond (originally), we're done now (no splitting)
                    if self.molecules[index].bonds[bond_index].count > 0 {
                        continue;
                    }

                    // We're going to modify the molecule, clone it and we'll put the new one in later
                    let mut src = self.molecules[index].clone();
                    src.bonds.remove(bond_index);

                    // Now starting at the origin in src, remove anything we can get to in dst
                    // TODO: This can be factored into a function
                    let mut connected_elements = Vec::new();
                    let mut connected_bonds = Vec::new();

                    let mut todo = vec![Point::ZERO];
                    let mut done = vec![];

                    while let Some(pt) = todo.pop() {
                        done.push(pt);

                        for (i, element) in src.elements.iter().enumerate() {
                            if element.offset == pt && !connected_elements.contains(&i) {
                                connected_elements.push(i);
                            }
                        }

                        for (i, src_bond) in src.bonds.iter().enumerate() {
                            if src_bond.a == pt {
                                if !connected_bonds.contains(&i) {
                                    connected_bonds.push(i);
                                }

                                if !(todo.contains(&src_bond.b) || done.contains(&src_bond.b)) {
                                    todo.push(src_bond.b);
                                }
                            }

                            if src_bond.b == pt {
                                if !connected_bonds.contains(&i) {
                                    connected_bonds.push(i);
                                }

                                if !(todo.contains(&src_bond.a) || done.contains(&src_bond.a)) {
                                    todo.push(src_bond.a);
                                }
                            }
                        }
                    }

                    // If we actually keep all of the elements, we don't need to modify src/dst at all
                    if connected_elements.len() == src.elements.len() {
                        self.molecules[index] = src;
                        continue;
                    }

                    // Otherwise, we're going to create a second molecule
                    // This is basically a partition (src gets connected, dst gets everything else)
                    let mut dst = src.clone();
                    connected_elements.sort();
                    for i in connected_elements.iter().rev() {
                        dst.elements.remove(*i);
                    }
                    connected_bonds.sort();
                    for i in connected_bonds.iter().rev() {
                        dst.bonds.remove(*i);
                    }

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

                    // Now update the offset in dst so that one of the elements is at 0,0
                    // Not strictly necessary, but I think it will make rotation easier later?
                    let new_zero = dst.elements[0].offset;
                    dst.offset = dst.offset + new_zero;

                    for element in &mut dst.elements {
                        element.offset = element.offset - new_zero;
                    }
                    for dst_bond in &mut dst.bonds {
                        dst_bond.a = dst_bond.a - new_zero;
                        dst_bond.b = dst_bond.b - new_zero;
                    }

                    // We now have two molecules, replace src and add dst
                    self.molecules[index] = src;
                    self.molecules.push(dst);
                }
                ModifierKind::Strengthen => {
                    // Verify we have enough free electrons
                    if self.molecules[index].elements[el_a_index].free_electrons == 0
                        || self.molecules[index].elements[el_b_index].free_electrons == 0
                    {
                        continue;
                    }

                    // If so, strengthen the bond and take the electrons
                    self.molecules[index].elements[el_a_index].free_electrons -= 1;
                    self.molecules[index].elements[el_b_index].free_electrons -= 1;
                    self.molecules[index].bonds[bond_index].count += 1;
                }
                ModifierKind::Rotate => todo!(),
            }
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

        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules[0].offset, Point { x: 1, y: 0 });
    }

    #[test]
    fn test_bump_into_wall() {
        use super::*;

        let (map, mut state) = Map::load("-H#");

        assert!(!state.try_move(&map, 0, Point { x: 1, y: 0 }));
    }

    #[test]
    fn test_move_push() {
        use super::*;

        let (map, mut state) = Map::load("Ee-#");

        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules.len(), 2);
        assert_eq!(state.molecules[0].offset, Point { x: 1, y: 0 });
        assert_eq!(state.molecules[1].offset, Point { x: 2, y: 0 });
    }

    #[test]
    fn test_move_and_bond() {
        use super::*;

        let (map, mut state) = Map::load("H-h#");

        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].offset, Point { x: 1, y: 0 });
        assert_eq!(state.molecules[0].elements.len(), 2);
    }

    #[test]
    fn test_push_and_bond() {
        use super::*;

        let (map, mut state) = Map::load("Hh-#");

        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].offset, Point { x: 1, y: 0 });
        assert_eq!(state.molecules[0].elements.len(), 2);
    }

    #[test]
    fn test_no_move_push_into_wall() {
        use super::*;

        let (map, mut state) = Map::load("Ee#");

        assert!(!state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules[0].offset, Point { x: 0, y: 0 });
        assert_eq!(state.molecules[1].offset, Point { x: 1, y: 0 });
    }

    #[test]
    fn test_push_unbound() {
        use super::*;

        let (map, mut state) = Map::load("-H-\no--\n-h-\n---");

        assert!(state.try_move(&map, 0, Point { x: 0, y: 1 }));
        assert!(state.try_move(&map, 0, Point { x: 0, y: 1 }));

        println!("{}", state.stringify(&map));

        // Two molecules
        assert_eq!(state.molecules.len(), 2);

        // First is OH with a free still on the O
        assert_eq!(state.molecules[0].offset, Point { x: 1, y: 2 });
        assert_eq!(state.molecules[0].elements[0].free_electrons, 0);
        assert_eq!(state.molecules[0].elements[1].free_electrons, 1);

        // Second is h with a free still open
        assert_eq!(state.molecules[1].offset, Point { x: 1, y: 3 });
        assert_eq!(state.molecules[1].elements[0].free_electrons, 1);
    }

    #[test]
    fn test_move_across_splitter() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
H - - -
   /
h - - -",
        );

        // Move once, still together
        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 2);

        // Move again, now we're split
        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        println!("{}", state.stringify(&map));
        assert_eq!(state.molecules.len(), 2);
        assert_eq!(state.molecules[0].elements.len(), 1);
        assert_eq!(state.molecules[1].elements.len(), 1);

        // Molecule 0 is the one that moved
        assert_eq!(state.molecules[0].offset, Point { x: 2, y: 0 });
        assert_eq!(state.molecules[1].offset, Point { x: 1, y: 1 });

        // Move once more, doesn't rejoin or anything
        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules.len(), 2);
        assert_eq!(state.molecules[0].offset, Point { x: 3, y: 0 });
        assert_eq!(state.molecules[1].offset, Point { x: 1, y: 1 });
    }

    #[test]
    fn test_partial_split() {
        // Map will look like:
        // O-H - -
        // |  +
        // H - - -

        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
O h - - -
     /
h - - - -",
        );

        // Move twice, still together
        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 3);

        // Move again, now we're split
        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules.len(), 2);
        assert_eq!(state.molecules[0].elements.len(), 2);
        assert_eq!(state.molecules[1].elements.len(), 1);

        // Molecule 0 is the one that moved
        assert_eq!(state.molecules[0].offset, Point { x: 3, y: 0 });
        assert_eq!(state.molecules[1].offset, Point { x: 2, y: 1 });

        // Move down, rejoin
        assert!(state.try_move(&map, 0, Point { x: 0, y: 1 }));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 3);
        assert_eq!(state.molecules[0].offset, Point { x: 3, y: 1 });
    }

    #[test]
    fn test_small_key() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
x - h - x 
   / /
x - N - - 
   / /
x - h - x 
",
        );

        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));

        // We should have moved just the N split from the two h
        assert_eq!(state.molecules.len(), 3);
        assert_eq!(state.molecules[0].elements.len(), 1);
        assert_eq!(state.molecules[0].offset, Point { x: 3, y: 1 });
        assert_eq!(state.molecules[1].elements.len(), 1);
        assert_eq!(state.molecules[2].elements.len(), 1);
    }

    #[test]
    fn test_doubler() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
O - - -
   + 
o - - -",
        );

        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));

        // We should have moved the O and o together with a single bond
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 2);
        assert_eq!(state.molecules[0].elements[0].free_electrons, 1);
        assert_eq!(state.molecules[0].elements[1].free_electrons, 1);
        assert_eq!(state.molecules[0].bonds.len(), 1);
        assert_eq!(state.molecules[0].bonds[0].count, 1);

        // After second move, we should have a double bond and no more free electrons
        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 2);
        assert_eq!(state.molecules[0].elements[0].free_electrons, 0);
        assert_eq!(state.molecules[0].elements[1].free_electrons, 0);
        assert_eq!(state.molecules[0].bonds.len(), 1);
        assert_eq!(state.molecules[0].bonds[0].count, 2);
    }

    #[test]
    fn test_cats_cradle() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
- - c
   + 
n - N",
        );

        assert!(state.try_move(&map, 0, Point { x: -1, y: 0 }));

        // Should result in one molecule with the right C/N double bonded and a single N/N bond
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 3);
        assert_eq!(state.molecules[0].elements[0].free_electrons, 0); // N has 0/3 free left
        assert_eq!(state.molecules[0].elements[1].free_electrons, 2); // c has 2/4 free left
        assert_eq!(state.molecules[0].elements[2].free_electrons, 2); // Other n has 2/3 free left
        assert_eq!(state.molecules[0].bonds.len(), 2);
        assert_eq!(state.molecules[0].bonds[0].count, 2); // C/N double bond
        assert_eq!(state.molecules[0].bonds[1].count, 1); // N/N single bond
    }

    #[test]
    fn test_split_multiple() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
C - - -
 + + /
c - - -",
        );

        // First move double bond
        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules[0].bonds.len(), 1);
        assert_eq!(state.molecules[0].bonds[0].count, 2);

        // Second move triple bond
        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules[0].bonds.len(), 1);
        assert_eq!(state.molecules[0].bonds[0].count, 3);

        // Third move makes it a double again
        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
        assert_eq!(state.molecules[0].bonds.len(), 1);
        assert_eq!(state.molecules[0].bonds[0].count, 2);

        // Going back will make it a single
        assert!(state.try_move(&map, 0, Point { x: -1, y: 0 }));
        assert_eq!(state.molecules[0].bonds.len(), 1);
        assert_eq!(state.molecules[0].bonds[0].count, 1);
    }

    #[test]
    fn test_splitter_second_join() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
- -
 /
o o
   
N n",
        );

        // Move one up, should break the o-o bond, but keep one molecule
        assert!(state.try_move(&map, 0, Point { x: 0, y: -1 }));
        assert_eq!(state.molecules.len(), 1);
    }

    #[test]
    fn test_split_join() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
C n n
 + /
- - -");
        
        // Verify that we have one molecule with two single bonds
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].bonds.len(), 2);
        assert_eq!(state.molecules[0].bonds[0].count, 1);
        assert_eq!(state.molecules[0].bonds[1].count, 1);

        // Move down one, should double the C/N bond and split (plus leave) the other N
        assert!(state.try_move(&map, 0, Point { x: 0, y: 1 }));
        assert_eq!(state.molecules.len(), 2);
        assert_eq!(state.molecules[0].offset, Point { x: 0, y: 1 });
        assert_eq!(state.molecules[0].bonds.len(), 1);
        assert_eq!(state.molecules[0].bonds[0].count, 2);
        assert_eq!(state.molecules[1].offset, Point { x: 2, y: 0 });
        assert_eq!(state.molecules[1].bonds.len(), 0);
    }
    
    #[test]
    fn test_join_split() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
C n n
 / +
- - -");

        // Verify that we have one molecule with two single bonds
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].bonds.len(), 2);
        assert_eq!(state.molecules[0].bonds[0].count, 1);
        assert_eq!(state.molecules[0].bonds[1].count, 1);

        // Move down one, only the C should move leaving the N/N where it was
        assert!(state.try_move(&map, 0, Point { x: 0, y: 1 }));
        assert_eq!(state.molecules.len(), 2);
        assert_eq!(state.molecules[0].offset, Point { x: 0, y: 1 });
        assert_eq!(state.molecules[0].elements.len(), 1);
        assert_eq!(state.molecules[0].bonds.len(), 0);
        assert_eq!(state.molecules[1].offset, Point { x: 1, y: 0 });
        assert_eq!(state.molecules[1].bonds.len(), 1);
        assert_eq!(state.molecules[1].bonds[0].count, 1);
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
            Step::North => Point { x: 0, y: -1 },
            Step::South => Point { x: 0, y: 1 },
            Step::East => Point { x: 1, y: 0 },
            Step::West => Point { x: -1, y: 0 },
        }
    }
}

impl State<Map, Step> for LocalState {
    fn is_valid(&self, map: &Map) -> bool {
        if self.is_solved(map) {
            return true;
        }

        // If we have any splitters, anything goes :)
        if !map.modifiers.is_empty() {
            return true;
        }

        // The primary molecule must have free electrons
        if !(self.molecules[0].is_helium() || self.molecules[0].free_electrons() > 0) {
            return false;
        }

        // Any non-primary (other than helium) must have free electrons
        if self
            .molecules
            .iter()
            .skip(1)
            .any(|m| !m.is_helium() && m.free_electrons() == 0)
        {
            return false;
        }

        return true;
    }

    fn is_solved(&self, _global: &Map) -> bool {
        self.molecules
            .iter()
            .all(|m| m.is_helium() || m.free_electrons() == 0)
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
        self.molecules
            .iter()
            .map(|m| m.free_electrons() as i64)
            .sum()
    }

    fn stringify(&self, map: &Map) -> String {
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
            for element in &molecule.elements {
                let offset = molecule.offset + element.offset;
                let mut c: char = element.kind.into();
                if i > 0 {
                    c = c.to_ascii_lowercase();
                }
                grid[2 * offset.y as usize][2 * offset.x as usize] = c;
            }

            // Add the bonds
            for bond in &molecule.bonds {
                let real_a = molecule.offset + bond.a;
                let real_b = molecule.offset + bond.b;
                assert!(real_a.manhattan_distance(real_b) == 1);

                // Offset from src to dst
                let dx = real_b.x - real_a.x;
                let dy = real_b.y - real_a.y;

                // Use abs to choose horizontal or vertical; count to choose single, double, or triple
                let c = match (dx.abs(), dy.abs(), bond.count) {
                    (1, 0, 1) => SINGLE_HORIZONTAL,
                    (0, 1, 1) => SINGLE_VERTICAL,
                    (1, 0, 2) => DOUBLE_HORIZONTAL,
                    (0, 1, 2) => DOUBLE_VERTICAL,
                    (1, 0, 3) => TRIPLE_HORIZONTAL,
                    (0, 1, 3) => TRIPLE_VERTICAL,
                    _ => panic!("invalid bond in {molecule:?}: {real_a:?} {real_b:?}"),
                };

                // Place the bond on the offset grid
                grid[(dy + 2 * real_a.y) as usize][(dx + 2 * real_a.x) as usize] = c;
            }
        }

        // Add splitters
        for modifier in &map.modifiers {
            let c: char = modifier.kind.into();
            grid[1 + 2 * modifier.location.y as usize][1 + 2 * modifier.location.x as usize] = c;
        }

        let mut output = String::new();
        for y in 0..(2 * map.height) {
            for x in 0..(2 * map.width) {
                output.push(grid[y][x]);
            }
            output.push('\n');
        }
        output
    }
}

fn solve(input: &str) -> Option<String> {
    let (map, molecules) = Map::load(input);

    log::info!("Initial state:\n{}", molecules.stringify(&map));

    let mut solver = Solver::new(map.clone(), molecules.clone());

    while let Some(state) = solver.next() {
        if solver.states_checked() % 100000 != 0 {
            continue;
        }
        log::info!("{solver}, state:\n{}", state.stringify(&map));
    }

    let solution = solver.get_solution();
    if solution.is_none() {
        return None;
    }
    let solution = solution.unwrap();

    log::info!(
        "Solved after {} states in {} seconds:\n{}",
        solver.states_checked(),
        solver.time_spent(),
        solution.stringify(&map),
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

    test! {test_01_01, "01 - Yellow", "01 - Let's Go.txt", "WWDDWWDD"}
    test! {test_01_02, "01 - Yellow", "02 - Cell.txt", "DWDDAA"}
    test! {test_01_03, "01 - Yellow", "03 - Loop.txt", "SSDDWWSSAAAWWWWD"}
    test! {test_01_04, "01 - Yellow", "04 - Monty Hall.txt", "WSDDWSAAWWWDAAA"}
    test! {test_01_05, "01 - Yellow", "05 - Knot.txt", "WDASASWADDDD"}
    test! {test_01_06, "01 - Yellow", "06 - Push.txt", "DDSDWD"}
    test! {test_01_07, "01 - Yellow", "07 - Structure.txt", "DSSAWAA"}
    test! {test_01_08, "01 - Yellow", "08 - Lotus.txt", "SASSDDSAS"}
    test! {test_01_09, "01 - Yellow", "09 - Coyote.txt", "WWAASASSAA"}
    test! {test_01_10, "01 - Yellow", "10 - Roadrunner.txt", "WDWDSDDWAAASSD"}

    test! {test_02_01, "02 - Orange", "01 - Suit.txt", "SSDDWSSDSAS"}
    test! {test_02_02, "02 - Orange", "02 - Chimney.txt", "WDAWSSAAWDDW"}
    test! {test_02_03, "02 - Orange", "03 - Heart.txt", "WWDWASAAWDAWWDSDDW"}
    test! {test_02_04, "02 - Orange", "04 - Factory.txt", "DDSSAAAAAWWWSDDDDDWWAWAS"}
    test! {test_02_05, "02 - Orange", "05 - Scoop.txt", "DAASSDSDDAWAAA"}
    test! {test_02_06, "02 - Orange", "06 - Block.txt", "ASSDDSDDWWWADSSAAAW"}
    test! {test_02_07, "02 - Orange", "07 - Chandelier.txt", "DSSSAWDWDSSAWWWA"}
    test! {test_02_08, "02 - Orange", "08 - Creature.txt", "WDDWWAAWASDDDDSSAASSAWAAW"}
    test! {test_02_09, "02 - Orange", "09 - Rosie.txt", "WAWWDDAASASSDDDD"}
    test! {test_02_10, "02 - Orange", "10 - Kruskal.txt", "ASAWDWDDSSSWAAWW"}

    test! {test_03_01, "03 - Gray", "01 - Helium.txt", "WDDDDDSAAWASSDSSS"}
    test! {test_03_02, "03 - Gray", "02 - Tee.txt", "WWASWDDDSAAWASDSSADSSWWW"}
    test! {test_03_03, "03 - Gray", "03 - Freedom.txt", "AAWAWDD"}
    test! {test_03_04, "03 - Gray", "04 - Against the Wall.txt", "WAASWAWWDSASDDSSSAWWSAAWWDDWWDDSSAAWAWASASDDDD"}
    test! {test_03_05, "03 - Gray", "05 - Pathways.txt", "AWWDSASDSSSDDWWAASSAAWW"}
    test! {test_03_06, "03 - Gray", "06 - Three Doors.txt", "DAWWSAAWWAWDDSDSSAASDSDWWWSDDWSDSSAWWAAAAWWDD"}
    test! {test_03_07, "03 - Gray", "07 - Cloud.txt", "DWWASDDSWWASSD"}
    test! {test_03_08, "03 - Gray", "08 - Planning.txt", "WDSDSAASSAWDWASAAWDDDDSDDW"}
    test! {test_03_09, "03 - Gray", "09 - Out of the Way.txt", "AWDDDSAWAAAWWDDDDSWAAAASSDDWSDSDDWAAADDWWAAADSSASAW"}
    test! {test_03_10, "03 - Gray", "10 - Impasse.txt", "DWADDSASAAWADDSSWWDWWA"}
    test! {test_03_11, "03 - Gray", "11 - Fetch.txt", "WDWASSSADDASWWWWDSAWASDSDAA"}
    test! {test_03_12, "03 - Gray", "12 - Drill.txt", "AWWWWDSASSSDWDWAWWAASDSDASSDWAWWDSWWDDSAS"}

    test! {test_04_01, "04 - Red", "01 - Split.txt", "DDWWDSSDDDAAWA"}
    test! {test_04_02, "04 - Red", "02 - Lock.txt", "DDWDSAWWWWAADDSSSA"}
    test! {test_04_03, "04 - Red", "03 - Push Up.txt", "DAAWDSDWA"}
    test! {test_04_04, "04 - Red", "04 - Out of Reach.txt", "WDDDDSAWDSAAWW"}
    test! {test_04_05, "04 - Red", "05 - Small Key.txt", "ASDWDWASSAWSDWDDD"}
    test! {test_04_06, "04 - Red", "06 - Anxiety.txt", "WAASAWDWWSDSDW"}
    test! {test_04_07, "04 - Red", "07 - Wingman.txt", "WAWWSDDDSASDSAAAAWW"}
    test! {test_04_08, "04 - Red", "08 - Wing.txt", "DDAAWWSDWWDAA"}
    test! {test_04_09, "04 - Red", "09 - Anvil.txt", "SDSSDSAADWWWAASSDDDWWADSSAAAWSA"}
    test! {test_04_10, "04 - Red", "10 - Cottage.txt", "SWWWWDDSDSAWWAAAASASDWWDDSSSS"}
    test! {test_04_11, "04 - Red", "11 - Hanoi.txt", "SAAWDDDDSAAWAASDDDDWAAAASDDDDWAAAA"}
    test! {test_04_12, "04 - Red", "12 - Trap.txt", "DSDWAASDWWASDSDDDDDWSAWWDSASAAA"}

    test! {test_05_01, "05 - Green", "01 - Breathe.txt", "WWDSDSWAW"}
    test! {test_05_02, "05 - Green", "02 - Doubling.txt", "WAWWASSDW"}
    test! {test_05_03, "05 - Green", "03 - Koan.txt", "WAWAADDSASWDW"}
    test! {test_05_04, "05 - Green", "04 - Together.txt", "SSWDDWWSSAWA"}
    test! {test_05_05, "05 - Green", "05 - Agoraphobia.txt", "SWDWDDSSSA"}
    test! {test_05_06, "05 - Green", "06 - Forethought.txt", "WDWWSSDDDAAASA"}
    test! {test_05_07, "05 - Green", "07 - Apartment.txt", "SDDWWDDDSDWWWWAAAADDDDW"}
    test! {test_05_08, "05 - Green", "08 - Lettuce.txt", "WDDWAWSSA"}
    test! {test_05_09, "05 - Green", "09 - Landing Pad.txt", "SDDDSSWWWWSSAAAAWSS"}
    test! {test_05_10, "05 - Green", "10 - Matchmaker.txt", "DDWWWDASAAASDDDSS"}
    test! {test_05_11, "05 - Green", "11 - Cat's Cradle.txt", "DWWDAAASSWWDDDSSSAAA"}

    test! {test_06_01, "06 - Dark Green", "01 - Papers Please.txt", "DADDDDDDDASWWASAAA"}
    test! {test_06_02, "06 - Dark Green", "02 - Airplane.txt", "SSWWWAASDDDWDSAAWWSWDSSSADWWDDDSAAWASSASDWWAWAAASDDWDSSDSAWWWWDADSSSSS"} // Slow
    test! {test_06_03, "06 - Dark Green", "03 - Mine Field.txt", "SDDDDSSSAWAWAWDSAAWWSDSDSDSAAAAAA"}
    // test! {test_06_04, "06 - Dark Green", "04 - Workshop.txt", ""}
    test! {test_06_05, "06 - Dark Green", "05 - Dissection.txt", "AASSDDAWADWASDDDDD"}
    test! {test_06_06, "06 - Dark Green", "06 - Casket.txt", "DDAASSDDDWAWSAWSD"}
    test! {test_06_07, "06 - Dark Green", "07 - Three Body Problem.txt", "DDAWDSAWWDDASSAAWWS"}
    test! {test_06_08, "06 - Dark Green", "08 - Cat.txt", "WAADDDAADAAWSSDDWD"}
    test! {test_06_09, "06 - Dark Green", "09 - Halves.txt", "ADDSAWAAWDSDD"}
    test! {test_06_10, "06 - Dark Green", "10 - Arrow.txt", "SWWWADWWSWASSSSDWWWAA"}
    test! {test_06_11, "06 - Dark Green", "11 - Hallway.txt", "DDASAWSAAWDDDDSAWDSWAADDAASDWSAAWSWASDDDA"}

    test! {test_07_01, "07 - Dark Red", "01 - Plunge.txt", "SAWDSDWWSAWASSWDWDWASASDDWWAASDDSDDD"}
    test! {test_07_02, "07 - Dark Red", "02 - Compass.txt", "SSDDASAAWADSDS"}
}
