use std::io;
use std::ops::Add;
use std::ops::Sub;

use solver::{Point, Solver, State};

const SINGLE_HORIZONTAL: char = '-';
const SINGLE_VERTICAL: char = '|';
const DOUBLE_HORIZONTAL: char = '=';
const DOUBLE_VERTICAL: char = '‖';
const TRIPLE_HORIZONTAL: char = '≡';
const TRIPLE_VERTICAL: char = '⦀';

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
    allow_multiple: bool,
    solutions: Vec<String>,
}

impl Map {
    fn load(input: &str) -> (Map, LocalState) {
        let mut width = 0;
        let mut height = 0;

        let mut walls = Vec::new();
        let mut modifiers = Vec::new();

        let mut molecules = Vec::new();

        let mut solutions = Vec::new();

        // Version 1 files are just the grid of elements

        // Version 2 files start with v2 on the first line
        // After that, grid elements are spaced out with extra items in between, like:
        // H - -
        //    /
        // H - -

        let mut lines = input.lines().peekable();

        let mut multiple = false;

        // In v2, build a new input array as alternating lines
        let grid = if lines.peek().is_some_and(|l| l.starts_with("v2")) {
            // Look for additional options
            let options = lines.next().unwrap();
            if options.contains("multiple") {
                multiple = true;
            }

            let mut new_input = String::new();

            let mut y = 0;
            while let Some(line) = lines.next() {
                if line.starts_with("=") {
                    new_input.push_str(line);
                    new_input.push('\n');
                    continue;
                }

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
            // Lines starting with = represent solutions
            // They should only be at the end of the file
            if line.starts_with('=') {
                solutions.push(line[1..].to_string());
                continue;
            }

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
        'bonding: loop {
            for i in 0..molecules.len() {
                for j in 0..molecules.len() {
                    if i == j {
                        continue;
                    }

                    let mut primary = molecules[i].clone();
                    if primary.try_bond(Point::ZERO, &molecules[j]) {
                        molecules[i] = primary;
                        molecules[j].active = false;
                        continue 'bonding;
                    }
                }
            }

            break 'bonding;
        }

        // Remove molecules marked inactive
        molecules.retain(|m| m.active);

        // Modifiers have to be applied in a specific order: weaken, strengthen, rotate(?)
        modifiers.sort_by_key(|m| m.kind);

        (
            Map {
                width,
                height,
                walls,
                modifiers,
                allow_multiple: multiple,
                solutions,
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

impl Add<Point> for Bond {
    type Output = Bond;

    fn add(self, rhs: Point) -> Self::Output {
        Bond {
            a: self.a + rhs,
            b: self.b + rhs,
            count: self.count,
        }
    }
}

impl Sub<Point> for Bond {
    type Output = Bond;

    fn sub(self, rhs: Point) -> Self::Output {
        Bond {
            a: self.a - rhs,
            b: self.b - rhs,
            count: self.count,
        }
    }
}

// Represent elements joined together into a molecule
// Each element contains the Element, offset, and remaining free electrons
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct Molecule {
    active: bool,
    offset: Point,
    elements: Vec<Element>,
    bonds: Vec<Bond>,
}

impl Molecule {
    fn new(offset: Point, element: ElementKind) -> Molecule {
        Molecule {
            active: true,
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
        if !self.active {
            return false;
        }

        for element in &self.elements {
            let target = self.offset + element.offset + offset;
            if map.is_wall(target.x as usize, target.y as usize) {
                return true;
            }
        }

        false
    }

    // If the given molecule + an offset would intersect another molecule
    // If there's an intersection, return the intersecting point of each molecule (without offset)
    fn intersects(&self, offset: Point, other: &Molecule) -> Option<(Point, Point)> {
        if !self.active || !other.active {
            return None;
        }

        for element in &self.elements {
            for other_element in &other.elements {
                if self.offset + element.offset + offset == other.offset + other_element.offset {
                    return Some((element.offset, other_element.offset));
                }
            }
        }

        None
    }

    // Try to bond two molecules together
    // Offset is between the centers of the molecules
    // Updates and returns true if successful
    fn try_bond(&mut self, offset: Point, other: &Molecule) -> bool {
        if !self.active || !other.active {
            return false;
        }

        let mut bound = false;

        // Make local mutable copies
        let mut other = other.clone();

        // Go through each molecule pairwise
        for a in self.elements.iter_mut() {
            for b in other.elements.iter_mut() {
                let real_a = self.offset + offset + a.offset;
                let real_b = other.offset + b.offset;

                // Not adjacent
                if real_a.manhattan_distance(real_b) != 1 {
                    continue;
                }

                // Not enough free electrons
                if a.free_electrons == 0 || b.free_electrons == 0 {
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

    fn recenter(&mut self, new_zero: Point) {
        self.offset = self.offset + new_zero;

        for element in &mut self.elements {
            element.offset = element.offset - new_zero;
        }
        for dst_bond in &mut self.bonds {
            dst_bond.a = dst_bond.a - new_zero;
            dst_bond.b = dst_bond.b - new_zero;
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
            allow_multiple: false,
            solutions: vec![],
        };

        assert!(a.intersects_wall(Point { x: 1, y: 0 }, &map));
        assert!(!a.intersects_wall(Point { x: 1, y: 1 }, &map));
    }

    #[test]
    fn test_molecule_intersection() {
        let a = Molecule::new(Point::ZERO, ElementKind::Hydrogen);
        let b = Molecule::new(Point { x: 1, y: 0 }, ElementKind::Hydrogen);

        assert!(a.intersects(Point { x: 1, y: 0 }, &b).is_some());
        assert!(a.intersects(Point { x: 0, y: 0 }, &b).is_none());
        assert!(a.intersects(Point { x: 0, y: 1 }, &b).is_none());
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
        self.__try_move_recursive__(map, index, offset, true)
    }

    fn __try_move_recursive__(
        &mut self,
        map: &Map,
        index: usize,
        offset: Point,
        first: bool,
    ) -> bool {
        let original_molecules = self.molecules.clone();
        let mut moved_early = false;

        // Collect all map modifiers that we are trying to cross (this may take multiple passes)
        let mut modifiers_applied = Vec::new();
        loop {
            let mut bond_to_modify = None;

            // We want to handle bonds from the closest to origin first
            // This will help with cases with multiple rotators when we only want to hit the 'first' one
            let mut sorted_bonds = self.molecules[index]
                .bonds
                .iter()
                .enumerate()
                .collect::<Vec<_>>();

            sorted_bonds.sort_by_key(|(_, bond)| {
                bond.a.manhattan_distance(Point::ZERO) + bond.b.manhattan_distance(Point::ZERO)
            });

            'find_bond: for (bond_index, bond) in sorted_bonds {
                let mut sorted_modifiers = map.modifiers.clone();

                sorted_modifiers.sort_by_key(|m| {
                    let real_bond = *bond + self.molecules[index].offset;
                    real_bond.a.manhattan_distance(m.location)
                        + real_bond.b.manhattan_distance(m.location)

                    // TODO: Do we still need to sort by type?
                });

                for modifier in sorted_modifiers {
                    if modifiers_applied.contains(&modifier) {
                        continue;
                    }

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

                    // Try to performe a split at that bond
                    // If this returns None, it means we have a ring, so have nothing else to do
                    // Otherwise, update our molecule list
                    if let Some((part_a, part_b)) = self.split_at_bond(index, bond_index) {
                        self.molecules[index] = part_a;
                        self.molecules.push(part_b);
                    };
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
                ModifierKind::Rotate => {
                    // When rotating, the rest of the molecule will move as expected
                    // But the bond with the primary will 'wrap' around
                    // I expect I will get this wrong

                    // Split the molecule into two parts, (temporarily) removing the rotate bond
                    // The part with the old primary will move as expected
                    // The other part will be 'pulled' along the bond
                    // Both have to move without colliding
                    // And then we'll merge them and put the bond back

                    // Part b will move 'along' the path of the original bond
                    let parts = self.split_at_bond(index, bond_index);

                    // If it doesn't split, we have a ring; rotator won't work
                    if parts.is_none() {
                        self.molecules = original_molecules;
                        return false;
                    }

                    let (part_a, part_b) = parts.unwrap();

                    // Determine how which way part b will remove because we'll move a and b shortly

                    // Find the half of the bond that is still in a
                    let part_a_el_of_bond = part_a.offset
                        + part_a
                            .elements
                            .iter()
                            .find(|el| {
                                part_a.offset + el.offset == part_a.offset + bond.a
                                    || part_a.offset + el.offset == part_a.offset + bond.b
                            })
                            .expect("couldn't find bond in part a")
                            .offset;

                    let part_b_el_of_bond = part_b.offset
                        + part_b
                            .elements
                            .iter()
                            .find(|el| {
                                part_b.offset + el.offset == part_a.offset + bond.a
                                    || part_b.offset + el.offset == part_a.offset + bond.b
                            })
                            .expect("couldn't find bond in part b")
                            .offset;

                    // Determine which side of the rotate modifier that element is in
                    let left_side = part_a_el_of_bond.x == modifier.location.x;
                    let top_side = part_a_el_of_bond.y == modifier.location.y;
                    let moving_horizontal = offset.x != 0;

                    // And use that information to figure out how the b half of the bond will move
                    let new_b_offset = match (left_side, top_side, moving_horizontal) {
                        (true, true, true) => Point { x: 0, y: -1 },
                        (true, true, false) => Point { x: -1, y: 0 },
                        (true, false, true) => Point { x: 0, y: 1 },
                        (true, false, false) => Point { x: -1, y: 0 },
                        (false, true, true) => Point { x: 0, y: -1 },
                        (false, true, false) => Point { x: 1, y: 0 },
                        (false, false, true) => Point { x: 0, y: 1 },
                        (false, false, false) => Point { x: 1, y: 0 },
                    };

                    // Part a becomes the new primary, b is added
                    self.molecules[index] = part_a;
                    let part_b_index = self.molecules.len();
                    self.molecules.push(part_b);

                    // Part a contains the original primary, it moves in the original direction
                    // Disable collision checking with b to allow tail chasing
                    self.molecules[part_b_index].active = false;
                    let moved = self.__try_move_recursive__(map, index, offset, false);
                    if !moved {
                        self.molecules = original_molecules;
                        return false;
                    }
                    self.molecules[part_b_index].active = true;

                    // Part b contains the 'other' half which moves along the bond (as calculated earlier)
                    // Likewise, disable collision checking with a
                    self.molecules[index].active = false;
                    let moved = self.__try_move_recursive__(map, part_b_index, new_b_offset, false);
                    if !moved {
                        self.molecules = original_molecules;
                        return false;
                    }
                    self.molecules[index].active = true;

                    // Once they've both moved, make sure they're non intersecting
                    if self.molecules[index]
                        .intersects(Point::ZERO, &self.molecules[part_b_index])
                        .is_some()
                    {
                        self.molecules = original_molecules;
                        return false;
                    }

                    // Combine b into a
                    let mut part_a = self.molecules[index].clone();
                    let part_b = self.molecules[part_b_index].clone();

                    for element in part_b.elements {
                        part_a.elements.push(Element {
                            kind: element.kind,
                            offset: part_b.offset - part_a.offset + element.offset,
                            free_electrons: element.free_electrons,
                        });
                    }

                    for bond in part_b.bonds {
                        part_a.bonds.push(Bond {
                            a: part_b.offset - part_a.offset + bond.a,
                            b: part_b.offset - part_a.offset + bond.b,
                            count: bond.count,
                        });
                    }

                    // Create the new bond, this will be based on the part_a_el_of_bond and the new_b_offset
                    let new_bond = Bond {
                        a: part_a_el_of_bond + offset - part_a.offset,
                        b: part_b_el_of_bond + new_b_offset - part_a.offset,
                        count: bond.count,
                    };

                    // Validate that we didn't create a screwy bond
                    assert!(part_a.elements.iter().any(|el| el.offset == new_bond.a));
                    assert!(part_a.elements.iter().any(|el| el.offset == new_bond.b));

                    // Validate we don't have any overlapping bonds
                    assert!(!part_a.bonds.iter().enumerate().any(|(i, b1)| part_a
                        .bonds
                        .iter()
                        .enumerate()
                        .any(|(j, b2)| i != j
                            && ((b1.a == b2.a && b1.b == b2.b)
                                || (b1.a == b2.b && b1.b == b2.a)))));

                    part_a.bonds.push(new_bond);

                    self.molecules[index] = part_a;
                    self.molecules[part_b_index].active = false;

                    // We've already moved both parts, so don't move again
                    moved_early = true;

                    // HACK: Only apply one rotate per move
                    break;
                }
            }
        }

        if !moved_early {
            // Make sure we won't hit a wall
            // Check each moving molecule to see if it would hit a wall
            if self.molecules[index].intersects_wall(offset, map) {
                self.molecules = original_molecules;
                return false;
            }

            // Try to update each molecule we're pushing on
            'moving: loop {
                for other_index in 0..self.molecules.len() {
                    if other_index == index {
                        continue;
                    }

                    if let Some((_, other_intersection_offset)) =
                        self.molecules[index].intersects(offset, &self.molecules[other_index])
                    {
                        // HACK: Don't try to move the primary molecule more than once
                        // This can happen if the primary wraps around another, like:
                        /*
                        O-O-O
                        |   |
                        H h H
                        */
                        // This can happen for non-primaries, but I don't have a good way to fix that
                        // Assume if we made it this far, the primary *can* move, but don't actually do it here
                        if other_index == 0 {
                            continue;
                        }

                        // Recenter the other point so that the intersection is at 0,0
                        // This will properly handle a molecule being pushed across a split
                        // TODO: This won't properly work if multiple points are being pushed since only one is returned
                        self.molecules[other_index].recenter(other_intersection_offset);

                        let moved = self.__try_move_recursive__(map, other_index, offset, false);
                        if !moved {
                            self.molecules = original_molecules;
                            return false;
                        }

                        continue 'moving;
                    }
                }

                break 'moving;
            }

            // Verify that after pushing, we're no longer intersecting anything
            for i in 0..self.molecules.len() {
                if i == index {
                    continue;
                }

                if self.molecules[index]
                    .intersects(offset, &self.molecules[i])
                    .is_some()
                {
                    self.molecules = original_molecules;
                    return false;
                }
            }

            // Finally, update our own location
            self.molecules[index].offset = self.molecules[index].offset + offset;
        }

        // If we're the first call, re-apply bonds and remove inactive molecules
        if first {
            'bonding: loop {
                for i in 0..self.molecules.len() {
                    for j in 0..self.molecules.len() {
                        if i == j {
                            continue;
                        }

                        let mut primary = self.molecules[i].clone();
                        if primary.try_bond(Point::ZERO, &self.molecules[j]) {
                            self.molecules[i] = primary;
                            self.molecules[j].active = false;
                            continue 'bonding;
                        }
                    }
                }

                break 'bonding;
            }

            self.molecules.retain(|m| m.active);
        }

        true
    }

    fn split_at_bond(&mut self, index: usize, bond_index: usize) -> Option<(Molecule, Molecule)> {
        // Create a new primary half of the molecule
        let mut part_a = self.molecules[index].clone();
        part_a.bonds.remove(bond_index);

        // Flood fill to determine what elements and bonds are part of part A
        let mut connected_elements = Vec::new();
        let mut connected_bonds = Vec::new();
        let mut todo = vec![Point::ZERO];
        let mut done = vec![];

        // Keep going until we find all connected elements
        while let Some(pt) = todo.pop() {
            done.push(pt);

            // If any element matches the point, add to connected
            for (i, element) in part_a.elements.iter().enumerate() {
                if element.offset == pt && !connected_elements.contains(&i) {
                    connected_elements.push(i);
                }
            }

            // If any bond matches the point, add to connected
            // Then, if the other end isn't already accounted for, add it to the list to do
            for (i, src_bond) in part_a.bonds.iter().enumerate() {
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

        // If we are still connected to everything, this split will not create a new molecule
        // This can happen if you have a ring
        if connected_elements.len() == part_a.elements.len() {
            self.molecules[index] = part_a;
            return None;
        }

        // Create the second element
        // The first will remove anything that is *not* connected
        // The second (this one) will remove anything that *is* connected
        let mut part_b = part_a.clone();

        // Sort the elements and bonds and remove from the end so the indexes are correct
        connected_elements.sort();
        for i in connected_elements.iter().rev() {
            part_b.elements.remove(*i);
        }
        connected_bonds.sort();
        for i in connected_bonds.iter().rev() {
            part_b.bonds.remove(*i);
        }

        // Likewise, remove anything not in the list from the end to preserve indexes
        for i in (0..part_a.elements.len()).rev() {
            if !connected_elements.contains(&i) {
                part_a.elements.remove(i);
            }
        }
        for i in (0..part_a.bonds.len()).rev() {
            if !connected_bonds.contains(&i) {
                part_a.bonds.remove(i);
            }
        }

        // Part b no longer has an element at (0,0), choose one arbitrarily and recenter
        part_b.recenter(part_b.elements[0].offset);

        // Return the new parts
        Some((part_a, part_b))
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
- - -",
        );

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
- - -",
        );

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

    #[test]
    fn test_push_across_splitter() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
- - h - C
   / 
- - h - -",
        );

        assert!(state.try_move(&map, 0, Point { x: -1, y: 0 }));
        assert!(state.try_move(&map, 0, Point { x: -1, y: 0 }));

        // Should end up with a single molecule with two C/H bonds (one left, one down)
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 3);
    }

    #[test]
    fn test_push_across_splitter_inverse() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
- - h - -
   / 
- - h - C",
        );

        assert!(state.try_move(&map, 0, Point { x: -1, y: 0 }));
        assert!(state.try_move(&map, 0, Point { x: -1, y: 0 }));

        // Should end up with a single molecule with two C/H bonds (one left, one down)
        assert_eq!(state.molecules.len(), 1);
        assert_eq!(state.molecules[0].elements.len(), 3);
    }

    #[test]
    fn test_split_and_push() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
- - - h
   /
- c - O
   /
- - - h",
        );

        assert!(state.try_move(&map, 0, Point { x: -1, y: 0 }));
        assert!(state.try_move(&map, 0, Point { x: -1, y: 0 }));
    }

    #[test]
    fn test_rotate_simple() {
        use super::*;

        let (map, state) = Map::load(
            "\
v2
H - 
 @
h -",
        );

        // Order is:
        // - Primary hydrogen (absolute position)
        // - Other hydrogen (absolute position)
        // - The way the primary tries to move
        // - The ending relative offset of the secondary hydrogen
        let tests = [
            // Vertical, moving right, primary is on the top
            ((0, 0), (0, 1), (1, 0), (-1, 0)),
            // Vertical, moving right, primary is on the bottom
            ((0, 1), (0, 0), (1, 0), (-1, 0)),
            // Vertical, moving left, primary is on the top
            ((1, 0), (1, 1), (-1, 0), (1, 0)),
            // Vertical, moving left, primary is on the bottom
            ((1, 1), (1, 0), (-1, 0), (1, 0)),
            // Horizontal, moving down, primary is on the left
            ((0, 0), (1, 0), (0, 1), (0, -1)),
            // Horizontal, moving down, primary is on the right
            ((1, 0), (0, 0), (0, 1), (0, -1)),
            // Horizontal, moving up, primary is on the left
            ((0, 1), (1, 1), (0, -1), (0, 1)),
            // Horizontal, moving up, primary is on the right
            ((1, 1), (0, 1), (0, -1), (0, 1)),
        ];

        for (_i, (h0, h1, offset, end_offset)) in tests.into_iter().enumerate() {
            let h0: Point = h0.into();
            let h1: Point = h1.into();
            let offset: Point = offset.into();
            let end_offset: Point = end_offset.into();

            // Create the new test case from the original state
            let mut state = state.clone();
            state.molecules[0].offset = h0;
            state.molecules[0].elements[0].offset = Point::ZERO;
            state.molecules[0].elements[1].offset = h1 - h0;
            state.molecules[0].bonds[0].a = Point::ZERO;
            state.molecules[0].bonds[0].b = h1 - h0;

            // Try to move the requested direction
            assert!(state.try_move(&map, 0, offset));

            // Verify that the molecule itself and the secondary hydrogen moved as expected
            assert_eq!(state.molecules[0].offset, h0 + offset);
            assert_eq!(state.molecules[0].elements[0].offset, Point::ZERO);
            assert_eq!(state.molecules[0].elements[1].offset, end_offset);

            // Verify we didn't keep any extra molecules around or lose any bonds
            assert_eq!(state.molecules.len(), 1);
            assert_eq!(state.molecules[0].elements.len(), 2);
            assert_eq!(state.molecules[0].bonds.len(), 1);
        }
    }

    #[test]
    fn test_rotate_chain_succeed() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
- - h

- O o
   @
- - -",
        );

        assert!(state.try_move(&map, 0, Point { x: 0, y: 1 }));

        assert_eq!(state.molecules[0].offset, Point { x: 1, y: 2 });
        assert_eq!(state.molecules[0].elements[0].offset, Point { x: 0, y: 0 });
        assert_eq!(state.molecules[0].elements[1].offset, Point { x: 0, y: -1 });
        assert_eq!(state.molecules[0].elements[2].offset, Point { x: 0, y: -2 });
    }

    #[test]
    fn test_rotate_chain_fail() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
- - h

- o O
   @
- - -",
        );

        assert!(!state.try_move(&map, 0, Point { x: 0, y: 1 }));
    }

    #[test]
    fn test_rotate_chase_tail() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
- - - h 
 @
- H n n",
        );

        assert!(state.try_move(&map, 0, Point { x: -1, y: 0 }));
        assert!(state.try_move(&map, 0, Point { x: 0, y: -1 }));
        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));
    }

    #[test]
    fn test_rotate_tail_part() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
o H - -
 @ 
h - - -",
        );

        assert!(state.try_move(&map, 0, Point { x: 1, y: 0 }));

        assert_eq!(state.molecules[0].offset, Point { x: 2, y: 0 });
        assert_eq!(state.molecules[0].elements[0].offset, Point { x: 0, y: 0 });
        assert_eq!(state.molecules[0].elements[1].offset, Point { x: -1, y: 0 });
        assert_eq!(state.molecules[0].elements[2].offset, Point { x: -2, y: 0 });
    }

    #[test]
    fn test_rotate_in_a_corner() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
- c
   
- c
 @ 
- N",
        );
        assert!(state.try_move(&map, 0, Point { x: -1, y: 0 }));
        let actual = state.stringify(&map);

        let (map, state) = Map::load(
            "\
v2
- -
  
- c
 @
N c",
        );
        let expected = state.stringify(&map);

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_not_double_rotate() {
        use super::*;

        let (map, mut state) = Map::load(
            "\
v2
- n
 @
- c
  
- c
 @
- N",
        );
        assert!(state.try_move(&map, 0, Point { x: -1, y: 0 }));
        let actual = state.stringify(&map);

        let (map, state) = Map::load(
            "\
v2
- -
 @
- n
  
- c
 @
N c",
        );
        let expected = state.stringify(&map);

        assert_eq!(actual, expected);
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

        // If it's a board that allows multiples, anything goes
        if map.allow_multiple {
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
        self.molecules.iter().all(|m| m.free_electrons() == 0)
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
            if !molecule.active {
                continue;
            }

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

    log::info!(
        "Initial state (multiple={}):\n{}",
        map.allow_multiple,
        molecules.stringify(&map)
    );

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

    // If there is an arg, assume it's a test case and run it
    if std::env::args().len() > 1 {
        std::env::args().skip(1).for_each(|instructions| {
            let (map, mut state) = Map::load(&input);
            println!("Initial state:\n{}", state.stringify(&map));

            for c in instructions.chars() {
                let offset = match c {
                    'W' => Step::North,
                    'S' => Step::South,
                    'D' => Step::East,
                    'A' => Step::West,
                    _ => panic!("invalid instruction: {}", c),
                };

                assert!(state.try_move(&map, 0, offset.into()));
                println!("Move: {}, state:\n{}", c, state.stringify(&map));
                // for (i, molecule) in state.molecules.iter().enumerate() {
                //     println!("Molecule {}: {:?}", i, molecule);
                // }
                println!("Is solved? {}", state.is_solved(&map));
            }
        });
        return;
    }

    // Otherwise, try to solve the input
    let solution = solve(&input);

    if solution.is_none() {
        panic!("no solution found");
    }

    let solution = solution.unwrap();
    println!("{}", solution);
}

#[cfg(test)]
mod test_solutions {
    use std::{fs::File, io::Read, sync::mpsc, thread};

    use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

    use super::*;

    #[test]
    fn test_all_solutions() {
        // Timeout after 1 second or SOKOBOND_TEST_TIMEOUT if set
        let timeout = std::time::Duration::from_secs(
            std::env::var("SOKOBOND_TEST_TIMEOUT")
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(1),
        );

        // Collect all tests to run in order
        let mut test_files = Vec::new();
        for entry in std::fs::read_dir("data/sokobond").unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if !path.is_dir() {
                continue;
            }

            for entry in std::fs::read_dir(&path).unwrap() {
                let entry = entry.unwrap();
                let path = entry.path();
                if path.extension().unwrap() != "txt" {
                    continue;
                }

                test_files.push(path);
            }
        }
        test_files.sort();

        // Run each test with timeout
        #[derive(Debug, Clone, Eq, PartialEq)]
        enum TestResult {
            Success,
            NoSolution,
            InvalidSolution(String),
            TimedOut,
        }

        let results = test_files
            .par_iter()
            .map(move |path| {
                let mut file = File::open(&path).unwrap();
                let mut input = String::new();
                file.read_to_string(&mut input).unwrap();

                let (map, _) = Map::load(&input);

                let (tx, rx) = mpsc::channel();
                thread::spawn(move || {
                    let solution = solve(&input);
                    match tx.send(solution) {
                        Ok(_) => {}
                        Err(_) => {}
                    } // I don't actually care if this succeeds, but need to consume it
                });

                match rx.recv_timeout(timeout) {
                    Ok(solution) => {
                        if solution.is_none() {
                            log::debug!("No solution: {:?}", path);
                            return TestResult::NoSolution;
                        }
                        let solution = solution.unwrap();

                        if !map.solutions.contains(&solution) {
                            log::debug!("Invalid solution ({}): {:?}", solution, path);
                            return TestResult::InvalidSolution(solution);
                        }

                        log::debug!("Solved: {:?}", path);
                        return TestResult::Success;
                    }
                    Err(_) => {
                        log::debug!("Timed out: {:?}", path);
                        return TestResult::TimedOut;
                    }
                }
            })
            .collect::<Vec<_>>();

        // Print out the results
        if results.iter().any(|r| *r == TestResult::TimedOut) {
            println!("\nTimed out tests:");
            for (path, result) in test_files.iter().zip(results.iter()) {
                if *result == TestResult::TimedOut {
                    println!("  {:?}", path);
                }
            }
        }

        if results.iter().any(|r| *r == TestResult::NoSolution) {
            println!("\nUnsolved tests:");
            for (path, result) in test_files.iter().zip(results.iter()) {
                if *result == TestResult::NoSolution {
                    println!("  {:?}", path);
                }
            }
        }

        if results.iter().any(|r| {
            if let TestResult::InvalidSolution(_) = r {
                true
            } else {
                false
            }
        }) {
            println!("\nFailed tests:");
            for (path, result) in test_files.iter().zip(results.iter()) {
                if let TestResult::InvalidSolution(solution) = result {
                    println!("  {:?} -> {:?}", path, solution);
                }
            }
        }

        let perfect = results
            .iter()
            .all(|r| *r == TestResult::Success || *r == TestResult::TimedOut);
        if !perfect {
            println!();
        }
        assert!(perfect);
    }
}
