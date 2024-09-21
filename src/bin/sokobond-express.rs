use std::{collections::HashMap, io::Read, process::exit};

use fxhash::FxHashSet;
use solver::{Direction, Point, Solver, State};

use anyhow::{anyhow, Result};

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
enum Element {
    Hydrogen,
    Oxygen,
    Nitrogen,
    Carbon,
}

impl TryFrom<&str> for Element {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> std::result::Result<Self, Self::Error> {
        match value {
            "hydrogen" => Ok(Element::Hydrogen),
            "oxygen" => Ok(Element::Oxygen),
            "nitrogen" => Ok(Element::Nitrogen),
            "carbon" => Ok(Element::Carbon),
            _ => Err(anyhow!("unknown element: {}", value)),
        }
    }
}

impl Element {
    fn valence(&self) -> usize {
        match self {
            Element::Hydrogen => 1,
            Element::Oxygen => 2,
            Element::Nitrogen => 3,
            Element::Carbon => 4,
        }
    }

    fn as_char(&self) -> char {
        match self {
            Element::Hydrogen => 'H',
            Element::Oxygen => 'O',
            Element::Nitrogen => 'N',
            Element::Carbon => 'C',
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
struct ElementData {
    element: Element,
    offset: Point,
    electrons: usize,
    holes: usize,
}
impl ElementData {
    fn atom(el: Element) -> ElementData {
        ElementData {
            element: el,
            offset: Point { x: 0, y: 0 },
            electrons: el.valence(),
            holes: 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct BondData {
    i: usize,
    j: usize,
    strength: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Molecule {
    pt: Point,
    active: bool,

    elements: Vec<ElementData>,
    bonds: Vec<BondData>,
}

impl Molecule {
    fn available_bonds(&self) -> usize {
        self.elements.iter().map(|e| e.electrons + e.holes).sum()
    }

    fn contains(&self, pt: &Point) -> bool {
        self.elements
            .iter()
            .find(|e| self.pt + e.offset == *pt)
            .is_some()
    }

    fn maybe_bond(&self, other: &Molecule) -> Option<Molecule> {
        // At least one molecule must be active to bond
        if !(self.active || other.active) {
            return None;
        }

        // Look for the first pair of points that are 1 apart and both have at least one free electron
        for (i, &self_el) in self.elements.iter().enumerate() {
            for (j, &other_el) in other.elements.iter().enumerate() {
                let map_self_pt = self.pt + self_el.offset;
                let map_other_pt = other.pt + other_el.offset;

                if map_self_pt.manhattan_distance(map_other_pt) != 1 {
                    continue;
                }

                let bond_strength = self_el.electrons.min(other_el.electrons);
                if bond_strength == 0 {
                    continue;
                }

                // If we made it here, create a new molecule that merges the two

                // Combine the new elements by:
                // - Removing electrons used in the bond
                // - Updating the offset from other_el to center on self_el
                let mut elements = Vec::new();

                for (k, el) in self.elements.iter().enumerate() {
                    let mut new_el = el.clone();

                    if k == i {
                        new_el.electrons -= bond_strength;
                    }

                    elements.push(new_el);
                }

                for (k, el) in other.elements.iter().enumerate() {
                    let mut new_el = el.clone();

                    new_el.offset = new_el.offset - self.pt + other.pt;

                    if k == j {
                        new_el.electrons -= bond_strength;
                    }

                    elements.push(new_el);
                }

                // Bonds are copied + the new one
                let mut bonds = self.bonds.clone();
                bonds.extend(other.bonds.iter().cloned());

                bonds.push(BondData {
                    i,
                    j: j + self.elements.len(),
                    strength: bond_strength,
                });

                return Some(Molecule {
                    pt: self.pt,
                    elements,
                    bonds,
                    active: true,
                });
            }
        }

        // If we made it this far, no matches
        None
    }
}

#[derive(Debug, Clone)]
struct Global {
    width: usize,
    height: usize,

    is_loop: bool,

    spawns: Vec<Point>,  // Where atoms/electrons start
    cuts: Vec<Point>,    // Where there are holes in the map
    walls: Vec<Point>,   // Solid parts of the map things can't go through
    crosses: Vec<Point>, // Points where the track can cross

    exit: (Point, Direction),
}

impl Global {
    // Nothing is at the point, so things can be placed there
    // This includes potentially hanging out over empty space
    fn is_empty(&self, pt: &Point) -> bool {
        !self.walls.contains(pt)
    }

    // Points that the track can be extended onto
    // Iota special case: the exit is always walkable
    fn is_walkable(&self, pt: &Point) -> bool {
        self.exit.0 == *pt
            || pt.x >= 0
                && pt.x < self.width as isize
                && pt.y >= 0
                && pt.y < self.height as isize
                && !self.walls.contains(&pt)
                && !self.cuts.contains(&pt)
                && !self.spawns.contains(pt)
    }

    // Points that can contains a final molecule
    fn is_exitable(&self, pt: &Point) -> bool {
        pt.x >= 0
            && pt.x < self.width as isize
            && pt.y >= 0
            && pt.y < self.height as isize
            && !self.walls.contains(&pt)
            && !self.cuts.contains(&pt)
            && !self.crosses.contains(&pt)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Local {
    head: Point,
    track: Vec<Point>,
    molecules: Vec<Molecule>,
    electrons: Vec<(Point, usize)>,
}

impl Local {
    // Try to extend the track in the given direction
    fn maybe_move(&mut self, d: Direction, global: &Global) -> bool {
        let new_head = self.head + d.into();

        // If we're moving back onto the track... don't do that
        if self.track.contains(&new_head) {
            if global.crosses.contains(&new_head) {
                // The exception is a cross allows exactly two crossings
                // There's no way to get into a cross more than twice, so just ignore them
            } else if global.is_loop && new_head == global.exit.0 {
                // Another special case is Iota levels where the start and exit are the same
            } else {
                return false;
            }
        }

        // The global state has to allow moving to that point
        if !global.is_walkable(&new_head) {
            return false;
        }

        // If we're moving onto the exit, we have to be going in the correct direction
        if new_head == global.exit.0 {
            if d != global.exit.1 {
                return false;
            }
        }

        // If we're moving off a cross, we have to continue straight
        if global.crosses.contains(&self.head) {
            let delta_1 = new_head - self.head;
            let delta_2 = self.head - self.track[self.track.len() - 1];

            let d1: Direction = delta_1.try_into().unwrap();
            let d2: Direction = delta_2.try_into().unwrap();

            if d1 != d2 {
                return false;
            }
        }

        // We can't move off of the exit
        // Unless it's also the starting point (Iota)
        if self.head == global.exit.0 && !self.track.is_empty() {
            return false;
        }

        // Find the molecule we're moving
        // TODO: Consider extract this into a method
        let moving_index = self.molecules.iter().position(|m| m.contains(&self.head));
        if moving_index.is_none() {
            return false;
        }
        let moving_index = moving_index.unwrap();

        // Move it, verify nothing collided
        self.molecules[moving_index].pt = self.molecules[moving_index].pt + d.into();

        // Once a molecule moves, it's active
        self.molecules[moving_index].active = true;

        if self.has_collision(global) {
            return false;
        }

        // Update the track
        self.track.push(self.head);
        self.head = new_head;

        // Check all molecules/electrons for bonding
        self.apply_electrons();
        self.apply_bonds();
        self.apply_self_bonds();

        // Success!
        true
    }

    // Check collisions between the current state of molecules and electrons
    fn has_collision(&self, global: &Global) -> bool {
        for (i, m) in self.molecules.iter().enumerate() {
            for el in m.elements.iter() {
                let map_pt = m.pt + el.offset;

                if !global.is_empty(&map_pt) {
                    return true;
                }

                // Hit another molecule
                for (j, m2) in self.molecules.iter().enumerate() {
                    if i == j {
                        continue;
                    }

                    for el2 in m2.elements.iter() {
                        let map_pt2 = m2.pt + el2.offset;

                        if map_pt == map_pt2 {
                            return true;
                        }
                    }
                }

                // Hit an electron
                if self.electrons.iter().any(|(pt, _)| *pt == map_pt) {
                    return true;
                }
            }
        }

        false
    }

    // Try to bond any molecules that are adjacent
    fn apply_bonds(&mut self) {
        // Keep looping until the simulation 'settles' (finds no new bonds)
        'bonding: loop {
            // Find the first indexes i,j that bond and the molecule they'd create
            let mut to_bond = None;
            'find_bond: for (i, m1) in self.molecules.iter().enumerate() {
                for (j, m2) in self.molecules.iter().enumerate() {
                    if i == j {
                        continue;
                    }

                    if let Some(new_m) = m1.maybe_bond(m2) {
                        to_bond = Some((i, j, new_m));
                        break 'find_bond;
                    }
                }
            }

            // If we have some, remove the old molecules, add the new one, and continue
            // The i,j removal is intentionally ordered to avoid invalidating the indexes
            // If to_remove is none, that means we've settled
            if let Some((i, j, new_m)) = to_bond {
                self.molecules.remove(i.max(j));
                self.molecules.remove(i.min(j));
                self.molecules.push(new_m);
            } else {
                break 'bonding;
            }
        }
    }

    // Apply additional bonds within active molecules
    fn apply_self_bonds(&mut self) {
        'self_bonding: loop {
            let mut to_bond = None;

            'find_bond: for (i, m) in self.molecules.iter().enumerate() {
                for (j, el1) in m.elements.iter().enumerate() {
                    for (k, el2) in m.elements.iter().enumerate() {
                        if j == k {
                            continue;
                        }

                        if el1.offset.manhattan_distance(el2.offset) != 1 {
                            continue;
                        }

                        let bond_strength = el1.electrons.min(el2.electrons);
                        if bond_strength == 0 {
                            continue;
                        }

                        to_bond = Some((i, j, k, bond_strength));
                        break 'find_bond;
                    }
                }
            }

            if let Some((i, j, k, bond_strength)) = to_bond {
                let m = &mut self.molecules[i];

                // We already had a bond
                if let Some(bond_index) = m.bonds.iter().position(|b| b.i == j && b.j == k) {
                    m.bonds[bond_index].strength += bond_strength;
                    m.elements[j].electrons -= bond_strength;
                    m.elements[k].electrons -= bond_strength;
                } else if let Some(bond_index) = m.bonds.iter().position(|b| b.i == k && b.j == j) {
                    m.bonds[bond_index].strength += bond_strength;
                    m.elements[j].electrons -= bond_strength;
                    m.elements[k].electrons -= bond_strength;
                } else {
                    // Didn't have one, Create a new bond!
                    m.bonds.push(BondData {
                        i: j,
                        j: k,
                        strength: bond_strength,
                    });

                    m.elements[j].electrons -= bond_strength;
                    m.elements[k].electrons -= bond_strength;
                }
            } else {
                break 'self_bonding;
            }
        }
    }

    // Any atoms that are missing electrons can grab them
    fn apply_electrons(&mut self) {
        // Each electron can only be grabbed once per tick
        let mut grabbed = FxHashSet::default();
        
        // Keep looping until we settle
        'electroning: loop {
            // Find the atom j of molecule i with a hole adjacent to electron k
            let mut to_electron = None;
            'find_electron: for (i, m) in self.molecules.iter().enumerate() {
                for (j, el) in m.elements.iter().enumerate() {
                    if el.holes == 0 {
                        continue;
                    }

                    for (k, (pt, _)) in self.electrons.iter().enumerate() {
                        if grabbed.contains(&((i, j, k))) {
                            continue;
                        }

                        if pt.manhattan_distance(m.pt + el.offset) == 1 {
                            to_electron = Some((i, j, k));
                            break 'find_electron;
                        }
                    }
                }
            }

            // Only move one electron at a time
            if let Some((i, j, k)) = to_electron {
                self.molecules[i].elements[j].electrons += 1;
                self.molecules[i].elements[j].holes -= 1;
                self.electrons[k].1 -= 1;

                if self.electrons[k].1 == 0 {
                    self.electrons.remove(k);
                }

                grabbed.insert((i, j, k));
            } else {
                break 'electroning;
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Step;

impl State<Global, Step> for Local {
    fn is_valid(&self, global: &Global) -> bool {
        if self.is_solved(global) {
            return true;
        }

        // If the molecule on the track has no more free electrons and we need to connect to more, fail
        // This improved beta 5 from 340s to 145s
        // TODO: This might not be valid for all levels
        let track_molecule = self
            .molecules
            .iter()
            .find(|m| m.contains(&self.head))
            .unwrap();

        if track_molecule.available_bonds() == 0
            && self.molecules.iter().any(|m| m.available_bonds() > 0)
        {
            return false;
        }

        // Floodfill validator
        // TODO: Add a flag for this
        if true {
            // Determine all points the track can still reach
            let mut reachable = vec![];
            let mut to_check = vec![self.head];

            while let Some(pt) = to_check.pop() {
                reachable.push(pt);

                for d in [
                    Direction::Up,
                    Direction::Right,
                    Direction::Down,
                    Direction::Left,
                ] {
                    let new_pt = pt + d.into();

                    if !global.is_walkable(&new_pt) {
                        continue;
                    }

                    if reachable.contains(&new_pt) {
                        continue;
                    }

                    if to_check.contains(&new_pt) {
                        continue;
                    }

                    // TODO: Is this just disabling this check for iota levels?
                    if self.track.contains(&new_pt)
                        && !global.crosses.contains(&new_pt)
                        && !global.is_loop
                    {
                        continue;
                    }

                    to_check.push(new_pt);
                }
            }

            // The radius of the tracked molecule is the largest manhattan distance from the track
            // This is an intentional simplification: assume the furthest 'out' atom can pick up more molecules because it's faster to calculate
            let tracked_molecule = self
                .molecules
                .iter()
                .find(|m| m.contains(&self.head))
                .unwrap();

            let radius = tracked_molecule
                .elements
                .iter()
                .map(|el| {
                    self.head
                        .manhattan_distance(tracked_molecule.pt + el.offset)
                })
                .max()
                .unwrap();

            // All points of all molecules must be within radius + 1 (for a new bond) of a reachable point
            for m in self.molecules.iter() {
                for el in m.elements.iter() {
                    let map_pt = m.pt + el.offset;

                    if !reachable
                        .iter()
                        .any(|pt| pt.manhattan_distance(map_pt) <= radius + 1)
                    {
                        return false;
                    }
                }
            }

            // All electrons must also be within radius + 1
            for (pt, _) in self.electrons.iter() {
                if !reachable
                    .iter()
                    .any(|pt2| pt.manhattan_distance(*pt2) <= radius + 1)
                {
                    return false;
                }
            }

            // The entrypoint to the exit must be reachable
            if !reachable.contains(&(global.exit.0 - global.exit.1.into())) {
                return false;
            }
        }

        true
    }

    // Solved if:
    // - there's one remaining molecule
    // - the molecule is at the exit
    // - the molecule has no free electrons
    fn is_solved(&self, global: &Global) -> bool {
        let m0 = &self.molecules[0];

        self.molecules.len() == 1
            && self.head == global.exit.0
            && m0.available_bonds() == 0
            && m0
                .elements
                .iter()
                .all(|e| global.is_exitable(&(m0.pt + e.offset)))
    }

    fn next_states(&self, global: &Global) -> Option<Vec<(i64, Step, Local)>> {
        let mut next_states = Vec::new();

        for d in &[
            Direction::Up,
            Direction::Right,
            Direction::Down,
            Direction::Left,
        ] {
            // We cannot double back
            let previous_delta = self.head - *self.track.last().unwrap();
            let previous_d: Direction = previous_delta.into();

            if previous_d == d.flip() {
                continue;
            }

            let mut new_state = self.clone();
            if new_state.maybe_move(*d, global) {
                next_states.push((1, Step, new_state));
            }
        }

        if next_states.is_empty() {
            None
        } else {
            Some(next_states)
        }
    }

    fn heuristic(&self, _g: &Global) -> i64 {
        let mut score = 0;

        // Calculate the distance between each molecule and the nearest other
        // Because molecules expand, this will over estimate
        // Filter i < j avoids double counting, it at least needs to be !-
        score += self
            .molecules
            .iter()
            .enumerate()
            .map(|(i, m)| {
                self.molecules
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| i < *j)
                    .map(|(_, m2)| m.pt.manhattan_distance(m2.pt))
                    .min()
                    .unwrap_or_default()
            })
            .sum::<isize>();

        // Add the distance from any molecule to the exit
        score += self
            .molecules
            .iter()
            .map(|m| m.pt.manhattan_distance(_g.exit.0))
            .min()
            .unwrap();

        score as i64
    }

    fn stringify(&self, g: &Global) -> String {
        let mut chars = HashMap::new();

        // Anything added here might overlap, so draw most important things last

        // Add floors (width and height but not cuts)
        for y in 0..g.height as isize {
            for x in 0..g.width as isize {
                if g.cuts.contains(&Point { x, y }) {
                    continue;
                }

                chars.insert(Point { x, y }, '.');
            }
        }

        // Add spawns
        for p in g.spawns.iter() {
            // chars.insert(*p, '⚬'); // Removed for backwards compatibility with testit
            chars.insert(*p, ' ');
        }

        // Add walls
        for p in g.walls.iter() {
            chars.insert(*p, '#');
        }

        // Add (unused) crossovers
        for p in g.crosses.iter() {
            chars.insert(*p, '╬');
        }

        // Add the first part of the track
        if self.track.len() > 0 {
            let p0 = self.track[0];
            let p1 = if self.track.len() == 1 {
                self.head
            } else {
                self.track[1]
            };

            let delta = p1 - p0;
            let d: Direction = delta.try_into().unwrap();

            let c = match d {
                Direction::Up => '↑',
                Direction::Down => '↓',
                Direction::Left => '←',
                Direction::Right => '→',
            };

            chars.insert(p0, c);
        }

        let mut track_and_head = self.track.clone();
        track_and_head.push(self.head);

        if track_and_head.len() > 2 {
            for i in 0..(track_and_head.len() - 2) {
                let delta_1 = track_and_head[i + 1] - track_and_head[i];
                let delta_2 = track_and_head[i + 2] - track_and_head[i + 1];

                let d1: Direction = delta_1.try_into().unwrap();
                let d2: Direction = delta_2.try_into().unwrap();

                let c = match (d1, d2) {
                    (Direction::Up, Direction::Up) => '│',
                    (Direction::Down, Direction::Down) => '│',
                    (Direction::Left, Direction::Left) => '─',
                    (Direction::Right, Direction::Right) => '─',

                    (Direction::Up, Direction::Left) => '┐',
                    (Direction::Up, Direction::Right) => '┌',

                    (Direction::Down, Direction::Left) => '┘',
                    (Direction::Down, Direction::Right) => '└',

                    (Direction::Left, Direction::Up) => '└',
                    (Direction::Left, Direction::Down) => '┌',

                    (Direction::Right, Direction::Up) => '┘',
                    (Direction::Right, Direction::Down) => '┐',

                    _ => panic!("invalid track: {d1:?}, {d2:?}"),
                };

                chars.insert(self.track[i + 1], c);
            }
        }

        // Draw the exit
        {
            let (exit, d) = g.exit;
            let c = match d {
                Direction::Up => '⇧',
                Direction::Down => '⇩',
                Direction::Left => '⇦',
                Direction::Right => '⇨',
            };

            chars.insert(exit, c);
        }

        // Add any electrons
        for (p, c) in self.electrons.iter() {
            chars.insert(*p, match c {
                1 => '+',
                2 => '⧺',
                3 => '⧻',
                _ => unimplemented!("can only print 1-3 electrons")
            });
        }

        // Add all molecules to the map
        for m in self.molecules.iter() {
            for el in m.elements.iter() {
                chars.insert(m.pt + el.offset, el.element.as_char());
            }
        }

        let min_x = chars.keys().map(|p| p.x).min().unwrap();
        let max_x = chars.keys().map(|p| p.x).max().unwrap();
        let min_y = chars.keys().map(|p| p.y).min().unwrap();
        let max_y = chars.keys().map(|p| p.y).max().unwrap();

        let mut result = String::new();

        for y in min_y..=max_y {
            for x in min_x..=max_x {
                let c = chars.get(&Point { x, y }).unwrap_or(&' ');
                result.push(*c);
            }
            result.push('\n');
        }

        result
    }
}

fn load(input: &str) -> Result<(Global, Local)> {
    let mut global = Global {
        width: 0,
        height: 0,

        is_loop: false,

        spawns: Vec::new(),
        cuts: Vec::new(),
        walls: Vec::new(),
        crosses: Vec::new(),

        exit: (Point { x: 0, y: 0 }, Direction::Up),
    };

    let mut local = Local {
        head: Point { x: 0, y: 0 },
        track: Vec::new(),
        molecules: Vec::new(),
        electrons: Vec::new(),
    };

    let mut first_move = Direction::Down;

    for line in input.lines() {
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts = line.split_whitespace().collect::<Vec<_>>();
        match parts[0] {
            "map" => {
                assert!(parts.len() == 3, "map <width> <height>");
                global.width = parts[1].parse()?;
                global.height = parts[2].parse()?;
            }
            "cut" => {
                assert!(parts.len() == 3, "cut <x> <y>");
                let x: isize = parts[1].parse()?;
                let y: isize = parts[2].parse()?;
                global.cuts.push(Point { x: x - 1, y: y - 1 });
            }
            "cross" => {
                assert!(parts.len() == 3, "cross <x> <y>");
                let x: isize = parts[1].parse()?;
                let y: isize = parts[2].parse()?;
                global.crosses.push(Point { x: x - 1, y: y - 1 });
            }
            "wall" => {
                assert!(parts.len() == 3, "wall <x> <y>");
                let x: isize = parts[1].parse()?;
                let y: isize = parts[2].parse()?;
                global.walls.push(Point { x: x - 1, y: y - 1 });
            }
            "track" => {
                assert!(parts.len() == 4, "track <x> <y> <direction>");
                let x: isize = parts[1].parse()?;
                let y: isize = parts[2].parse()?;

                local.head = Point { x: x - 1, y: y - 1 };
                first_move = parts[3].try_into()?;
            }
            "exit" => {
                assert!(parts.len() == 4, "exit <x> <y> <direction>");
                let x: isize = parts[1].parse()?;
                let y: isize = parts[2].parse()?;
                let d = parts[3].try_into()?;

                global.exit = (Point { x: x - 1, y: y - 1 }, d);
            }
            "atom" => {
                assert!(parts.len() == 4 || parts.len() == 5, "atom <x> <y> <element> <holes>?");
                let x: isize = parts[1].parse()?;
                let y: isize = parts[2].parse()?;
                let mut element = ElementData::atom(parts[3].try_into()?);

                if parts.len() == 5 {
                    let holes: usize = parts[4].parse()?;
                    assert!(
                        holes <= element.electrons,
                        "atom: cannot have holes > valence"
                    );

                    element.electrons -= holes;
                    element.holes = holes;
                }

                local.molecules.push(Molecule {
                    pt: Point { x: x - 1, y: y - 1 },
                    elements: vec![element],
                    bonds: Vec::new(),
                    active: false,
                });
            }
            "anion" => {
                panic!("anion has been removed; use atom");
            }
            "electron" | "electrons" => {
                assert!(parts.len() == 3 || parts.len() == 4, "electron <x> <y> <count>?");
                let x: isize = parts[1].parse()?;
                let y: isize = parts[2].parse()?;
                let count = if parts.len() == 4 {
                    parts[3].parse()?
                } else {
                    1
                };

                local.electrons.push((Point { x: x - 1, y: y - 1 }, count));
            }
            _ => return Err(anyhow!("unknown command: {}", parts[0])),
        }
    }

    // Add a spawn where each molecule/atom is since the track can't go there
    // This has to be different from cuts, because atoms *can* end there
    for m in local.molecules.iter() {
        for el in m.elements.iter() {
            global.spawns.push(m.pt + el.offset);
        }
    }

    for (pt, _) in local.electrons.iter() {
        global.spawns.push(*pt);
    }

    // If the track is a loop, flag it
    if local.head == global.exit.0 {
        global.is_loop = true;
    }

    // TODO: Validity checks
    // - All points should be within the bounds, not on walls or cuts
    // - No molecules should be colliding

    // We have to perform the first move here
    let success = local.maybe_move(first_move, &global);
    assert!(success, "failed to perform first move");

    Ok((global, local))
}

fn main() -> Result<()> {
    env_logger::init();

    let mut stdin = std::io::stdin();
    let mut input = String::new();
    stdin.read_to_string(&mut input)?;

    let (global, local) = load(&input)?;

    log::info!("Initial state:\n{}", local.stringify(&global));

    // If we have a sequence in args, run it
    // If there is an arg, assume it's a test case and run it
    if std::env::args().len() > 1 {
        let mut local = local;

        std::env::args().skip(1).for_each(|directions| {
            println!("Initial state:\n{}", local.stringify(&global));

            for c in directions.chars() {
                let d: Direction = c.try_into().expect("invalid direction");

                let success = local.maybe_move(d, &global);
                if !success {
                    println!("Failed to move in direction: {}", c);
                    exit(1);
                }

                println!("Move: {}, state:\n{}", c, local.stringify(&global));
                println!("Is solved? {}", local.is_solved(&global));
            }
        });
        return Ok(());
    }

    // Otherwise, solve it
    let mut solver = Solver::new(global.clone(), local);

    while let Some(state) = solver.next() {
        if solver.states_checked() % 100000 == 0 {
            log::debug!("=== Current state: ===\n{}", state.stringify(&global));
            log::debug!("{solver}");
        }
    }
    let solution = solver.get_solution();

    if let Some(solution) = solution {
        println!("{}", solver.stringify(&solution));

        let mut track = solution.track.clone();
        track.push(solution.head);

        let output = track
            .iter()
            .skip(1)
            .zip(track.iter().skip(2))
            .map(|(a, b)| {
                let delta = *b - *a;
                let d: Direction = delta.try_into().unwrap();

                match d {
                    Direction::Up => 'U',
                    Direction::Down => 'D',
                    Direction::Left => 'L',
                    Direction::Right => 'R',
                }
            })
            .collect::<String>();
        eprintln!("{}", output);

    } else {
        println!("{solver}\nNo solution found");
        exit(1);
    }

    Ok(())
}
