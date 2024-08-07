use std::hash::Hash;
use std::io;

use fxhash::FxHashSet;
use lazy_static::lazy_static;
use serde::{Deserialize, Serialize};

use solver::{Point, Solver, State};

lazy_static! {
    static ref FLOODFILL_VALIDATOR: bool =
        std::env::var("COSMIC_EXPRESS_FLOODFILL_VALIDATOR").is_ok();

    static ref HEURISTIC_COUNT_ENTITIES: bool =
        std::env::var("COSMIC_EXPRESS_HEURISTIC_COUNT_ENTITIES").is_ok();
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize, Hash)]
enum Color {
    Purple,
    Orange,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct CosmicExpressGlobal {
    width: isize,
    height: isize,
    length: isize,

    entrances: Vec<Point>,
    exits: Vec<Point>,
    walls: Vec<Point>,
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct CosmicExpressLocal {
    path: Vec<Point>,
    seats: Vec<Option<Color>>,
    aliens: Vec<(Point, Color)>,
    houses: Vec<(Point, Color)>,
}

impl State<CosmicExpressGlobal, ()> for CosmicExpressLocal {
    fn is_valid(&self, g: &CosmicExpressGlobal) -> bool {
        if !*FLOODFILL_VALIDATOR {
            return true;
        }

        // Flood fill from the current head, stopping at all walls and current path
        // If we can't reach all aliens, houses, and at least one exit, it's invalid

        let mut reachable = FxHashSet::default();
        let mut to_check = vec![self.path.last().unwrap().clone()];

        while let Some(p) = to_check.pop() {
            if reachable.contains(&p) {
                continue;
            }
            reachable.insert(p);

            // Flood fill all empty points
            for neighbor in p.neighbors().into_iter() {
                if reachable.contains(&neighbor) || to_check.contains(&neighbor) {
                    continue;
                }

                if neighbor.x < 1 || neighbor.y < 1 || neighbor.x > g.width || neighbor.y > g.height
                {
                    continue;
                }

                if g.walls.contains(&neighbor) || self.path.contains(&neighbor) {
                    continue;
                }

                to_check.push(neighbor);
            }
        }

        // Expand all points by one
        let mut expanded = FxHashSet::default();
        for p in reachable.iter() {
            expanded.insert(*p);
            for neighbor in p.neighbors().into_iter() {
                expanded.insert(neighbor);
            }
        }
        let reachable = expanded;

        // All aliens and houses must be reachable
        if self.aliens.iter().any(|(p, _)| !reachable.contains(p)) {
            return false;
        }

        if self.houses.iter().any(|(p, _)| !reachable.contains(p)) {
            return false;
        }

        // At least one exit must be reachable
        if !g.exits.iter().any(|p| reachable.contains(p)) {
            return false;
        }

        true
    }

    fn is_solved(&self, g: &CosmicExpressGlobal) -> bool {
        self.aliens.len() == 0
            && self.houses.len() == 0
            && g.exits.contains(&self.path.last().unwrap())
    }

    fn next_states(&self, g: &CosmicExpressGlobal) -> Option<Vec<(i64, (), Self)>> {
        let mut result = Vec::new();

        // If we don't have a path yet, start with each entrance
        if self.path.len() == 0 {
            for entrance in g.entrances.iter() {
                let mut new_local = self.clone();
                new_local.path.push(*entrance);

                result.push((1, (), new_local));
            }
        }
        // If we do have a path, expand the last node
        else {
            'neighbor: for p in self.path.last().unwrap().neighbors() {
                // Validate that the next point is valid
                // Cannot leave the bounds unless on entrance or exit
                if !(g.entrances.contains(&p) || g.exits.contains(&p))
                    && (p.x < 1 || p.y < 1 || p.x > g.width || p.y > g.height)
                {
                    continue 'neighbor;
                }

                // Cannot move onto walls
                if g.walls.contains(&p) {
                    continue 'neighbor;
                }

                // Cannot visit the same tile more than once
                for p2 in self.path.iter() {
                    if &p == p2 {
                        continue 'neighbor;
                    }
                }

                // Assume we can move, create the new state
                let mut new_local = self.clone();
                new_local.path.push(p);

                // Update each seat
                for (seat_index, seat_point) in self.path.iter().rev().skip(1).enumerate() {
                    // If we're over the end of the seats, we're done
                    if seat_index >= new_local.seats.len() {
                        break;
                    }

                    let mut seat_contents = new_local.seats[seat_index];

                    // Full seats next to the correct house; drop it off
                    if let Some(seat_color) = seat_contents {
                        for (house_index, (house_point, house_color)) in
                            self.houses.iter().enumerate()
                        {
                            if seat_point.manhattan_distance(*house_point) == 1
                                && seat_color == *house_color
                            {
                                new_local.houses.remove(house_index);
                                new_local.seats[seat_index] = None;
                                seat_contents = None;
                                break;
                            }
                        }
                    }

                    // Empty seats next to an alien; pick it up
                    // TODO: Multiple loads are not possible
                    if seat_contents.is_none() {
                        for (alien_index, (alien_point, alien_color)) in
                            self.aliens.iter().enumerate()
                        {
                            if seat_point.manhattan_distance(*alien_point) == 1 {
                                new_local.aliens.remove(alien_index);
                                new_local.seats[seat_index] = Some(*alien_color);
                                break;
                            }
                        }
                    }
                }

                result.push((1, (), new_local));
            }
        }

        // We'll always have nodes, so always return
        // We're relying on is_valid to filter impossible states this time
        Some(result)
    }

    fn stringify(&self, g: &CosmicExpressGlobal) -> String {
        let output_width = g.width + 2;
        let output_height = g.height + 2;
        let mut output_chars = vec![' '; (output_width * output_height) as usize];

        // Fill non-borders (1 offset)
        for y in 0..g.height {
            for x in 0..g.width {
                output_chars[((y + 1) * output_width + x + 1) as usize] = '.';
            }
        }

        // Add entrance and exit(s)
        for p in g.entrances.iter() {
            output_chars[(p.y * output_width + p.x) as usize] = '[';
        }
        for p in g.exits.iter() {
            output_chars[(p.y * output_width + p.x) as usize] = ']';
        }

        // Add path characters
        for (i, p) in self.path.iter().enumerate() {
            let p_before = self.path.get(i - 1);
            let p_after = self.path.get(i + 1);

            if p_before.is_some() && p_after.is_some() {
                let deltas = (
                    (p.x - p_before.unwrap().x, p.y - p_before.unwrap().y),
                    (p_after.unwrap().x - p.x, p_after.unwrap().y - p.y),
                );

                let c = match deltas {
                    ((1, 0), (1, 0)) => '─',
                    ((-1, 0), (-1, 0)) => '─',

                    ((0, 1), (0, 1)) => '│',
                    ((0, -1), (0, -1)) => '│',

                    ((1, 0), (0, -1)) => '┘',
                    ((0, 1), (-1, 0)) => '┘',

                    ((0, 1), (1, 0)) => '└',
                    ((-1, 0), (0, -1)) => '└',

                    ((1, 0), (0, 1)) => '┐',
                    ((0, -1), (-1, 0)) => '┐',

                    ((0, -1), (1, 0)) => '┌',
                    ((-1, 0), (0, 1)) => '┌',

                    _ => 'X',
                };
                output_chars[(p.y * output_width + p.x) as usize] = c;
            } else {
                output_chars[(p.y * output_width + p.x) as usize] = 'O';
            }
        }

        // Add walls
        for p in g.walls.iter() {
            output_chars[(p.y * output_width + p.x) as usize] = '#';
        }

        // Add remaining aliens
        for (p, c) in self.aliens.iter() {
            output_chars[(p.y * output_width + p.x) as usize] = ('a' as u8 + *c as u8) as char;
        }

        // Add remaining houses
        for (p, c) in self.houses.iter() {
            output_chars[(p.y * output_width + p.x) as usize] = ('A' as u8 + *c as u8) as char;
        }

        // Convert to a string
        let mut output = String::new();
        for y in 0..output_height {
            for x in 0..output_width {
                output.push(output_chars[(y * output_width + x) as usize]);
            }
            output.push('\n');
        }
        output
    }

    fn heuristic(&self, global: &CosmicExpressGlobal) -> i64 {
        let mut heuristic = 0;

        if *HEURISTIC_COUNT_ENTITIES {
            heuristic += ((self.aliens.len() + self.houses.len()) as isize
                * global.width.max(global.height)) as i64;
        }

        heuristic
    }
}

// Used only for loading
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
enum Entity {
    Wall,
    Alien { color: Color },
    House { color: Color },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct CosmicExpressDefinition {
    width: isize,
    height: isize,
    length: isize,

    entrances: Vec<Point>,
    exits: Vec<Point>,

    entities: Vec<(Point, Entity)>,
}

fn main() {
    env_logger::init();

    let stdin = io::stdin().lock();
    let definition: CosmicExpressDefinition = serde_json::from_reader(stdin).unwrap();

    // Convert to a global and local

    // There is a hidden wall under all aliens and houses
    // TODO: Entities should have an is_wall property or something
    let walls = definition.entities.iter().map(|(p, _)| *p).collect();

    let global = CosmicExpressGlobal {
        width: definition.width,
        height: definition.height,
        length: definition.length,
        entrances: definition.entrances,
        exits: definition.exits,
        walls,
    };

    let seats = (0..definition.length).map(|_| None).collect();

    let aliens = definition
        .entities
        .iter()
        .filter_map(|(p, e)| {
            if let Entity::Alien { color } = e {
                Some((*p, *color))
            } else {
                None
            }
        })
        .collect();

    let houses = definition
        .entities
        .iter()
        .filter_map(|(p, e)| {
            if let Entity::House { color } = e {
                Some((*p, *color))
            } else {
                None
            }
        })
        .collect();

    let local = CosmicExpressLocal {
        path: Vec::new(),
        seats,
        aliens,
        houses,
    };

    let mut solver = Solver::new(global.clone(), local.clone());
    println!("{}", solver.stringify(&local));

    while let Some(state) = solver.next() {
        if solver.states_checked() % 100000 == 0 {
            println!("{}", state.stringify(&global));
            println!(
                "{} states, {} invalidated, {} seconds",
                solver.states_checked(),
                solver.states_invalidated(),
                solver.time_spent()
            );
        }
    }
    let solution = solver.get_solution();

    println!("{:?}", solution);
    if let Some(path) = solution {
        println!("{}", solver.stringify(&path));
    }

    println!(
        "{} states, {} seconds",
        solver.states_checked(),
        solver.time_spent()
    );
}
