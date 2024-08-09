use std::hash::Hash;
use std::io;

use fxhash::FxHashSet;
use lazy_static::lazy_static;
use serde::{Deserialize, Serialize};

use solver::{Point, Solver, State};

lazy_static! {
    static ref DEBUG_PRINT: bool = std::env::var("COSMIC_EXPRESS_DEBUG_PRINT").is_ok();
    static ref FLOODFILL_VALIDATOR: bool =
        std::env::var("COSMIC_EXPRESS_FLOODFILL_VALIDATOR").is_ok();
    static ref HEURISTIC_COUNT_ENTITIES: bool =
        std::env::var("COSMIC_EXPRESS_HEURISTIC_COUNT_ENTITIES").is_ok();
    static ref HEURISTIC_NEAREST_HOUSE: bool =
        std::env::var("COSMIC_EXPRESS_HEURISTIC_NEAREST_HOUSE").is_ok();
    static ref USE_CUSTOM_HASH: bool =
        std::env::var("COSMIC_EXPRESS_CUSTOM_HASH").is_ok();
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize, Hash, Ord, PartialOrd)]
enum Color {
    Purple,
    Orange,
    Green,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct CosmicExpressGlobal {
    width: isize,
    height: isize,
    length: isize,

    entrance: Point,
    exit: Point,
    walls: Vec<Point>,
}

#[derive(Clone, PartialEq, Eq, Debug)]
struct CosmicExpressLocal {
    // The path the train has taken so far
    path: Vec<Point>,

    // The current seats of the train
    // Seats contains the color of the alien in the seat
    // Goop is a flag on if a seat has been 'gooped' by a green alien
    seats: Vec<Option<Color>>,
    seat_goop: Vec<bool>,

    // Remaining aliens that haven't been picked up / houses that haven't been delivered to
    aliens: Vec<(Point, Color)>,
    houses: Vec<(Point, Color)>,
}

impl Hash for CosmicExpressLocal {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        if *USE_CUSTOM_HASH {
            // Getting to exactly the same points a different way counts as equal
            // So long as the last step (where we'll expand from) is the same
            let mut path = self.path.clone();
            let last = path.pop().unwrap();

            path.sort();
            path.push(last);
            path.hash(state);

            // We don't care about which aliens since they don't move, just the list of remaining points
            // Sort so order is preserved
            let mut entities = Vec::with_capacity(self.aliens.len() + self.houses.len());
            for (p, _) in self.aliens.iter() {
                entities.push(p);
            }
            for (p, _) in self.houses.iter() {
                entities.push(p);
            }
            entities.sort();
            entities.hash(state);

            self.seats.hash(state);
            self.seat_goop.hash(state);
        } else {
            // Default hashing
            self.path.hash(state);
            self.seats.hash(state);
            self.seat_goop.hash(state);
            self.aliens.hash(state);
            self.houses.hash(state);
        }
    }
}

impl State<CosmicExpressGlobal, ()> for CosmicExpressLocal {
    fn is_valid(&self, g: &CosmicExpressGlobal) -> bool {
        // Check goop, if we have no un-gooped seats and there are non-Green aliens left, it's invalid
        if self.seat_goop.iter().all(|&b| b)
            && self
                .aliens
                .iter()
                .filter(|(_, c)| *c != Color::Green)
                .count()
                > 0
        {
            return false;
        }

        // Flood fill from the current head, stopping at all walls and current path
        // If we can't reach all remaining aliens, houses, and the exit
        if *FLOODFILL_VALIDATOR {
            let max_size = (g.width * g.height) as usize;

            let mut reachable = FxHashSet::default();
            reachable.reserve(max_size);

            let mut to_check = Vec::with_capacity(max_size);
            to_check.push(*self.path.last().unwrap());

            // All points current under a seat are reachable
            // This took a while to track down
            self.path
                .iter()
                .rev()
                .take(1 + self.seats.len())
                .for_each(|p| {
                    reachable.insert(*p);
                });

            // Flood fill from the head of the current path
            while let Some(p) = to_check.pop() {
                reachable.insert(p);

                // Flood fill all empty points
                for neighbor in p.neighbors().into_iter() {
                    // Only check each point once
                    if reachable.contains(&neighbor) || to_check.contains(&neighbor) {
                        continue;
                    }

                    // Don't add points out of bounds (remember there's a border)
                    if neighbor.x < 1
                        || neighbor.y < 1
                        || neighbor.x > g.width
                        || neighbor.y > g.height
                    {
                        continue;
                    }

                    // Keep expanding empty points
                    if !(g.walls.contains(&neighbor) || self.path.contains(&neighbor)) {
                        to_check.push(neighbor);
                    }
                }
            }

            // Expand all points by one: any alien or house *adjacent* to a reachable point is also reachable
            let mut expanded = FxHashSet::default();
            expanded.reserve(max_size);

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
            if !reachable.contains(&g.exit) {
                return false;
            }
        }

        // All validators passed
        true
    }

    fn is_solved(&self, g: &CosmicExpressGlobal) -> bool {
        self.aliens.len() == 0 && self.houses.len() == 0 && self.path.last().unwrap() == &g.exit
    }

    fn next_states(&self, g: &CosmicExpressGlobal) -> Option<Vec<(i64, (), Self)>> {
        let mut result = Vec::with_capacity(4);

        'neighbor: for p in self.path.last().unwrap().neighbors() {
            // Validate that the next point is valid
            // Cannot leave the bounds unless on entrance or exit
            if !(g.entrance == p || g.exit == p)
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
                    for (house_index, (house_point, house_color)) in new_local.houses.iter().enumerate()
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
                if seat_contents.is_none() {
                    // Find all viable aliens
                    // This has to be done this way because if two try to load at once, none get to
                    let mut viable_aliens = Vec::new();
                    for (alien_index, (alien_point, alien_color)) in new_local.aliens.iter().enumerate()
                    {
                        if seat_point.manhattan_distance(*alien_point) == 1 {
                            // Non-green aliens will not try to sit in gooped seats
                            if alien_color != &Color::Green && new_local.seat_goop[seat_index] {
                                continue;
                            }

                            // Record that this alien can be loaded
                            viable_aliens.push((alien_index, *alien_color));
                        }
                    }

                    // If we found exactly one, load it
                    if viable_aliens.len() == 1 {
                        let (alien_index, alien_color) = viable_aliens[0];

                        new_local.aliens.remove(alien_index);
                        new_local.seats[seat_index] = Some(alien_color);

                        if alien_color == Color::Green {
                            new_local.seat_goop[seat_index] = true;
                        }
                    }
                }
            }

            result.push((1, (), new_local));
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

        // Add path characters
        for (i, p) in self.path.iter().enumerate() {
            let p_before = if i == 0 { 
                None 
            } else { 
                self.path.get(i - 1)
            };
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

        // Add entrance and exit(s) (overwrite path ends)
        output_chars[(g.entrance.y * output_width + g.entrance.x) as usize] = '[';
        output_chars[(g.exit.y * output_width + g.exit.x) as usize] = ']';

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

        if *HEURISTIC_NEAREST_HOUSE {
            let current_point = self.path.last().unwrap();

            // Distance from the path to the nearest remaining alien
            if let Some(distance) = self
                .aliens
                .iter()
                .map(|(alien_point, _)| alien_point.manhattan_distance(*current_point))
                .min()
            {
                heuristic += distance as i64;
            }

            // Distance from each alien to the closest matching house
            for (alien_point, alien_color) in self.aliens.iter() {
                let nearest_house = self
                    .houses
                    .iter()
                    .filter(|(_, house_color)| house_color == alien_color)
                    .map(|(house_point, _)| alien_point.manhattan_distance(*house_point))
                    .min();

                if let Some(distance) = nearest_house {
                    heuristic += distance as i64;
                }
            }

            // Also, distance to the exit if we're out of aliens
            heuristic += current_point.manhattan_distance(global.exit) as i64;
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

    assert!(definition.entrances.len() == 1);
    assert!(definition.exits.len() == 1);

    let global = CosmicExpressGlobal {
        width: definition.width,
        height: definition.height,
        length: definition.length,
        entrance: *definition.entrances.first().unwrap(),
        exit: *definition.exits.first().unwrap(),
        walls,
    };

    let seats = (0..definition.length).map(|_| None).collect();
    let seat_goop = (0..definition.length).map(|_| false).collect();

    let aliens: Vec<(Point, Color)> = definition
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

    let houses: Vec<(Point, Color)> = definition
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

    // Validity checks
    {
        // Counts of each color alien and house match
        let mut alien_colors = aliens.iter().map(|(_, c)| c).collect::<Vec<_>>();
        alien_colors.sort();

        let mut house_colors = houses.iter().map(|(_, c)| c).collect::<Vec<_>>();
        house_colors.sort();

        assert_eq!(alien_colors, house_colors, "Alien and house counts don't match: {alien_colors:?} {house_colors:?}");

        // No entities (from the original definition) overlap
        let mut points = FxHashSet::default();
        for (p, _) in definition.entities.iter() {
            assert!(points.insert(*p), "Entities overlap at {p:?}");
        }
        
        // The entrances and exits are all unique
        let mut points = FxHashSet::default();
        for p in definition.entrances.iter().chain(definition.exits.iter()) {
            assert!(points.insert(*p), "Entrance/exit overlap {p:?}");
        }
    }

    let local = CosmicExpressLocal {
        path: vec![global.entrance],
        seats,
        seat_goop,
        aliens,
        houses,
    };

    let mut solver = Solver::new(global.clone(), local.clone());

    if *DEBUG_PRINT {
        println!("{}", solver.stringify(&local));
    }

    while let Some(state) = solver.next() {
        if *DEBUG_PRINT && solver.states_checked() % 100000 == 0 {
            println!("{}", state.stringify(&global));
            println!(
                "{} states, {} invalidated, {} queued, {} seconds",
                solver.states_checked(),
                solver.states_invalidated(),
                solver.in_queue(),
                solver.time_spent()
            );
        }
    }
    let solution = solver.get_solution();

    if let Some(solution) = solution {
        println!("{}", solver.stringify(&solution));
    } else {
        println!("No solution found");
    }

    if *DEBUG_PRINT {
        println!(
            "{} states, {} invalidated, {} seconds",
            solver.states_checked(),
            solver.states_invalidated(),
            solver.time_spent()
        );
    }
}
