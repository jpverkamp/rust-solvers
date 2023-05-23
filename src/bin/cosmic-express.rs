use rand::seq::SliceRandom;
use rand::thread_rng;
use std::hash::Hash;
use std::io;

use serde::{Deserialize, Serialize};

use solver::{Solver, State};

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
struct Point {
    x: i32,
    y: i32,
}

impl Point {
    fn neighbors(&self) -> Vec<Point> {
        let mut points = vec![
            Point {
                x: self.x + 1,
                y: self.y,
            },
            Point {
                x: self.x,
                y: self.y + 1,
            },
            Point {
                x: self.x,
                y: self.y - 1,
            },
            Point {
                x: self.x - 1,
                y: self.y,
            },
        ];

        points.shuffle(&mut thread_rng());

        return points;
    }

    fn is_neighbor(&self, other: &Point) -> bool {
        self.x == other.x && (self.y - other.y).abs() == 1
            || self.y == other.y && (self.x - other.x).abs() == 1
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
enum Color {
    Purple,
    Orange,
}

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
enum Entity {
    Wall,
    Alien { color: Color },
    House { color: Color },
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct CosmicExpressGlobal {
    width: i32,
    height: i32,
    length: i32,

    entrances: Vec<Point>,
    exits: Vec<Point>,

    entities: Vec<(Point, Entity)>,
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct Path(Vec<Point>);

impl State<CosmicExpressGlobal> for Path {
    fn is_valid(&self, _g: &CosmicExpressGlobal) -> bool {
        return true;
    }

    fn is_solved(&self, g: &CosmicExpressGlobal) -> bool {
        if self.0.len() == 0 {
            return false;
        }

        // The current alien sitting in each seat
        let mut seats: Vec<Option<Entity>> = Vec::new();
        for _ in 0..=g.length {
            seats.push(None);
        }

        // The aliens that have not been picked up yet
        let mut aliens: Vec<&(Point, Entity)> = g
            .entities
            .iter()
            .filter(|(_, e)| matches!(e, Entity::Alien { .. }))
            .collect();

        // The houses that have not been dropped off too yet
        let mut houses: Vec<&(Point, Entity)> = g
            .entities
            .iter()
            .filter(|(_, e)| matches!(e, Entity::House { .. }))
            .collect();

        // Loop over points on the path
        // Note: Each car is one back on the path (if such a thing is possible)
        for (engine_i, _) in self.0.iter().enumerate() {
            for car_i in 0..g.length {
                // If the car isn't on the track yet, ignore it
                if (engine_i as i32 - car_i as i32) < 1 {
                    continue;
                }

                // Figure out which car we're dealing with
                let car_p = self.0.get(engine_i as usize - car_i as usize - 1).unwrap();
                match seats.get(car_i as usize).unwrap() {
                    // If there's an alien in that car already
                    Some(Entity::Alien { color: alien_color }) => {
                        let mut matches = Vec::new();

                        // Find matching houses
                        for (i, (p, e)) in houses.iter().enumerate() {
                            match e {
                                Entity::House { color: house_color } => {
                                    if car_p.is_neighbor(p) && alien_color == house_color {
                                        matches.push(i);
                                    }
                                }
                                _ => {
                                    panic!(
                                        "Something that isn't a house ({:?}) ended up in houses",
                                        e
                                    );
                                }
                            }
                        }

                        // If we have one, remove from houses
                        // TODO: Multiple houses?
                        if matches.len() > 0 {
                            houses.remove(*matches.get(0).unwrap());
                            seats[car_i as usize] = None;
                        }
                    }
                    // Seat is empty
                    _ => {
                        let mut matches = Vec::new();

                        for (i, (p, _)) in aliens.iter().enumerate() {
                            if car_p.is_neighbor(p) {
                                matches.push(i);
                            }
                        }

                        // Can only load exactly one alien
                        // Put it in the seat and remove it from the aliens map
                        if matches.len() == 1 {
                            let alien_i = *matches.get(0).unwrap();
                            seats[car_i as usize] = Some(aliens[alien_i].1.to_owned());
                            aliens.remove(alien_i);
                        }
                    }
                }
            }
        }

        // If we are at the end of the path, check that we have no more aliens or houses
        // TODO: currently assumes aliens and houses are 1:1
        if aliens.len() > 0 || houses.len() > 0 {
            return false;
        }

        // If we get to the end and are adjacent to an exit, return true
        for exit in g.exits.iter() {
            if exit == self.0.last().unwrap() {
                return true;
            }
        }

        // TODO: this
        return false;
    }

    fn next_states(&self, g: &CosmicExpressGlobal) -> Option<Vec<(i64, Path)>> {
        let mut result = Vec::new();

        // If we don't have a path yet, start with each entrance
        if self.0.len() == 0 {
            for entrance in g.entrances.iter() {
                let mut new_path = Path(Vec::new());
                new_path.0.push(entrance.clone());
                result.push((0, new_path));
            }
        }
        // If we do have a path, expand the last node
        else {
            'neighbor: for p in self.0.last().unwrap().neighbors() {
                // Validate that the next point is valid
                // Cannot leave the bounds unless on entrance or exit
                if !(g.entrances.contains(&p) || g.exits.contains(&p))
                    && (p.x < 1 || p.y < 1 || p.x > g.width || p.y > g.height)
                {
                    continue 'neighbor;
                }

                // Cannot intersect any entity
                for (ep, _) in g.entities.iter() {
                    if &p == ep {
                        continue 'neighbor;
                    }
                }

                // Cannot visit the same tile more than once
                for p2 in self.0.iter() {
                    if &p == p2 {
                        continue 'neighbor;
                    }
                }

                let mut new_path = self.clone();
                new_path.0.push(p);
                result.push((1, new_path));
            }
        }

        // We'll always have nodes, so always return
        // We're relying on is_valid to filter impossible states this time
        Some(result)
    }

    fn display(&self, g: &CosmicExpressGlobal) {
        for y in 1..=g.height {
            for x in 1..=g.width {
                let p = Point { x, y };

                if self.0.contains(&p) {
                    let mut pi = 0;
                    for (i, p2) in self.0.iter().enumerate() {
                        if &p == p2 {
                            pi = i;
                        }
                    }

                    let p_before = self.0.get((pi - 1) as usize);
                    let p_after = self.0.get((pi + 1) as usize);

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
                        print!("{}", c);
                    } else {
                        print!("X");
                    }
                } else {
                    let mut printed = false;

                    for (ep, e) in g.entities.iter() {
                        if p != *ep {
                            continue;
                        }

                        match e {
                            Entity::Alien { color: c } => {
                                print!("{}", ('a' as u8 + *c as u8) as char)
                            }
                            Entity::House { color: c } => {
                                print!("{}", ('A' as u8 + *c as u8) as char)
                            }
                            Entity::Wall => print!("*"),
                        }

                        printed = true;
                    }

                    if !printed {
                        print!(".");
                    }
                }
            }
            println!();
        }
    }

    fn heuristic(&self, global: &CosmicExpressGlobal) -> i64 {
        let p = self.0.last().unwrap();
        let mut min_d = i64::MAX;

        for exit in global.exits.iter() {
            let d = (p.x - exit.x).abs() + (p.y - exit.y).abs();
            min_d = min_d.min(d.into());
        }

        return min_d;
    }
}

fn main() {
    env_logger::init();

    let stdin = io::stdin().lock();
    let global: CosmicExpressGlobal = serde_json::from_reader(stdin).unwrap();
    let initial_path = Path(Vec::new());

    println!("Global State: {:#?}", global);

    let mut solver: Solver<CosmicExpressGlobal, Path> = Solver::new(global.clone(), initial_path.clone());
    solver.display(&initial_path);

    while let Some(state) = solver.next() {
        if solver.states_checked() % 10000 == 0 {
            state.display(&global);
            println!(
                "{} states, {} seconds",
                solver.states_checked(),
                solver.time_spent()
            );
        }
        
    }
    let solution = solver.get_solution();

    println!("{:?}", solution);
    if let Some(path) = solution {
        solver.display(&path);
    }

    println!(
        "{} states, {} seconds",
        solver.states_checked(),
        solver.time_spent()
    );
}
