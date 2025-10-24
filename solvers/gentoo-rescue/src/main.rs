mod model;
mod simulation;

use direction::Direction;
use solver::{Solver, State};
use std::io::Read;

use crate::model::map::Map;

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    let map = Map::from(input.as_str());

    let tracing_enabled = std::env::var("RUST_TRACE").is_ok();
    if tracing_enabled {
        tracing_subscriber::fmt()
            .pretty()
            .without_time()
            .with_max_level(tracing::Level::DEBUG)
            .init();
    } else {
        env_logger::init();
    }

    log::info!("Initial state:\n{}", map.stringify(&()));

    // If we have args, run each of those as a solution
    if std::env::args().len() > 1 {
        for arg in std::env::args().skip(1) {
            let mut map = map.clone();

            for line in arg.lines() {
                if line.is_empty() {
                    continue;
                }

                // All lines should be "{row} {column} {color} {kind} {movements...}"
                let parts: Vec<_> = line.split_ascii_whitespace().collect();
                if parts.len() == 4 {
                    continue; // If we don't move the first critter
                }
                assert_eq!(parts.len(), 5);

                let row = parts[0]
                    .parse::<isize>()
                    .expect("Lines must start with {row} {col}");
                let col = parts[1]
                    .parse::<isize>()
                    .expect("Lines must start with {row} {col}");

                let index = map
                    .critters
                    .iter()
                    .position(|c| c.location.x == col - 1 && c.location.y == row - 1)
                    .expect("No critter at {row} {col}");
                map.active_critter = index;
                println!(
                    "=== Switched to critter {index}: {} ===",
                    map.critters[map.active_critter]
                );

                for c in parts[4].chars() {
                    let d = match c {
                        'U' => Direction::Up,
                        'D' => Direction::Down,
                        'L' => Direction::Left,
                        'R' => Direction::Right,
                        _ => panic!("Unknown movement char {c}"),
                    };
                    println!("=== Moving {d:?} ===");
                    match map.try_move(d) {
                        Some((new_map, _)) => {
                            map = new_map;
                            println!("{}", map.stringify(&()));
                        }
                        None => {
                            panic!("Failed to move");
                        }
                    }
                    println!();
                }
            }
        }

        return;
    }

    // Otherwise, run the solver
    let mut solver = Solver::new((), map.clone());

    while let Some(state) = solver.next() {
        if solver.states_checked() % 100_000 == 0 {
            log::debug!("\n{}", state.stringify(&()));
            log::debug!("{solver}");
        }
    }

    if let Some(solution) = solver.get_solution() {
        log::info!("{solution:?}");

        let path = solver.path(&map, &solution).unwrap();
        print!("{start_critter:30}", start_critter = map.critters[0]);

        for step in path {
            match step {
                simulation::Step::SwitchCritter { critter } => {
                    print!("\n{critter:<30}");
                }
                simulation::Step::Move {
                    direction,
                    new_critter,
                } => {
                    print!(
                        "{}",
                        match direction {
                            direction::Direction::Up => 'U',
                            direction::Direction::Down => 'D',
                            direction::Direction::Left => 'L',
                            direction::Direction::Right => 'R',
                        }
                    );
                    if let Some(critter) = new_critter {
                        print!("\n{critter:30}");
                    }
                }
            }
        }
        println!();
    } else {
        panic!("No solution found");
    }
}
