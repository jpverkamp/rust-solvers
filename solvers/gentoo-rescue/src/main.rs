mod model;
mod simulation;

use direction::Direction;
use solver::{Solver, State};
use std::io::Read;

use crate::model::map::Map;

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    let original_map = Map::from(input.as_str());

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

    log::info!("Initial state:\n{}", original_map.stringify(&()));

    // If we have args, run each of those as a solution
    if std::env::args().len() > 1 {
        for arg in std::env::args().skip(1) {
            let mut map = original_map.clone();

            for line in arg.lines() {
                if line.is_empty() || line.starts_with('#') {
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
                    .position(|c| c.location() == (col - 1, row - 1).into())
                    .expect("No critter at {row} {col}");

                log::info!("Switching to critter {index}: {}", map.critters[index]);

                for c in parts[4].chars() {
                    let d = match c {
                        'U' => Direction::Up,
                        'D' => Direction::Down,
                        'L' => Direction::Left,
                        'R' => Direction::Right,
                        _ => panic!("Unknown movement char {c}"),
                    };
                    log::info!("Moving {} {d:?}", map.critters[index]);

                    // This clone is done to reset the teleporter state
                    let mut new_map = map.clone();
                    if new_map.try_move(index, d, true) {
                        map = new_map;
                        println!("{}", map.stringify(&()));
                        for critter in map.critters.iter() {
                            println!("{critter}");
                        }
                    } else {
                        panic!("Failed to move");
                    }
                    println!();
                }
            }
        }

        return;
    }

    // Otherwise, run the solver
    let mut solver = Solver::new((), original_map.clone());

    while let Some(state) = solver.next() {
        if solver.states_checked() % 100_000 == 0 {
            log::debug!("\n{}", state.stringify(&()));
            log::debug!("{solver}");
        }
    }

    if let Some(solution) = solver.get_solution() {
        log::info!("{solution:?}");

        let path = solver.path(&original_map, &solution).unwrap();
        let mut last_critter_index = usize::MAX;
        let mut map = original_map.clone();

        for step in path {
            match step {
                simulation::Step::Move {
                    critter_index,
                    direction,
                } => {
                    // println!(" PRE {step:?}, {c}", c = map.critters[critter_index]);
                    // dbg!(&step, &map);

                    if critter_index != last_critter_index {
                        if last_critter_index != usize::MAX {
                            println!();
                        }
                        print!("{critter:30}", critter = map.critters[critter_index]);
                        last_critter_index = critter_index;
                    }

                    // This clone is done to reset the teleporter state
                    let mut new_map = map.clone();
                    if new_map.try_move(critter_index, direction, true) {
                        map = new_map;
                    } else {
                        panic!("Simulation failed");
                    }

                    // println!("POST {step:?}, {c}", c = map.critters[critter_index]);

                    print!(
                        "{}",
                        match direction {
                            direction::Direction::Up => 'U',
                            direction::Direction::Down => 'D',
                            direction::Direction::Left => 'L',
                            direction::Direction::Right => 'R',
                        }
                    );
                }
            }
        }
        println!();
    } else {
        panic!("No solution found");
    }
}
