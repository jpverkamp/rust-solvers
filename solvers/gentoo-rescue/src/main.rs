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

                tracing::debug!("Processing line: {line}");

                // All lines should be "{row} {column} ... {movements...}"
                // The middle could have color, kind, (optional carrying)
                // Or you can leave it out
                let parts: Vec<_> = line.split_ascii_whitespace().collect();

                if parts[0].eq_ignore_ascii_case("swap") {
                    let swap_index = parts[1]
                        .parse::<usize>()
                        .expect("Swap lines should be 'swap {index}'");
                    let swap = &map.swaps[swap_index];
                    log::info!("Using swap {swap_index}: {swap:?}");
                    map.try_swap(swap_index);
                } else {
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

                    for c in parts.iter().last().unwrap().chars() {
                        let d = match c {
                            'U' => Direction::Up,
                            'D' => Direction::Down,
                            'L' => Direction::Left,
                            'R' => Direction::Right,
                            _ => panic!("Unknown movement char {c}"),
                        };
                        log::info!("Moving {} {d:?}", map.critters[index]);

                        if !map.try_move(index, d, true) {
                            panic!("Failed to move");
                        }
                    }
                }

                println!("{}", map.stringify(&()));
                for critter in map.critters.iter() {
                    println!("{critter}");
                }
                println!();
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
                    if critter_index != last_critter_index {
                        if last_critter_index != usize::MAX {
                            println!();
                        }
                        print!("{critter:30}", critter = map.critters[critter_index]);
                        last_critter_index = critter_index;
                    }

                    if !map.try_move(critter_index, direction, true) {
                        panic!("Simulation failed");
                    }

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
                simulation::Step::Swap(swap_index) => {
                    // Find the critter at the swap point
                    let critter_index = map
                        .critters
                        .iter()
                        .position(|c| c.location() == map.swaps[swap_index].critter.location())
                        .expect("No critter at swap point");

                    let previous_critter = map.critters[critter_index];

                    if !map.try_swap(swap_index) {
                        panic!("Simulation failed");
                    }

                    println!();
                    println!(
                        "Swap {swap_index} ({previous_critter} -> {critter})",
                        critter = map.critters[critter_index]
                    );
                    last_critter_index = usize::MAX;
                }
            }
        }
        println!();
    } else {
        panic!("No solution found");
    }
}
