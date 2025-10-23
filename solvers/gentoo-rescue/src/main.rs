mod model;
mod simulation;

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
            .without_time()
            .with_max_level(tracing::Level::DEBUG)
            .init();
    } else {
        env_logger::init();
    }

    log::info!("Initial state:\n{}", map.stringify(&()));

    // If we have args, run each of those as a solution
    if std::env::args().len() > 1 {
        for sequence in std::env::args().skip(1) {
            todo!("{sequence:?}")
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
        print!("{}\t", map.critters[0]);

        for step in path {
            match step {
                simulation::Step::SwitchCritter { critter } => {
                    print!("\n{critter}\t");
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
                        print!("\n{critter}\t");
                    }
                }
            }
        }
        println!();
    } else {
        println!("No solution found");
    }
}
