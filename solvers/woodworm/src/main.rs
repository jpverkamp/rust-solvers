use itertools::Itertools;
use std::io::Read;

use direction::Direction;
use solver::{Solver, State};

mod global;
mod local;
mod step;
mod state;

use global::Global;

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    let global = Global::from(input.as_str());
    let local = global.make_local();

    let tracing_enabled = std::env::var("woodworm_TRACE").is_ok();

    if tracing_enabled {
        tracing_subscriber::fmt()
            .without_time()
            .with_max_level(tracing::Level::DEBUG)
            .init();
    } else {
        env_logger::init();
    }

    log::info!("Initial state:\n{}", local.stringify(&global));

    // If we have args, run each of those as a solution
    if std::env::args().len() > 1 {
        for sequence in std::env::args().skip(1) {
            let mut local = global.make_local();

            for d in sequence.chars() {
                let d = match d {
                    'U' | 'u' => Direction::Up,
                    'D' | 'd' => Direction::Down,
                    'L' | 'l' => Direction::Left,
                    'R' | 'r' => Direction::Right,
                    ' ' => continue,
                    _ => panic!("Invalid direction: {}", d),
                };

                match local.step(d, &global) {
                    Ok(new_state) => {
                        log::debug!("Step: {d:?}");
                        log::debug!("{}", new_state.stringify(&global));
                        log::debug!("{new_state}");
                        local = new_state;
                    }
                    Err(e) => {
                        log::error!("Invalid step in direction {:?}: {}", d, e);
                        break;
                    }
                }
            }
        }

        return;
    }

    // Otherwise, run the solver
    let mut solver = Solver::new(global.clone(), local.clone());

    while let Some(state) = solver.next() {
        if solver.states_checked() % 100_000 == 0 {
            log::debug!("\n{}", state.stringify(&global));
            log::debug!("{solver}");
        }
    }
    let solution = solver.get_solution();

    if let Some(solution) = solution {
        log::debug!("{solver}");
        log::info!("{}", solver.stringify(&solution));

        let path = solver.path(&local, &solution).unwrap();
        log::info!("Path: {:?}", path);

        let sequence = path
            .iter()
            .map(|d| match d {
                Direction::Up => 'U',
                Direction::Down => 'D',
                Direction::Left => 'L',
                Direction::Right => 'R',
            })
            .collect::<String>();

        let sequence = sequence
            .chars()
            .chunks(5)
            .into_iter()
            .map(|chunk| chunk.collect::<String>())
            .collect::<Vec<_>>()
            .join(" ");

        println!("Sequence: {}", sequence);
    } else {
        println!("No solution found");
        std::process::exit(1);
    }
}
