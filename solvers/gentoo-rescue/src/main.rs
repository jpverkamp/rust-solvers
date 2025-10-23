mod global;
mod local;
mod state;

use std::io::Read;

use global::Global;
use solver::{Solver, State};

use crate::global::Critter;

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    let global = Global::from(input.as_str());
    let local = global.make_local();

    let tracing_enabled = std::env::var("SOLVER_TRACE").is_ok();
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
            todo!("{sequence:?}")
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

    if let Some(solution) = solver.get_solution() {
        log::info!("{solution:?}");

        let path = solver.path(&local, &solution).unwrap();
        let mut last_index = usize::MAX;

        for (i, d) in path {
            if last_index != i {
                let critter = global.critter(i);
                print!(
                    "\n{}\t{:?}\t{:?}\t",
                    Critter::index_char(i),
                    critter.color,
                    critter.kind,
                );
                last_index = i;
            }
            let c = match d {
                direction::Direction::Up => "U",
                direction::Direction::Down => "D",
                direction::Direction::Left => "L",
                direction::Direction::Right => "R",
            };
            print!("{c}");
        }
        println!();
    } else {
        println!("No solution found");
    }
}
