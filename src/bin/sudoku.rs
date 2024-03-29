use std::fmt;
use std::io::{self, BufRead};

use solver::{Solver, State};

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
struct Play {
    x: u8,
    y: u8,
    value: u8,
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
struct Sudoku {
    board: [[u8; 9]; 9],
}

impl fmt::Display for Sudoku {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s: String = String::new();

        s += "sudoku <\n";
        for x in 0..9 {
            s += "  ";
            for y in 0..9 {
                if self.board[x][y] == 0 {
                    s += " ";
                } else {
                    // This seems weird
                    s = format!("{}{}", s, self.board[x][y]);
                }

                if y == 2 || y == 5 {
                    s += "|";
                }
            }
            s += "\n";

            if x == 2 || x == 5 {
                s += "  ---+---+---\n";
            }
        }
        s += ">";

        write!(f, "{}", s)
    }
}

impl<G> State<G, Play> for Sudoku {
    fn is_valid(&self, _: &G) -> bool {
        // TODO: fix this
        return true;
    }

    fn is_solved(&self, _: &G) -> bool {
        for x in 0..9 {
            for y in 0..9 {
                if self.board[x][y] == 0 {
                    return false;
                }
            }
        }

        return true;
    }

    fn next_states(&self, _: &G) -> Option<Vec<(i64, Play, Sudoku)>> {
        let mut states = Vec::new();

        // Find the next empty square
        for x in 0..9 {
            for y in 0..9 {
                if self.board[x][y] == 0 {
                    // Try each value
                    'duplicate: for value in 1..=9 {
                        for other in 0..9 {
                            // Already used in this row or column, skip
                            if self.board[x][other] == value {
                                continue 'duplicate;
                            }

                            if self.board[other][y] == value {
                                continue 'duplicate;
                            }

                            // Already used in this 3x3 block, skip
                            if self.board[x / 3 * 3 + other / 3][y / 3 * 3 + other % 3] == value {
                                continue 'duplicate;
                            }
                        }

                        // let step = Step { x: x as u8, y: y as u8, value: value };

                        // Valid so far, so generate a new board using that value
                        let mut next = self.clone();
                        let play = Play { x: x as u8, y: y as u8, value: value };
                        next.board[x][y] = value;
                        states.push((1, play, next));
                    }

                    // Return the possible next states for this empty square
                    // If we didn't generate any states, something went wrong
                    if states.len() == 0 {
                        return None;
                    } else {
                        return Some(states);
                    }
                }
            }
        }

        // If we made it here, there are no empty squares, is_solved should be true
        // (We shouldn't make it here)
        None
    }

    fn heuristic(&self, _global: &G) -> i64 {
        let mut empty = 0;

        for x in 0..9 {
            for y in 0..9 {
                if self.board[x][y] == 0 {
                    empty += 1;
                }
            }
        }

        return empty;
    }

    fn display(&self, _global: &G) {
        todo!()
    }
}

fn main() {
    env_logger::init();

    let mut board = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0],
    ];

    for (i, line) in io::stdin().lock().lines().enumerate() {
        for (j, c) in line.unwrap().as_bytes().iter().enumerate() {
            board[i][j] = c - b'0';
        }
    }

    let initial_state = Sudoku { board };
    println!("initial: {}", initial_state);

    let mut solver = Solver::new((), initial_state);

    while let Some(_state) = solver.next() {
        // println!("state: {}", _state);
        // println!(
        //     "{} states, {} seconds",
        //     solver.states_checked(),
        //     solver.time_spent()
        // );
    }

    println!("solution: {}", solver.get_solution().unwrap());

    println!(
        "{} states, {} seconds",
        solver.states_checked(),
        solver.time_spent()
    );
}
