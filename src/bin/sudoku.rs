use std::fmt;
use std::io::{self, BufRead};

use solver::{State, SearchMode, Solver};

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
struct Sudoku {
    board: [[u8; 9]; 9]
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

impl State for Sudoku {
    fn is_valid(self) -> bool {
        // TODO: fix this
        return true;
    }
    
    fn is_solved(self) -> bool {
        for x in 0..9 {
            for y in 0..9 {
                if self.board[x][y] == 0 {
                    return false;
                }
            }
        }
        
        return true;
    }

    fn next_states(self) -> Option<Vec<Sudoku>> {
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
                        next.board[x][y] = value;
                        states.push(next);
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
    
    {
        println!("\nBreadth first:");
        let mut solver = Solver::new(&initial_state);
        solver.set_mode(SearchMode::BreadthFirst);
        match solver.next() {
            Some(solution) => println!("solution: {}", solution),
            None => println!("no solution found! :(")
        };
        println!("{} states, {} seconds", solver.states_checked(), solver.time_spent());
    }

    {
        println!("\nDepth first:");
        let mut solver = Solver::new(&initial_state);
        solver.set_mode(SearchMode::DepthFirst);
        match solver.next() {
            Some(solution) => println!("solution: {}", solution),
            None => println!("no solution found! :(")
        };
        println!("{} states, {} seconds", solver.states_checked(), solver.time_spent());
    }

    {
        println!("\nCheck all (loop):");
        let mut solver = Solver::new(&initial_state);
        
        loop {
            match solver.next() {
                Some(solution) => println!("solution: {}", solution),
                None => { break; }
            };
        }

        println!("{} states, {} seconds", solver.states_checked(), solver.time_spent());
    }

    {
        println!("\nCheck all (for):");
        let mut solver = Solver::new(&initial_state);
        
        for solution in &mut solver {
            println!("solution: {}", solution);
        }

        println!("{} states, {} seconds", solver.states_checked(), solver.time_spent());
    }
    
    
}