use std::collections::{HashSet, VecDeque};
use std::hash::Hash;
use std::time::Instant;

pub trait State: Clone + Eq + Hash {
    fn next_states(self) -> Option<Vec<Self>> where Self: Sized + Copy;
    
    fn is_valid(self) -> bool;
    fn is_solved(self) -> bool;
}

#[derive(Debug)]
pub enum SearchMode {
    BreadthFirst,
    DepthFirst,
}

#[derive(Debug)]
pub struct Solver<S: State> {
    to_check: VecDeque<S>,
    checked: HashSet<S>,
    search_mode: SearchMode,
    time_spent: f32,
}

impl<S> Solver<S> where S: State + Copy {
    pub fn new(initial_state: &S) -> Solver<S> {
        Solver {
            to_check: VecDeque::from([initial_state.clone()]),
            checked: HashSet::new(),
            search_mode: SearchMode::DepthFirst,
            time_spent: 0 as f32,
        }
    }

    pub fn set_mode(&mut self, new_mode: SearchMode) {
        self.search_mode = new_mode;
    }

    pub fn states_checked(&self) -> usize {
        self.checked.len()
    }

    pub fn time_spent(&self) -> f32 {
        self.time_spent
    }
}

// Iterate through the solver's given solutions until all are returned
impl<S> Iterator for Solver<S> where S: State + Copy {
    type Item = S;

    fn next(&mut self) -> Option<Self::Item> {
        let start = Instant::now();

        while !self.to_check.is_empty() {
            let current_state = match self.search_mode {
                SearchMode::BreadthFirst => self.to_check.pop_front().unwrap(),
                SearchMode::DepthFirst => self.to_check.pop_back().unwrap()
            };
            self.checked.insert(current_state);
            
            // If this is a valid solution, return it from our iterator
            if current_state.is_solved() {
                self.time_spent += start.elapsed().as_secs_f32();
                return Some(current_state);
            }
            
            // Otherwise:
            // - If we have next states, push those onto the queue
            // - If not, just drop this state (and effectively backtrack)
            if let Some(ls) = current_state.next_states() {
                for next_state in ls {
                    if self.checked.contains(&next_state) { continue; }
                    if !next_state.is_valid() { continue; }
        
                    self.to_check.push_back(next_state);
                }
            }
        }

        // If we make it here, return None to stop iterator
        self.time_spent += start.elapsed().as_secs_f32();
        return None;
    }

}    
