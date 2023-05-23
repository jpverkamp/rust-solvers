use core::fmt::Debug;
use std::collections::{HashSet, VecDeque};
use std::hash::Hash;
use std::time::Instant;

pub trait State<G>: Clone + Eq + Hash {
    fn next_states(&self, global: &G) -> Option<Vec<Self>>
    where
        Self: Sized;

    fn is_valid(&self, global: &G) -> bool;
    fn is_solved(&self, global: &G) -> bool;

    fn display(&self, global: &G);
}

#[derive(Debug)]
pub enum SearchMode {
    BreadthFirst,
    DepthFirst,
}

#[derive(Debug)]
pub struct Solver<G, S: State<G>> {
    global_state: G,
    to_check: VecDeque<S>,
    checked: HashSet<S>,
    search_mode: SearchMode,
    time_spent: f32,
}

impl<G, S> Solver<G, S>
where
    S: State<G>,
{
    pub fn new(global_state: G, initial_state: &S) -> Solver<G, S> {
        Solver {
            global_state: global_state,
            to_check: VecDeque::from([initial_state.clone()]),
            checked: HashSet::new(),
            search_mode: SearchMode::DepthFirst,
            time_spent: 0 as f32,
        }
    }

    pub fn set_mode(&mut self, new_mode: SearchMode) -> &mut Self {
        self.search_mode = new_mode;
        self
    }

    pub fn states_checked(&self) -> usize {
        self.checked.len()
    }

    pub fn time_spent(&self) -> f32 {
        self.time_spent
    }

    pub fn display(&self, state: &S)
    where
        S: State<G>,
    {
        state.display(&self.global_state)
    }
}

// Iterate through the solver's given solutions until all are returned
impl<G, S> Iterator for Solver<G, S>
where
    S: State<G> + Debug,
{
    type Item = S;

    fn next(&mut self) -> Option<Self::Item> {
        let start = Instant::now();
        let mut iter = 0;

        while !self.to_check.is_empty() {
            iter += 1;

            let current_state = &match self.search_mode {
                SearchMode::BreadthFirst => self.to_check.pop_front().unwrap(),
                SearchMode::DepthFirst => self.to_check.pop_back().unwrap(),
            };
            // self.checked.insert(current_state.clone());

            // DEBUG
            if iter % 1000 == 0 {
                tracing::debug!(
                    "\n[DEBUG] iter: {}, queue: {}, checked: {}, time: {}",
                    iter,
                    self.to_check.len(),
                    self.checked.len(),
                    start.elapsed().as_secs_f32()
                );
                // self.display(current_state);
            }
            // /DEBUG

            // If this is a valid solution, return it from our iterator
            if current_state.is_solved(&self.global_state) {
                self.time_spent += start.elapsed().as_secs_f32();
                return Some(current_state.clone());
            }

            // Otherwise:
            // - If we have next states, push those onto the queue
            // - If not, just drop this state (and effectively backtrack)
            if let Some(ls) = current_state.next_states(&self.global_state) {
                for next_state in ls {
                    if self.checked.contains(&next_state) {
                        continue;
                    }
                    if !next_state.is_valid(&self.global_state) {
                        continue;
                    }

                    self.to_check.push_back(next_state);
                }
            }
        }

        // If we make it here, return None to stop iterator
        self.time_spent += start.elapsed().as_secs_f32();
        return None;
    }
}
