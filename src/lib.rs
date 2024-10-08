use core::fmt::Debug;
use fxhash::FxHashMap;
// use priority_queue::PriorityQueue;
use keyed_priority_queue::KeyedPriorityQueue as PriorityQueue;
use std::fmt::Display;
use std::hash::Hash;
use std::time::Instant;

pub mod direction;
pub mod point;

pub use direction::Direction;
pub use point::Point;

/// A trait for a state in a search problem
///
/// `G` should store global state which is shared between all states and does not have to be cloned
/// `S` should be the type of a single 'step' for reporting the final path
pub trait State<G, S>: Clone + Eq + Hash {
    /// Given the current state, return a vector of all possible next states with the cost for each single step
    /// We will call is_valid on all states, so these do not *have* to be valid
    fn next_states(&self, global: &G) -> Option<Vec<(i64, S, Self)>>
    where
        Self: Sized;

    /// Given the current state, return an estimate to the solution from this state
    /// In order to find an optimal solution, this should not overestimate
    fn heuristic(&self, global: &G) -> i64;

    /// Is this state valid? If not, we will not consider it or any state from here
    fn is_valid(&self, global: &G) -> bool;

    /// Is this state a valid solution? If so, this will end the search (given A* it should be optimal)
    fn is_solved(&self, global: &G) -> bool;

    /// Display this state for debugging purposes
    fn stringify(&self, global: &G) -> String;
}

/// Create a solver for the current problem; this will store current states etc
///
/// Global state is shared between all states and shouldn't change
/// Local state is what we're searching through and should only store mutable values
#[derive(Debug)]
pub struct Solver<GlobalState, LocalState: State<GlobalState, Step>, Step> {
    // A link to the shared global state
    global_state: GlobalState,

    // A* parameters
    solution: Option<LocalState>,
    to_check: PriorityQueue<LocalState, i64>,
    steps: FxHashMap<LocalState, (Step, LocalState)>,
    distances: FxHashMap<LocalState, i64>,

    // Debug values
    states_checked: usize,
    states_invalidated: usize,
    states_time_pruned: usize,
    time_spent: f32,
}

/// Implement a generic A* solver given a Global and local State
///
/// Global state is shared between all states and shouldn't change
/// Local state is what we're searching through and should only store mutable values
impl<GlobalState, LocalState: State<GlobalState, Step>, Step: Copy>
    Solver<GlobalState, LocalState, Step>
{
    /// Initialize a new solver with a global and local state
    pub fn new(
        global_state: GlobalState,
        initial_state: LocalState,
    ) -> Solver<GlobalState, LocalState, Step> {
        let mut to_check = PriorityQueue::new();
        to_check.push(initial_state.clone(), 0);

        let steps = FxHashMap::default();

        let mut distances = FxHashMap::default();
        distances.insert(initial_state.clone(), 0);

        Solver {
            global_state,

            solution: None,
            to_check,
            steps,
            distances,

            states_checked: 0,
            states_invalidated: 0,
            states_time_pruned: 0,
            time_spent: 0 as f32,
        }
    }

    /// Return the final solution (if found)
    pub fn get_solution(&self) -> Option<LocalState> {
        self.solution.clone()
    }

    /// Determine how many states we've checked (mostly for debug output)
    pub fn states_checked(&self) -> usize {
        self.states_checked
    }

    pub fn states_invalidated(&self) -> usize {
        self.states_invalidated
    }

    pub fn in_queue(&self) -> usize {
        self.to_check.len()
    }

    /// Determine how long we've been searching (mostly for debug output)
    pub fn time_spent(&self) -> f32 {
        self.time_spent
    }

    /// Display the current state of the solver (mostly for debug output)
    pub fn stringify(&self, state: &LocalState) -> String {
        state.stringify(&self.global_state)
    }

    // Return the path between two states (if one exists)
    pub fn path(&self, src: &LocalState, dst: &LocalState) -> Option<Vec<Step>> {
        let mut path = Vec::new();
        let mut current = dst.clone();

        while current != *src {
            let (step, next) = self.steps.get(&current).unwrap().clone();
            path.push(step);
            current = next;
        }

        path.reverse();
        Some(path)
    }
}

impl<GlobalState, LocalState: State<GlobalState, Step> + Debug, Step> Display
    for Solver<GlobalState, LocalState, Step>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Solver<time={}, checked={}, queue={}, invalid={}, time_pruned={}>",
            self.time_spent,
            self.states_checked,
            self.to_check.len(),
            self.states_invalidated,
            self.states_time_pruned,
        )
    }
}

/// Iterate through the search space for the current state
///
/// This will:
/// - take the next state to investigate
/// - if it's a solution, return it directly
/// - otherwise, find neighbors from this state:
///     - skip any that are invalid
///     - calculate heuristic and priority queue the new state
/// - update debug statistics
/// - yield the state investigated
impl<GlobalState, LocalState: State<GlobalState, Step> + Debug, Step> Iterator
    for Solver<GlobalState, LocalState, Step>
{
    type Item = LocalState;

    /// Return the next state in the search, or None if we've exhausted all states
    fn next(&mut self) -> Option<Self::Item> {
        // We have nothing left to search or we've already returned a valid solution
        if self.solution.is_some() || self.to_check.is_empty() {
            return None;
        }

        // Update debug parameters
        let start = Instant::now();
        self.states_checked += 1;

        // Unwrap the next state to investigate
        let (current_state, _) = self.to_check.pop().unwrap();

        // If this is a valid solution, return it directly and mark solved
        if current_state.is_solved(&self.global_state) {
            self.time_spent += start.elapsed().as_secs_f32();
            self.solution = Some(current_state.clone());

            return Some(current_state.clone());
        }

        // Otherwise, iterate and add neighbors to the queue with priority based on heuristic
        if let Some(ls) = current_state.next_states(&self.global_state) {
            for (step_distance, step, next_state) in ls {
                // Skip invalidate states
                if !next_state.is_valid(&self.global_state) {
                    self.states_invalidated += 1;
                    continue;
                }

                // The estimated score is distance to current + step + heuristic
                let current_distance = self.distances[&current_state];
                let estimated_distance =
                    current_distance + step_distance + next_state.heuristic(&self.global_state);

                let best_next_distance =
                    self.distances.get(&next_state).cloned().unwrap_or(i64::MAX);

                // If we've already found a better path, ignore this one
                // if estimated_distance >= best_next_distance {
                if current_distance + step_distance >= best_next_distance {
                    self.states_time_pruned += 1;
                    continue;
                }

                // Update the distance map with the known distance
                self.distances
                    .insert(next_state.clone(), current_distance + step_distance);

                // Otherwise, record this step and add to queue
                self.steps
                    .insert(next_state.clone(), (step, current_state.clone()));

                // Insert the next state with estimated distance
                // Or if was already queued, updated to lower priority
                // 'Biggest' priority will be popped first, so we negate the estimated distance 
                let old_priority = self.to_check.get_priority(&next_state);
                if old_priority.is_none() || -estimated_distance < *old_priority.unwrap() {
                    self.to_check.push(next_state, -estimated_distance);
                }
            }
        }

        // Now update debugger and return the state we checked
        self.time_spent += start.elapsed().as_secs_f32();
        return Some(current_state.clone());
    }
}
