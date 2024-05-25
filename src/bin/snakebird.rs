use anyhow::{anyhow, Result};
use fxhash::FxHashMap;
use std::{
    io::{BufRead, Read},
    ops::Add,
};

use solver::{Solver, State};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Tile {
    Empty,
    Wall,
    Exit,
    Spike,
}

#[derive(Debug, Clone)]
struct Global {
    width: usize,
    height: usize,
    tiles: FxHashMap<Point, Tile>,
    solutions: Vec<String>,
}

impl Global {
    // Read a global + local from a Readable
    fn read<R: Read + BufRead>(reader: &mut R) -> Result<(Global, Local)> {
        let mut width = 0;
        let mut height = 0;

        let mut tiles = FxHashMap::default();

        let mut snakes = vec![
            vec![], // 0..=9
            vec![], // a..=a
            vec![], // A..=Z
            vec![], // {|}~
        ];

        let mut fruit = Vec::new();

        for (y, line) in reader.lines().enumerate() {
            let y = y as isize;
            let line = line?;

            if line.is_empty() {
                break;
            }

            if width == 0 {
                width = line.len();
            } else if line.len() != width {
                return Err(anyhow!("Inconsistent line length"));
            }

            height += 1;

            for (x, c) in line.chars().enumerate() {
                let x = x as isize;
                let p = Point { x, y };

                match c {
                    // Known static tiles
                    ' ' | '-' => {}
                    '#' => {
                        tiles.insert(p, Tile::Wall);
                    }
                    '=' => {
                        tiles.insert(p, Tile::Exit);
                    }
                    '!' => {
                        tiles.insert(p, Tile::Spike);
                    }
                    // Fruit is stored in the local state
                    '+' => {
                        fruit.push(p);
                    }
                    // Currently known types of snakes stored in local state
                    '0'..='9' => {
                        snakes[0].push((c, p));
                    }
                    'a'..='z' => {
                        snakes[1].push((c, p));
                    }
                    'A'..='Z' => {
                        snakes[2].push((c, p));
                    }
                    '{'..='~' => {
                        snakes[3].push((c, p));
                    }
                    _ => return Err(anyhow!("Invalid character {c} in map")),
                }
            }
        }

        // Remove empty snakes, sort points by character, convert to snake struct
        snakes.retain(|snake| !snake.is_empty());
        let snakes = snakes
            .iter()
            .map(|points| {
                let mut points = points.clone();
                points.sort();
                Snake {
                    head: points.first().unwrap().0,
                    points: points.iter().map(|(_, point)| *point).collect::<Vec<_>>(),
                }
            })
            .collect::<Vec<_>>();
        let local = Local { snakes, fruit };

        // Any remaining lines are solutions
        let mut solutions = Vec::new();
        for line in reader.lines() {
            let line = line?.trim().to_string();
            if line.is_empty() {
                break;
            }
            solutions.push(line);
        }

        Ok((
            Global {
                width,
                height,
                tiles,
                solutions,
            },
            local,
        ))
    }

    // Get the tile at a point
    fn tile(&self, point: Point) -> Tile {
        self.tiles.get(&point).copied().unwrap_or(Tile::Empty)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Point {
    x: isize,
    y: isize,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Snake {
    head: char,
    points: Vec<Point>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Local {
    snakes: Vec<Snake>,
    fruit: Vec<Point>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Direction {
    Up,
    Down,
    Left,
    Right,
}

impl Add<Direction> for Point {
    type Output = Point;

    fn add(self, direction: Direction) -> Point {
        match direction {
            Direction::Up => Point {
                x: self.x,
                y: self.y - 1,
            },
            Direction::Down => Point {
                x: self.x,
                y: self.y + 1,
            },
            Direction::Left => Point {
                x: self.x - 1,
                y: self.y,
            },
            Direction::Right => Point {
                x: self.x + 1,
                y: self.y,
            },
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Step {
    snake: usize,
    direction: Direction,
}

impl Local {
    fn try_move(&mut self, global: &Global, index: usize, direction: Direction) -> bool {
        let head = self.snakes[index].points.first().unwrap().clone();
        let new_head = head + direction;

        // Cannot move into a wall
        if global.tile(new_head) == Tile::Wall {
            return false;
        }

        // Cannot move into a spike
        if global.tile(new_head) == Tile::Spike {
            return false;
        }

        // Cannot move into yourself
        if self.snakes[index].points.contains(&new_head) {
            return false;
        }

        // If the currently moving snake moves onto the exit, check if it can leave
        // Otherwise, we're going to move as normal; potentially eating fruit
        if global.tile(new_head) == Tile::Exit {
            // Check for exit

            // Cannot exit (or even move onto exit) until all fruit is eaten
            if !self.fruit.is_empty() {
                return false;
            }

            // Otherwise, snake literally exits
            // Bye bye snake
            self.snakes.remove(index);
        } else {
            // Update the snake
            self.snakes[index].points.insert(0, new_head);

            // Potentially eat fruit; don't remove tail if fruit was eaten
            if let Some(fruit_index) = self.fruit.iter().position(|fruit| fruit == &new_head) {
                self.fruit.remove(fruit_index);
            } else {
                self.snakes[index].points.pop();
            }
        }

        // Push any other snakes
        // If any snake is pushed into a wall, the whole move is invalid
        if !self.try_push(global, index, direction, new_head) {
            return false;
        }

        // Apply gravity to all snakes; if any snakes fall out fo the world, the whole move is invalid
        if !self.try_gravity(global) {
            return false;
        }

        // If we made it here, the move was successful
        true
    }

    fn try_push(
        &mut self,
        global: &Global,
        index: usize,
        direction: Direction,
        new_head: Point,
    ) -> bool {
        // Attempt to push any snakes in the way
        let mut snake_pushing_indexes = Vec::new();
        let mut snake_pushing_points = vec![new_head];

        'daisies_remain: loop {
            for (other_index, _) in self.snakes.iter().enumerate() {
                // TODO: What if we have a weird loop where we're pushing ourselves?
                if index == other_index {
                    continue;
                }

                // Already pushing this snake
                if snake_pushing_indexes.contains(&other_index) {
                    continue;
                }

                // If any pushing point is in the new snake, it's getting pushed too
                if snake_pushing_points
                    .iter()
                    .any(|p| self.snakes[other_index].points.contains(&p))
                {
                    snake_pushing_indexes.push(other_index);

                    snake_pushing_points.extend(
                        self.snakes[other_index]
                            .points
                            .iter()
                            .map(|p| *p + direction),
                    );

                    continue 'daisies_remain;
                }
            }

            break;
        }

        // No snake pushing points can hit anything (other than the original head)
        // TODO: Can you push a snake into the exit?
        if snake_pushing_points
            .iter()
            .skip(1)
            .any(|p| global.tile(*p) != Tile::Empty)
        {
            return false;
        }

        // If we have snakes to push, move them all
        for snake_index in snake_pushing_indexes.iter() {
            for point in self.snakes[*snake_index].points.iter_mut() {
                *point = *point + direction;
            }
        }

        true
    }

    fn try_gravity(&mut self, global: &Global) -> bool {
        // Apply gravity to all snakes
        'still_falling: loop {
            // Calculate all snakes directly support
            let mut supported_indexes = Vec::new();
            
            'finding_support: loop {
                // If all snakes are supported, we're done with both loops
                if supported_indexes.len() == self.snakes.len() {
                    break 'still_falling;
                }

                for (index, _) in self.snakes.iter().enumerate() {
                    // Skip snakes we've already supported (but only inner loop)
                    if supported_indexes.contains(&index) {
                        continue;
                    }

                    // Supported by a wall
                    if self.snakes[index].points.iter().any(|point| {
                        global.tile(*point + Direction::Down) == Tile::Wall
                    }) {
                        supported_indexes.push(index);
                        continue 'finding_support;
                    }

                    // Supported by fruit? (we can walk across fruit)
                    if self.snakes[index].points.iter().any(|point| {
                        self.fruit.contains(&(*point + Direction::Down))
                    }) {
                        supported_indexes.push(index);
                        continue 'finding_support;
                    }

                    // Supported by another supported snake
                    for other_index in supported_indexes.iter() {
                        if self.snakes[index].points.iter().any(|point| {
                            self.snakes[*other_index].points.contains(&(*point + Direction::Down))
                        }) {
                            supported_indexes.push(index);
                            continue 'finding_support;
                        }
                    }
                }

                // No more supported snakes found
                break 'finding_support;
            }

            // Otherwise, all non supported snakes fall by one
            for index in 0..self.snakes.len() {
                if supported_indexes.contains(&index) {
                    continue;
                }

                for point in self.snakes[index].points.iter_mut() {
                    *point = *point + Direction::Down;
                }
            }

            // If any snakes fall out of the world, the whole move is invalid
            if self.snakes.iter().any(|snake| {
                snake.points.iter().all(|point| point.y > global.height as isize)
            }) {
                return false;
            }
        }

        true
    }
}

#[cfg(test)]
mod local_move_tests {
    use super::*;

    fn test_world_1() -> (Global, Local) {
        let input = "\
----
10--
###-";
        Global::read(&mut std::io::Cursor::new(input)).unwrap()
    }

    fn test_world_2() -> (Global, Local) {
        let input = "\
----
10#-
###-";
        Global::read(&mut std::io::Cursor::new(input)).unwrap()
    }

    fn test_big_world_1() -> (Global, Local) {
        let input = "\
------------
------------
------------
3210--------
######---###";
        Global::read(&mut std::io::Cursor::new(input)).unwrap()
    }

    fn test_push_world_1() -> (Global, Local) {
        let input = "\
--a--
10b--
####-";
        Global::read(&mut std::io::Cursor::new(input)).unwrap()
    }

    #[test]
    fn test_move() {
        let (global, mut local) = test_world_1();

        assert!(local.try_move(&global, 0, Direction::Right));
        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 2, y: 1 }, Point { x: 1, y: 1 }]
        );
    }

    #[test]
    fn test_move_up() {
        let (global, mut local) = test_world_1();

        assert!(local.try_move(&global, 0, Direction::Up));
        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 1, y: 0 }, Point { x: 1, y: 1 }]
        );
    }

    #[test]
    fn test_right_up() {
        let (global, mut local) = test_world_1();

        assert!(local.try_move(&global, 0, Direction::Right));
        assert!(local.try_move(&global, 0, Direction::Up));
        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 2, y: 0 }, Point { x: 2, y: 1 }]
        );
    }

    #[test]
    fn test_not_into_ground() {
        let (global, mut local) = test_world_1();

        assert!(!local.try_move(&global, 0, Direction::Down));
    }

    #[test]
    fn test_not_backtrack() {
        let (global, mut local) = test_world_1();

        assert!(!local.try_move(&global, 0, Direction::Left));
    }

    #[test]
    fn test_not_into_wall() {
        let (global, mut local) = test_world_2();

        assert!(!local.try_move(&global, 0, Direction::Right));
    }

    #[test]
    fn test_fall() {
        let (global, mut local) = test_world_1();

        // Move on platform
        assert!(local.try_move(&global, 0, Direction::Right));

        // Hang off platform
        assert!(local.try_move(&global, 0, Direction::Right));

        // Fall off
        assert!(!local.try_move(&global, 0, Direction::Right));
    }

    #[test]
    fn test_move_when_upright_fall() {
        let (global, mut local) = test_world_1();

        // Stand up
        assert!(local.try_move(&global, 0, Direction::Up));

        // Move over
        assert!(local.try_move(&global, 0, Direction::Right));

        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 2, y: 1 }, Point { x: 1, y: 1 }]
        );
    }

    #[test]
    fn test_turn_around() {
        let (global, mut local) = test_world_1();

        // Stand up
        assert!(local.try_move(&global, 0, Direction::Up));

        // Turn and fall back left
        assert!(local.try_move(&global, 0, Direction::Left));

        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 0, y: 1 }, Point { x: 1, y: 1 }]
        );
    }

    #[test]
    fn test_big_curl_up() {
        let (global, mut local) = test_big_world_1();

        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Left));

        assert_eq!(
            local.snakes[0].points,
            vec![
                Point { x: 2, y: 2 },
                Point { x: 3, y: 2 },
                Point { x: 3, y: 3 },
                Point { x: 2, y: 3 },
            ]
        );
    }

    #[test]
    fn test_big_curl_up_too_much() {
        let (global, mut local) = test_big_world_1();

        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Left));
        assert!(!local.try_move(&global, 0, Direction::Down));
    }

    #[test]
    fn test_big_cross_chasm() {
        let (global, mut local) = test_big_world_1();

        for _ in 0..8 {
            assert!(local.try_move(&global, 0, Direction::Right));
        }

        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Left));

        assert_eq!(
            local.snakes[0].points,
            vec![
                Point { x: 10, y: 2 },
                Point { x: 11, y: 2 },
                Point { x: 11, y: 3 },
                Point { x: 10, y: 3 },
            ]
        );
    }

    #[test]
    fn test_big_stand() {
        let (global, mut local) = test_big_world_1();

        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Up));

        assert_eq!(
            local.snakes[0].points,
            vec![
                Point { x: 3, y: 0 },
                Point { x: 3, y: 1 },
                Point { x: 3, y: 2 },
                Point { x: 3, y: 3 }
            ]
        );
    }

    #[test]
    fn test_big_stand_and_fall() {
        let (global, mut local) = test_big_world_1();

        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Up));

        assert!(local.try_move(&global, 0, Direction::Right));
        assert_eq!(local.snakes[0].points[0], Point { x: 4, y: 1 });

        assert!(local.try_move(&global, 0, Direction::Right));
        assert_eq!(local.snakes[0].points[0], Point { x: 5, y: 2 });

        assert!(local.try_move(&global, 0, Direction::Right));
        assert_eq!(local.snakes[0].points[0], Point { x: 6, y: 3 });

        assert_eq!(
            local.snakes[0].points,
            vec![
                Point { x: 6, y: 3 },
                Point { x: 5, y: 3 },
                Point { x: 4, y: 3 },
                Point { x: 3, y: 3 }
            ]
        );
    }

    #[test]
    fn test_push() {
        let (global, mut local) = test_push_world_1();

        assert!(local.try_move(&global, 0, Direction::Right));
        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 2, y: 1 }, Point { x: 1, y: 1 }]
        );
        assert_eq!(
            local.snakes[1].points,
            vec![Point { x: 3, y: 0 }, Point { x: 3, y: 1 }]
        );
    }

    #[test]
    fn test_push_off() {
        let (global, mut local) = test_push_world_1();

        assert!(local.try_move(&global, 0, Direction::Right));
        assert!(!local.try_move(&global, 0, Direction::Right));
    }
}

impl State<Global, Step> for Local {
    fn is_valid(&self, global: &Global) -> bool {
        self.snakes.iter().all(|snake|
        snake.points.iter().all(|point| point.y <= global.height as isize))
    }

    fn is_solved(&self, _map: &Global) -> bool {
        // All snakes exited and fruit eaten
        self.snakes.is_empty() && self.fruit.is_empty()
    }

    fn next_states(&self, global: &Global) -> Option<Vec<(i64, Step, Local)>> {
        let mut next_states = Vec::new();

        // Try to move each snake in each direction
        for (snake_index, _) in self.snakes.iter().enumerate() {
            for direction in [
                Direction::Up,
                Direction::Down,
                Direction::Left,
                Direction::Right,
            ]
            .iter()
            {
                // Generate a potential new local state
                let step = Step {
                    snake: snake_index,
                    direction: *direction,
                };
                let mut new_local = self.clone();

                // Try to move the snake
                if !new_local.try_move(global, snake_index, *direction) {
                    continue;
                }

                // If the move was successful, add it to the list of next states
                next_states.push((1, step, new_local));
            }
        }

        if next_states.is_empty() {
            return None;
        }
        Some(next_states)
    }

    fn heuristic(&self, _global: &Global) -> i64 {
        0
    }

    fn stringify(&self, global: &Global) -> String {
        let mut output = String::new();

        for y in 0..global.height {
            let y = y as isize;
            for x in 0..global.width {
                let x = x as isize;
                let mut c = match global.tile(Point { x, y }) {
                    Tile::Empty => ' ',
                    Tile::Wall => '█',
                    Tile::Exit => '⊛',
                    Tile::Spike => '✴',
                };

                for (i, snake) in self.snakes.iter().enumerate() {
                    for (j, pt) in snake.points.iter().enumerate() {
                        if pt.x == x && pt.y == y {
                            c = match i {
                                0 => ('0' as u8 + j as u8) as char,
                                1 => ('a' as u8 + j as u8) as char,
                                2 => ('A' as u8 + j as u8) as char,
                                3 => ('{' as u8 + j as u8) as char,
                                _ => unreachable!(),
                            };
                        }
                    }
                }

                if self.fruit.contains(&Point { x, y }) {
                    c = '';
                }

                output.push(c);
            }
            output.push('\n');
        }

        output
    }
}

fn solve(global: Global, local: Local) -> Option<(Solver<Global, Local, Step>, Local)> {
    let mut solver = Solver::new(global.clone(), local);
    while let Some(state) = solver.next() {
        if solver.states_checked() % 100000 != 0 {
            continue;
        }
        log::info!("{solver}, state:\n{}", state.stringify(&global));
    }

    let solution = solver.get_solution();
    if solution.is_none() {
        log::error!(
            "No solution found after {} states in {} seconds",
            solver.states_checked(),
            solver.time_spent(),
        );
        return None;
    }
    let solution = solution.unwrap();

    log::info!(
        "Solved after {} states in {} seconds:\n{}",
        solver.states_checked(),
        solver.time_spent(),
        solution.stringify(&global),
    );

    Some((solver, solution))
}

fn stringify_solution(solver: &Solver<Global, Local, Step>, initial_state: &Local, solved_state: &Local) -> String {
    let mut last_moved_snake = '\0';

    let mut path = String::new();
    for step in solver.path(initial_state, solved_state).unwrap().iter() {
        let snake = match step.snake {
            0 => '0',
            1 => 'a',
            2 => 'A',
            3 => '{',
            _ => unreachable!(),
        };
        if snake != last_moved_snake {
            path.push(snake);
            last_moved_snake = snake;
        }

        path.push(match step.direction {
            Direction::Up => '↑',
            Direction::Down => '↓',
            Direction::Left => '←',
            Direction::Right => '→',
        });
    }
    path
}

fn main() -> Result<()> {
    env_logger::init();

    let stdin = std::io::stdin();
    let mut reader = stdin.lock();

    let (global, local) = Global::read(&mut reader).unwrap();

    log::info!("Initial state:\n{}", local.stringify(&global));

    // If there is an arg, assume it's a test case and run it
    if std::env::args().len() > 1 {
        std::env::args().skip(1).for_each(|instructions| {
            let mut local = local.clone();
            let mut current_snake_index = 0;

            for (step, c) in instructions.chars().enumerate() {
                println!("\n=== Step {}: {} ===", step, c);

                if "0aA{".contains(c) {
                    current_snake_index = local.snakes.iter().position(|snake| snake.head == c).unwrap();
                    println!("Switched to snake {} ({})", c, local.snakes[current_snake_index].head);
                } else {
                    let direction = match c {
                        '↑' => Direction::Up,
                        '↓' => Direction::Down,
                        '←' => Direction::Left,
                        '→' => Direction::Right,
                        _ => panic!("invalid instruction: {}", c),
                    };

                    println!("Moving snake {} ({}) {:?}", local.snakes[current_snake_index].head, c, direction);
                    assert!(local.try_move(&global, current_snake_index, direction));
                    println!("After move:\n{}", local.stringify(&global));
                    println!("Is valid? {}", local.is_valid(&global));
                    println!("Is solved? {}", local.is_solved(&global));
                }
            }
        });
        return Ok(());
    }

    // Otherwise, try to find a new solution
    if let Some((solver, solution)) = solve(global.clone(), local.clone()) {
        let path = stringify_solution(&solver, &local, &solution);
        println!("{path}");

        // Check against known solutions
        if !global.solutions.iter().any(|s| s == &path) {
            log::warn!("Solution does not match known solution")
        }

        return Ok(());
    } 
    
    Err(anyhow!("No solution found"))
}


#[cfg(test)]
mod test_solutions {
    use std::{fs::File, io::Read, sync::mpsc, thread};
    use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

    use super::*;

    #[test]
    fn test_all_solutions() {
        // Timeout after 1 second or SNAKEBIRD_TEST_TIMEOUT if set
        let timeout = std::time::Duration::from_secs(
            std::env::var("SNAKEBIRD_TEST_TIMEOUT")
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(1),
        );

        // Collect all tests to run in order
        let mut test_files = Vec::new();
        for entry in std::fs::read_dir("data/snakebird").unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if !path.is_dir() {
                continue;
            }

            for entry in std::fs::read_dir(&path).unwrap() {
                let entry = entry.unwrap();
                let path = entry.path();
                if path.extension().unwrap() != "txt" {
                    continue;
                }

                test_files.push(path);
            }
        }
        test_files.sort();
        println!("Running {} tests", test_files.len());

        // Run each test with timeout
        #[derive(Debug, Clone, Eq, PartialEq)]
        enum TestResult {
            Success,
            NoSolution,
            InvalidSolution(String),
            TimedOut,
        }

        let results = test_files
            .par_iter()
            .map(move |path| {
                let mut file = File::open(&path).unwrap();
                let mut input = String::new();
                file.read_to_string(&mut input).unwrap();

                let (global, local) = Global::read(&mut input.as_bytes()).unwrap();

                let (tx, rx) = mpsc::channel();

                let solver_global = global.clone();
                let solver_local = local.clone();

                thread::spawn(move || {
                    let solution = solve(solver_global, solver_local);
                    match tx.send(solution) {
                        Ok(_) => {}
                        Err(_) => {}
                    } // I don't actually care if this succeeds, but need to consume it
                });

                match rx.recv_timeout(timeout) {
                    Ok(solution) => {
                        if solution.is_none() {
                            log::debug!("No solution: {:?}", path);
                            return TestResult::NoSolution;
                        }

                        let (solver, solution) = solution.unwrap();
                        let path = stringify_solution(&solver, &local, &solution);

                        if !global.solutions.contains(&path) {
                            log::debug!("Invalid solution: {:?}", path);
                            return TestResult::InvalidSolution(path);
                        }

                        log::debug!("Solved: {:?}", path);
                        return TestResult::Success;
                    }
                    Err(_) => {
                        log::debug!("Timed out: {:?}", path);
                        return TestResult::TimedOut;
                    }
                }
            })
            .collect::<Vec<_>>();

        // Print out the results
        if results.iter().any(|r| *r == TestResult::TimedOut) {
            println!("\nTimed out tests:");
            for (path, result) in test_files.iter().zip(results.iter()) {
                if *result == TestResult::TimedOut {
                    println!("  {:?}", path);
                }
            }
        }

        if results.iter().any(|r| *r == TestResult::NoSolution) {
            println!("\nUnsolved tests:");
            for (path, result) in test_files.iter().zip(results.iter()) {
                if *result == TestResult::NoSolution {
                    println!("  {:?}", path);
                }
            }
        }

        if results.iter().any(|r| {
            if let TestResult::InvalidSolution(_) = r {
                true
            } else {
                false
            }
        }) {
            println!("\nFailed tests:");
            for (path, result) in test_files.iter().zip(results.iter()) {
                if let TestResult::InvalidSolution(solution) = result {
                    println!("  {:?} -> {:?}", path, solution);
                }
            }
        }

        let perfect = results
            .iter()
            .all(|r| *r == TestResult::Success || *r == TestResult::TimedOut);
        if !perfect {
            println!();
        }
        assert!(perfect);
    }
}
