use anyhow::{anyhow, Result};
use fxhash::FxHashMap;
use std::{
    io::{BufRead, Read},
    ops::{Add, Mul, Sub},
};

use solver::{Solver, State};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Point {
    x: isize,
    y: isize,
}

impl Point {
    fn manhattan_distance(&self, other: Point) -> isize {
        (self.x - other.x).abs() + (self.y - other.y).abs()
    }
}

impl Add<Point> for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        Point {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Sub<Point> for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        Point {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Mul<isize> for Point {
    type Output = Point;

    fn mul(self, other: isize) -> Point {
        Point {
            x: self.x * other,
            y: self.y * other,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Direction {
    Up,
    Down,
    Left,
    Right,
}

impl From<Direction> for Point {
    fn from(direction: Direction) -> Point {
        match direction {
            Direction::Up => Point { x: 0, y: -1 },
            Direction::Down => Point { x: 0, y: 1 },
            Direction::Left => Point { x: -1, y: 0 },
            Direction::Right => Point { x: 1, y: 0 },
        }
    }

}

impl Direction {
    fn all() -> Vec<Direction> {
        vec![Direction::Up, Direction::Down, Direction::Left, Direction::Right]
    }

    fn flip(&self) -> Direction {
        match self {
            Direction::Up => Direction::Down,
            Direction::Down => Direction::Up,
            Direction::Left => Direction::Right,
            Direction::Right => Direction::Left,
        }
    }
}

impl TryFrom<&str> for Direction {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        match value {
            "up" => Ok(Direction::Up),
            "down" => Ok(Direction::Down),
            "left" => Ok(Direction::Left),
            "right" => Ok(Direction::Right),
            _ => Err(anyhow!("Invalid direction: {value}")),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Tile {
    Empty,
    Flat(usize),
    Slope(usize, Direction),
    Water(usize),
    Sand(usize),
}

#[derive(Debug, Clone)]
struct Global {
    width: usize,
    height: usize,
    tiles: Vec<Tile>,
    flag: Point,
    solutions: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Card {
    Move(usize),
    Jump(usize),
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct Local {
    ball: Point,
    cards: Vec<Card>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Step {
    card: Card,
    direction: Direction,
}

impl Global {
    // Read a global + local from a Readable
    fn read<R: Read + BufRead>(reader: &mut R) -> Result<(Global, Local)> {
        // First line is size of map, {width}x{height}
        let mut lines = reader.lines();

        let size = lines
            .next()
            .expect("Missing size line")
            .expect("Missing size line")
            .split('x')
            .map(|v| v.parse().expect("Invalid size"))
            .collect::<Vec<usize>>();
        assert!(size.len() == 2, "Invalid size line");

        dbg!(&size);

        let width: usize = size[0];
        let height: usize = size[1];
        
        let mut tiles = vec![Tile::Empty; width * height];
        let mut flag = None;
        let mut ball = None;
        let mut solutions = Vec::new();

        // Read tiles (plus flag and ball)
        for y in 0..height {
            let mut line = lines
                .next()
                .expect("Missing tile line")
                .expect("Missing tile line");
            
            for (x, mut definition) in line.split_ascii_whitespace().enumerate() {
                // Empty space, already handled
                if definition == "x" {
                    continue;
                }

                let mut part;

                let mut height = 0;
                let mut is_flag = false;
                let mut is_ball = false;
                let mut is_slope = false;
                let mut slope_direction = Direction::Up;

                while !definition.is_empty() {
                    if definition.contains(',') {
                        (part, definition) = definition.split_once(',').unwrap();
                    } else {
                        part = definition;
                        definition = "";
                    }

                    // Contains a height definition
                    if let Ok(new_height) = part.parse::<usize>() {
                        height = new_height;
                        continue;
                    }

                    // Contains a flag ... flag
                    if part == "flag" {
                        is_flag = true;
                        continue;
                    }

                    // Contains the starting ball
                    if part == "ball" {
                        is_ball = true;
                        continue;
                    }

                    // Contains a slope
                    if let Ok(new_direction) = Direction::try_from(part) {
                        is_slope = true;
                        slope_direction = new_direction;
                        continue;
                    }

                    // Not something we know yet
                    return Err(anyhow!("Unknown tile definition `{part}` in `{line}`"));
                }
            
                if is_slope {
                    tiles[y * width + x] = Tile::Slope(height, slope_direction);
                } else {
                    tiles[y * width + x] = Tile::Flat(height);
                }

                if is_flag {
                    assert!(flag.is_none(), "Multiple flags in map");
                    flag = Some(Point { x: x as isize, y: y as isize });
                }

                if is_ball {
                    assert!(ball.is_none(), "Multiple balls in map");
                    ball = Some(Point { x: x as isize, y: y as isize });
                }
            }
        }

        // Read an empty line before cards
        lines.next();

        // Now read cards
        let mut cards = Vec::new();
        let line = lines
            .next()
            .expect("Missing card line")
            .expect("Missing card line");

        for card in line.split_ascii_whitespace() {
            assert!(card.len() == 2, "Invalid card: {card}");

            let count = card[0..1].parse().expect("Invalid card count");
            cards.push(match card.chars().nth(1).expect("No card type") {
                '-' => Card::Move(count),
                '/' => Card::Jump(count),
                _ => panic!("Invalid card type: {card}"),
            });
        }

        Ok((
            Global {
                width,
                height,
                tiles,
                flag: flag.expect("No flag in map"),
                solutions,
            },
            Local {
                ball: ball.expect("No ball in map"),
                cards,
            },
        ))
    }
}

impl Global {
    fn tile_at(&self, point: Point) -> Tile {
        // Points out of bounds are always empty
        if point.x < 0 || point.y < 0 || point.x >= self.width as isize || point.y >= self.height as isize {
            return Tile::Empty;
        }

        self.tiles[point.y as usize * self.width + point.x as usize]
    }
}

impl Local {
    fn try_move(&mut self, global: &Global, direction: Direction, strength: usize) -> bool {
        // No more moving to do, we're done
        if strength == 0 {
            return true;
        }

        let current_height = match global.tile_at(self.ball) {
            Tile::Empty => unreachable!(),
            Tile::Flat(height) => height,
            Tile::Slope(height, _) => height,
            Tile::Water(_) => todo!(),
            Tile::Sand(_) => todo!(),
        };
        

        let next_point = self.ball + Point::from(direction);
        let next_tile = global.tile_at(next_point);

        // Trying to move into empty space/out of bounds
        if let Tile::Empty = next_tile {
            return false;
        }

        // TODO: Slopes
        if let Tile::Slope(_, slope_direction) = next_tile {
            todo!();
        }

        // Normal flat tile
        if let Tile::Flat(height) = next_tile {
            // On the same level, just move
            if height == current_height {
                self.ball = next_point;
                return self.try_move(global, direction, strength - 1);
            }

            // New tile is higher, bounce
            // This effectively reverses direction and moves 'back' to the same tile
            if height > current_height {
                return self.try_move(global, direction.flip(), strength - 1);
            }

            // New tile is lower, fall?
            if height < current_height {
                todo!()
            }

            unreachable!();
        }
        
        // Normal flat tile, recur
        self.ball = self.ball + direction.into();
        self.try_move(global, direction, strength - 1)
    }

    fn try_jump(&mut self, global: &Global, direction: Direction, strength: usize) -> bool {
        let next_point = self.ball + Point::from(direction) * strength as isize;
        let next_tile = global.tile_at(next_point);

        // Trying to jump into empty space
        if let Tile::Empty = next_tile {
            return false;
        }

        // Otherwise, it always works
        self.ball = next_point;
        true
    }
}

impl State<Global, Step> for Local {
    fn next_states(&self, global: &Global) -> Option<Vec<(i64, Step, Self)>>
    where
        Self: Sized {
        
        let mut next_states = Vec::new();

        for (i, card) in self.cards.iter().enumerate() {
            for direction in Direction::all() {
                let mut next_state = self.clone();
                next_state.cards.remove(i);

                if !match card {
                    Card::Move(strength) => next_state.try_move(global, direction, *strength),
                    Card::Jump(strength) => next_state.try_jump(global, direction, *strength),
                } {
                    // Invalid state, try next direction
                    continue;
                }

                next_states.push((1, Step { card: *card, direction }, next_state));
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

    fn is_valid(&self, _global: &Global) -> bool {
        true
    }

    fn is_solved(&self, global: &Global) -> bool {
        self.ball == global.flag
    }

    fn stringify(&self, global: &Global) -> String {
        let mut output = String::new();

        output.push_str("Map:\n");

        for y in 0..global.height {
            for x in 0..global.width {
                if self.ball.x == x as isize && self.ball.y == y as isize {
                    output.push('(');
                } else if global.flag.x == x as isize && global.flag.y == y as isize {
                    output.push('[');
                } else {
                    output.push(' ');
                }


                let tile = global.tiles[y * global.width + x];
                output.push(match tile {
                    Tile::Empty => ' ',
                    Tile::Flat(height) => std::char::from_digit(height as u32, 10).unwrap(),
                    Tile::Slope(_, Direction::Up) => '^',
                    Tile::Slope(_, Direction::Right) => '>',
                    Tile::Slope(_, Direction::Down) => 'v',
                    Tile::Slope(_, Direction::Left) => '<',
                    Tile::Water(_) => todo!(),
                    Tile::Sand(_) => todo!(),
                });

                if self.ball.x == x as isize && self.ball.y == y as isize {
                    output.push(')');
                } else if global.flag.x == x as isize && global.flag.y == y as isize {
                    output.push(']');
                } else {
                    output.push(' ');
                }
            }

            output.push('\n');
        }

        output.push_str(&format!("Cards: {:?}", self.cards));

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

fn stringify_solution(
    solver: &Solver<Global, Local, Step>,
    initial_state: &Local,
    solved_state: &Local,
) -> String {
    let mut path = String::new();
    
    let mut path = String::new();
    for step in solver.path(initial_state, solved_state).unwrap().iter() {
        path.push_str(&format!("{step:?}\n"));
    }

    path
}

fn main() -> Result<()> {
    env_logger::init();

    let stdin = std::io::stdin();
    let mut reader = stdin.lock();

    let (global, local) = Global::read(&mut reader).unwrap();
    dbg!(&global, &local);

    log::info!("Initial state:\n{}", local.stringify(&global));

    // If there is an arg, assume it's a test case and run it
    if std::env::args().len() > 1 {
        std::env::args().skip(1).for_each(|instructions| {
            todo!();
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
    use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
    use std::{fs::File, io::Read, sync::mpsc, thread};

    use super::*;

    #[test]
    fn test_all_solutions() {
        // Timeout after 1 second or GOLF_PEAKS_TEST_TIMEOUT if set
        let timeout = std::time::Duration::from_secs(
            std::env::var("GOLF_PEAKS_TEST_TIMEOUT")
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(1),
        );

        // Collect all tests to run in order
        let mut test_files = Vec::new();

        let folders = ["data/GOLF_PEAKS", "data/GOLF_PEAKS/primer"];

        for folder in folders.iter() {
            for entry in std::fs::read_dir(folder).unwrap() {
                let entry = entry.unwrap();
                let path = entry.path();
                if !path.is_dir() {
                    continue;
                }

                for entry in std::fs::read_dir(&path).unwrap() {
                    let entry = entry.unwrap();
                    let path = entry.path();

                    if path.extension().is_none() || path.extension().unwrap() != "txt" {
                        continue;
                    }

                    test_files.push(path);
                }
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
