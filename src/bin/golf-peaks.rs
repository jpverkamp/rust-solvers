use anyhow::{anyhow, Result};
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
    #[allow(dead_code)]
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
        vec![
            Direction::Up,
            Direction::Down,
            Direction::Left,
            Direction::Right,
        ]
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
enum AngleType {
    TopLeft,
    TopRight,
    BottomLeft,
    BottomRight,
}

impl AngleType {
    fn try_reflect(&self, direction: Direction) -> Option<Direction> {
        match self {
            AngleType::TopLeft => match direction {
                Direction::Up => Some(Direction::Right),
                Direction::Left => Some(Direction::Down),
                _ => None,
            },
            AngleType::TopRight => match direction {
                Direction::Up => Some(Direction::Left),
                Direction::Right => Some(Direction::Down),
                _ => None,
            },
            AngleType::BottomLeft => match direction {
                Direction::Down => Some(Direction::Right),
                Direction::Left => Some(Direction::Up),
                _ => None,
            },
            AngleType::BottomRight => match direction {
                Direction::Down => Some(Direction::Left),
                Direction::Right => Some(Direction::Up),
                _ => None,
            },
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Tile {
    Empty,
    Flat(usize),
    Slope(usize, Direction),
    Angle(usize, AngleType),
    Sand(usize),
    Quicksand(usize),
    Water(usize),
    Spring(usize),
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
enum CardStep {
    Move(usize),
    Jump(usize),
    None,
}

const CARD_MAX_STEPS: usize = 3;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Card([CardStep; CARD_MAX_STEPS]);

impl From<&str> for Card {
    fn from(value: &str) -> Card {
        let mut card = Card([CardStep::None; CARD_MAX_STEPS]);
        let mut index = 0;

        // As soon as the string is empty, we're done
        let mut value = value;
        while !value.is_empty() {
            // The first part will be a number of strength
            let strength = value
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .collect::<String>();
            let rest = &value[strength.len()..];

            let strength = strength.parse().expect("Invalid card strength");

            // The next part will be a card type
            match rest.chars().next().expect("Invalid card type") {
                '-' => card.0[index] = CardStep::Move(strength),
                '/' => card.0[index] = CardStep::Jump(strength),
                _ => panic!("Invalid card type"),
            }

            // Advance rest + 1 and index
            value = &rest[1..];
            index += 1;
            if index >= CARD_MAX_STEPS {
                panic!("Too many steps in card, max is {CARD_MAX_STEPS}");
            }
        }

        card
    }
}

impl From<Card> for String {
    fn from(value: Card) -> Self {
        let mut output = String::new();

        for step in value.0.iter() {
            match step {
                CardStep::Move(strength) => output.push_str(&format!("{strength}-")),
                CardStep::Jump(strength) => output.push_str(&format!("{strength}/")),
                CardStep::None => break,
            }
        }

        output
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct Local {
    ball: Point,
    cards: Vec<Card>,
    last_safe: Point,
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

        let width: usize = size[0];
        let height: usize = size[1];

        let mut tiles = vec![Tile::Empty; width * height];
        let mut flag = None;
        let mut ball = None;
        let mut solutions = Vec::new();

        // Read tiles (plus flag and ball)
        for y in 0..height {
            let line = lines
                .next()
                .expect("Missing tile line")
                .expect("Missing tile line");

            for (x, mut definition) in line.split_ascii_whitespace().enumerate() {
                // Empty space, already handled
                if definition == "x" {
                    continue;
                }

                let mut part;

                enum ParsedType { 
                    Flat,
                    Flag,
                    Ball,
                    Slope,
                    Angle,
                    Sand,
                    Quicksand,
                    Water,
                    Spring
                }

                let mut height = 0;
                let mut current_type = ParsedType::Flat;
                let mut which_angle = AngleType::TopLeft;
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
                    } else if let Ok(new_direction) = Direction::try_from(part) {
                        current_type = ParsedType::Slope;
                        slope_direction = new_direction;
                    } else {
                        match part {
                            "flag" => current_type = ParsedType::Flag,
                            "ball" => current_type = ParsedType::Ball,
                            "sand" => current_type = ParsedType::Sand,
                            "quick" => current_type = ParsedType::Quicksand,
                            "water" => current_type = ParsedType::Water,
                            "boing" => current_type = ParsedType::Spring,
                            "/tl" => {
                                current_type = ParsedType::Angle;
                                which_angle = AngleType::TopLeft;
                            },
                            "/tr" => {
                                current_type = ParsedType::Angle;
                                which_angle = AngleType::TopRight;
                            },
                            "/bl" => {
                                current_type = ParsedType::Angle;
                                which_angle = AngleType::BottomLeft;
                            },
                            "/br" => {
                                current_type = ParsedType::Angle;
                                which_angle = AngleType::BottomRight;
                            },
                            _ => panic!("Unknown tile definition `{part}` in `{line}`")
                        }
                    }
                }

                tiles[y * width + x] = match current_type {
                    ParsedType::Flat => Tile::Flat(height),
                    ParsedType::Flag => {
                        assert!(flag.is_none(), "Multiple flags in map");
                        flag = Some(Point {
                            x: x as isize,
                            y: y as isize,
                        });
                        Tile::Flat(height)
                    },
                    ParsedType::Ball => {
                                            assert!(ball.is_none(), "Multiple balls in map");
                        ball = Some(Point {
                            x: x as isize,
                            y: y as isize,
                        });
                        Tile::Flat(height)
                    },
                    ParsedType::Slope => Tile::Slope(height, slope_direction),
                    ParsedType::Angle => Tile::Angle(height, which_angle),
                    ParsedType::Sand => Tile::Sand(height),
                    ParsedType::Quicksand => Tile::Quicksand(height),
                    ParsedType::Water => Tile::Water(height),
                    ParsedType::Spring => Tile::Spring(height),
                };
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
            cards.push(Card::from(card));
        }

        // Read another empty line then any solutions
        lines.next();

        for line in lines {
            let line = line.expect("Missing solution line").trim().to_string();
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
                flag: flag.expect("No flag in map"),
                solutions,
            },
            Local {
                ball: ball.expect("No ball in map"),
                cards,
                last_safe: ball.expect("No ball in map"),
            },
        ))
    }
}

impl Global {
    fn tile_at(&self, point: Point) -> Tile {
        // Points out of bounds are always empty
        if point.x < 0
            || point.y < 0
            || point.x >= self.width as isize
            || point.y >= self.height as isize
        {
            return Tile::Empty;
        }

        self.tiles[point.y as usize * self.width + point.x as usize]
    }
}

impl Local {
    fn try_card(&mut self, global: &Global, card: Card, direction: Direction) -> bool {
        let mut direction = direction;
        log::debug!("try_card({:?}, {card:?}, {direction:?})", self.ball);

        for card_step in card.0.iter() {
            log::debug!("try_card({:?}, {card:?}, {direction:?}), step={card_step:?}", self.ball);

            let success = match card_step {
                CardStep::Move(strength) => self.try_move(global, direction, *strength),
                CardStep::Jump(strength) => self.try_jump(global, direction, *strength),
                CardStep::None => break,
            };
            if !success {
                return false;
            }

            // Halfway through a move, if we're on a slope change direction to match it
            // This comes up in 3-10; I don't think I like it
            if let Tile::Slope(_, slope_direction) = global.tile_at(self.ball) {
                direction = slope_direction;
            }

            // If after any step, we're on water, reset to last safe tile and stop move
            // This shouldn't apply after move (it's handled in try_move), but might after jump or slide
            if let Tile::Water(_) = global.tile_at(self.ball) {
                self.ball = self.last_safe;
                return true;
            }
        }

        if !self.try_slopes(global) {
            return false;
        }

        // On quicksand, we fail
        if let Tile::Quicksand(_) = global.tile_at(self.ball) {
            return false;
        }

        true
    }

    fn try_move(&mut self, global: &Global, direction: Direction, strength: usize) -> bool {
        let current_tile = global.tile_at(self.ball);

        // If we're on a flat/safe tile, mark this as the last safe spot
        if let Tile::Flat(_) | Tile::Angle(_, _) | Tile::Sand(_) = current_tile {
            self.last_safe = self.ball;
        }

        // No more moving to do, we're done
        if strength == 0 {
            log::debug!("try_move({:?}, {direction:?}, {strength}), base case", self.ball);
            return true;
        }

        let current_height = match current_tile {
            Tile::Empty => unreachable!(),
            
            Tile::Flat(height)
            | Tile::Angle(height, _)
            | Tile::Sand(height)
            | Tile::Quicksand(height)
            | Tile::Water(height) 
            | Tile::Spring(height) => height,
            
            Tile::Slope(height, slope_direction) => {
                if direction == slope_direction {
                    height
                } else if direction == slope_direction.flip() {
                    height + 1 // TODO: Is this correct?
                } else {
                    return false; // TODO: Support this?
                }
            }
        };

        log::debug!(
            "try_move({self:?}, {direction:?}, {strength}), height={current_height}, tile={:?}",
            global.tile_at(self.ball)
        );

        // Cannot slide out of sand, just stop moving
        // But return true, this isn't an error, just stopping
        if let Tile::Sand(_) = current_tile {
            return true;
        }

        // If we're currently on an angled tile, we need to reflect
        // When this recurs, we'll be moving 'out' of the tile so won't trigger it twice
        // TODO: Can you fall onto an angled tile?
        if let Tile::Angle(_, a_type) = current_tile {
            if let Some(new_direction) = a_type.try_reflect(direction) {
                return self.try_move(global, new_direction, strength);
            }
        }

        let next_point = self.ball + Point::from(direction);
        let mut next_tile = global.tile_at(next_point);

        // Trying to move into empty space/out of bounds
        if let Tile::Empty = next_tile {
            return false;
        }

        // Trying to move into/up/down a slope
        if let Tile::Slope(height, slope_direction) = next_tile {
            // Sliding down at the proper height
            if height == current_height && slope_direction == direction {
                self.ball = next_point;
                return self.try_move(global, direction, strength - 1);
            }

            // Sliding up at the proper height
            if height == current_height + 1 && slope_direction == direction.flip() {
                self.ball = next_point;
                return self.try_move(global, direction, strength - 1);
            }

            // Dropping down onto a slope
            if height < current_height {
                self.ball = next_point;
                return self.try_move(global, direction, strength - 1);
            }

            // Anything else, we'll fall through as a normal flat tile
            next_tile = Tile::Flat(height);
        }

        // Normal flat tile (or equivalent, like quicksand or springs)
        match next_tile {
            Tile::Flat(height) | Tile::Quicksand(height) => {
                // On the same level, just move
                if height == current_height {
                    self.ball = next_point;
                    return self.try_move(global, direction, strength - 1);
                }

                // New tile is higher, bounce
                // This effectively reverses direction and moves 'back' to the same tile
                if height > current_height {
                    // Bouncing off a spring is a weird special case
                    if let Tile::Spring(_) = current_tile {
                        return self.try_jump(global, direction.flip(), strength - 1);
                    }

                    return self.try_move(global, direction.flip(), strength - 1);
                }

                // New tile is lower, fall onto that level and continue
                if height < current_height {
                    self.ball = next_point;
                    return self.try_move(global, direction, strength - 1);
                }

                unreachable!();
            }
            _ => {}
        }

        // Trying to move onto water, fall back to last safe tile and end move
        if let Tile::Water(_) = next_tile {
            self.ball = self.last_safe;
            return true;
        }

        // If we move onto a spring, treat the rest of the move as a jump
        if let Tile::Spring(_) = next_tile {
            self.ball = self.ball + direction.into();
            return self.try_jump(global, direction, strength - 1);
        }

        // Normal flat tile, recur
        self.ball = self.ball + direction.into();
        self.try_move(global, direction, strength - 1)
    }

    fn try_jump(&mut self, global: &Global, direction: Direction, strength: usize) -> bool {
        log::debug!("try_jump({:?}, {direction:?}, {strength})", self.ball);

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

    fn try_slopes(&mut self, global: &Global) -> bool {
        let current_tile = global.tile_at(self.ball);
        if current_tile == Tile::Empty {
            return false;
        }

        log::debug!("try_slopes({:?}) @ {current_tile:?}", self.ball);

        // Slopes apply a single tile move than recur
        if let Tile::Slope(_, slope_direction) = current_tile {
            if !self.try_move(global, slope_direction, 1) {
                return false;
            }

            return self.try_slopes(global);
        }

        // Any non-slopes just don't slide
        true
    }
}

impl State<Global, Step> for Local {
    fn next_states(&self, global: &Global) -> Option<Vec<(i64, Step, Self)>>
    where
        Self: Sized,
    {
        let mut next_states = Vec::new();

        for (i, card) in self.cards.iter().enumerate() {
            'next_direction: for direction in Direction::all() {
                let mut next_state = self.clone();
                next_state.cards.remove(i);

                if !next_state.try_card(global, *card, direction) {
                    // Invalid state, try next direction
                    continue 'next_direction;
                }

                next_states.push((
                    1,
                    Step {
                        card: *card,
                        direction,
                    },
                    next_state,
                ));
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

                    Tile::Angle(_, AngleType::TopLeft) => '◤',
                    Tile::Angle(_, AngleType::TopRight) => '◥',
                    Tile::Angle(_, AngleType::BottomLeft) => '◣',
                    Tile::Angle(_, AngleType::BottomRight) => '◢',

                    Tile::Sand(_) => '▒',
                    Tile::Quicksand(_) => '▓',
                    Tile::Water(_) => '≈',
                    Tile::Spring(_) => '⌉',
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

        output.push_str(&format!(
            "Cards: {}",
            self.cards
                .iter()
                .map(|c| String::from(*c))
                .collect::<Vec<String>>()
                .join(" ")
        ));

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

    for step in solver.path(initial_state, solved_state).unwrap().iter() {
        let c = String::from(step.card);
        let d = match step.direction {
            Direction::Up => '↗',
            Direction::Down => '↙',
            Direction::Left => '↖',
            Direction::Right => '↘',
        };

        path.push_str(&format!("{c}{d} "));
    }

    path.trim().to_string()
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

            for step in instructions.split_ascii_whitespace() {
                log::info!("Step: {step} on {}", local.stringify(&global));

                // All of the unicode chars except the last encode a card
                let mut chars = step.chars();
                let d = match chars.next_back().expect("No direction") {
                    '↗' => Direction::Up,
                    '↙' => Direction::Down,
                    '↖' => Direction::Left,
                    '↘' => Direction::Right,
                    _ => panic!("Invalid direction"),
                };
                let card = Card::from(chars.as_str());

                assert!(
                    local.cards.contains(&card),
                    "Card doesn't exist in step `{step}`"
                );

                if !local.try_card(&global, card, d) {
                    panic!("Failed to move by {card:?} in `{step}`")
                }

                log::info!("After: {}", local.stringify(&global));
            }

            assert!(local.ball == global.flag, "Ball not at flag")
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

        let folders = ["data/golf-peaks"];

        for folder in folders.iter() {
            for entry in std::fs::read_dir(folder).unwrap() {
                let entry = entry.unwrap();
                let path = entry.path();

                if path.extension().is_none() || path.extension().unwrap() != "txt" {
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
