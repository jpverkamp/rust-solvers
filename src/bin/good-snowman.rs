use std::fmt;
use std::io;
use std::ops::Add;
use std::ops::Sub;

use solver::{Solver, State};

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
struct Point(i8, i8);

impl Add<Point> for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        Point(self.0 + other.0, self.1 + other.1)
    }
}

impl Sub<Point> for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        Point(self.0 - other.0, self.1 - other.1)
    }

}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
enum Stack {
    Small = 1,
    Medium = 2,
    MediumSmall = 3,
    Large = 4,
    LargeSmall = 5,
    LargeMedium = 6,
    LargeMediumSmall = 7,
}

impl From<u8> for Stack {
    fn from(value: u8) -> Self {
        match value {
            1 => Stack::Small,
            2 => Stack::Medium,
            3 => Stack::MediumSmall,
            4 => Stack::Large,
            5 => Stack::LargeSmall,
            6 => Stack::LargeMedium,
            7 => Stack::LargeMediumSmall,
            _ => panic!("Invalid stack size"),
        }
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
enum Location {
    Empty,
    Snow,
    Wall,
    Snowman(Stack),
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct Map {
    width: u8,
    height: u8,
    player: Point,
    data: Vec<Location>,
}

impl Map {
    fn get(&self, p: Point) -> Option<Location> {
        self.get_xy(p.0, p.1)
    }

    fn set(&mut self, p: Point, l: Location) {
        self.data[(p.1 * (self.width as i8) + p.0) as usize] = l;
    }

    fn get_xy(&self, x: i8, y: i8) -> Option<Location> {
        if x < 0 || x >= self.width as i8 || y < 0 || y >= self.height as i8 {
            return None;
        }

        Some(self.data[(y * self.width as i8 + x) as usize])
    }
}

impl From<&str> for Map {
    fn from(value: &str) -> Self {
        use Location::*;

        let lines = value.lines().collect::<Vec<&str>>();
        let height = lines.len() as u8;
        let width = lines[0].len() as u8;

        let mut data = Vec::new();
        let mut player = Point(0, 0);

        for (y, line) in lines.iter().enumerate() {
            for (x, c) in line.chars().enumerate() {
                match c {
                    '-' => data.push(Empty),
                    '*' => data.push(Snow),
                    'X' => data.push(Wall),
                    '#' => {
                        data.push(Empty);
                        player = Point(x as i8, y as i8);
                    }
                    _ => {
                        let value = c.to_digit(10).unwrap() as u8;
                        data.push(Snowman(value.into()));
                    }
                }
            }
        }

        Map {
            width,
            height,
            player,
            data,
        }
    }
}

impl fmt::Display for Map {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        use Location::*;

        let mut s: String = String::new();

        s += "map:\n";
        for y in 0..(self.height as i8) {
            s += "  ";
            for x in 0..(self.width as i8) {
                if x == self.player.0 && y == self.player.1 {
                    s += "#";
                } else {
                    match self.get_xy(x, y).unwrap() {
                        Empty => s += "-",
                        Snow => s += "*",
                        Wall => s += "X",
                        Snowman(stack) => s += (stack as u8).to_string().as_str(),
                    }
                }
            }
            s += "\n";
        }

        write!(f, "{}", s)
    }
}

impl<G> State<G> for Map {
    fn is_valid(&self, g: &G) -> bool {
        use Location::*;
        use Stack::*;

        // If we're solved, we're valid
        if self.is_solved(g) {
            return true;
        }

        // We must have
        // Small >= Large
        // Small + Medium >= Large (since small can become medium)
        let mut small_balls = 0;
        let mut medium_balls = 0;
        let mut large_balls = 0;

        for y in 0..self.height {
            for x in 0..self.width {
                match self.get(Point(x as i8, y as i8)) {
                    Some(Snowman(any)) => {
                        small_balls += if (any as u8) & 1 != 0 { 1 } else { 0 };
                        medium_balls += if (any as u8) & 2 != 0 { 1 } else { 0 };
                        large_balls += if (any as u8) & 4 != 0 { 1 } else { 0 };
                    }
                    _ => continue,
                }
            }
        }

        if small_balls < large_balls {
            return false;
        }
        if small_balls + medium_balls < large_balls {
            return false;
        }

        // This overestimates
        // Technically a snowman is 'moveable' if you can either stack them 
        // Or the target snowman is itself moveable
        // Which we can calculate, but I haven't yet
        let is_moveable = |l: Location| {
            match l {
                Empty => true,
                Snow => true,
                Snowman(_) => true,
                _ => false,
            }
        };

        // Any non-large in a corner is trapped (and thus invalid)
        for y in 0..self.height {
            for x in 0..self.width {
                match self.get(Point(x as i8, y as i8)) {
                    Some(Snowman(LargeMediumSmall)) => continue,
                    Some(Snowman(LargeMedium)) => continue,
                    Some(Snowman(Large)) => continue,
                    Some(Snowman(_)) => {},
                    _ => continue,
                }

                let north = self.get(Point(x as i8, y as i8 - 1)).unwrap_or(Wall);
                let south = self.get(Point(x as i8, y as i8 + 1)).unwrap_or(Wall);
                let east = self.get(Point(x as i8 + 1, y as i8)).unwrap_or(Wall);
                let west = self.get(Point(x as i8 - 1, y as i8)).unwrap_or(Wall);

                if !is_moveable(north) && !is_moveable(west) {
                    return false;
                }
                if !is_moveable(north) && !is_moveable(east) {
                    return false;
                }
                if !is_moveable(south) && !is_moveable(west) {
                    return false;
                }
                if !is_moveable(south) && !is_moveable(east) {
                    return false;
                }
            }
        }

        // Otherwise, assume valid 
        return true;
    }

    fn is_solved(&self, _: &G) -> bool {
        use Location::*;
        use Stack::*;

        for y in 0..self.height {
            for x in 0..self.width {
                match self.get(Point(x as i8, y as i8)) {
                    Some(Snowman(LargeMediumSmall)) => continue,
                    Some(Snowman(_)) => return false,
                    _ => continue,
                }
            }
        }

        return true;
    }

    fn next_states(&self, _: &G) -> Option<Vec<(i64, Map)>> {
        use Location::*;
        use Stack::*;

        let mut states = Vec::new();

        let moves = [
            Point(0, -1),
            Point(0, 1),
            Point(1, 0),
            Point(-1, 0),
        ];

        for delta in moves.into_iter() {
            let next = self.player + delta;
            let next_2 = next + delta;

            let target = self.get(next);
            let target_2 = self.get(next_2);

            // Can't move out of bounds
            if target.is_none() {
                continue;
            }

            // Can't move into a wall
            if let Wall = target.unwrap() {
                continue;
            }

            // The target is empty, just move into it
            if let Empty = target.unwrap() {
                let mut new_state = self.clone();
                new_state.player = next;
                states.push((1, new_state));
                continue;
            }

            // Likewise onto snow
            if let Snow = target.unwrap() {
                let mut new_state = self.clone();
                new_state.player = next;
                states.push((1, new_state));
                continue;
            }

            // Try to push a snowman together or apart
            if let Snowman(any) = target.unwrap() {
                // Single balls can be pushed into empty spaces, snow (growing), and valid snowmen
                if target_2.is_some() && (any == Small || any == Medium || any == Large) {
                    // Empty spaces always work
                    if let Empty = target_2.unwrap() {
                        let mut new_state = self.clone();
                        new_state.player = next;
                        new_state.set(next, Empty);
                        new_state.set(next_2, Snowman(any));
                        states.push((1, new_state));
                        continue;
                    }

                    // Snow works and grows the ball
                    // Large stays large
                    // Snow is always removed (even if large -> large)
                    if let Snow = target_2.unwrap() {
                        let new_size = match any {
                            Small => Medium,
                            Medium => Large,
                            Large => Large,
                            _ => continue,
                        };

                        let mut new_state = self.clone();
                        new_state.player = next;
                        new_state.set(next, Empty);
                        new_state.set(next_2, Snowman(new_size));
                        states.push((1, new_state));
                        continue;
                    }

                    // Single balls can be pushed onto other smaller balls
                    if let Snowman(other) = target_2.unwrap() {
                        let combined = match (any, other) {
                            (Small, Medium) => Some(MediumSmall),
                            (Small, Large) => Some(LargeSmall),
                            (Medium, Large) => Some(LargeMedium),
                            (Small, LargeMedium) => Some(LargeMediumSmall),
                            _ => None,
                        };

                        if let Some(combined) = combined {
                            let mut new_state = self.clone();
                            new_state.player = next;
                            new_state.set(next, Empty);
                            new_state.set(next_2, Snowman(combined));
                            states.push((1, new_state));
                            continue;
                        }
                    }
                }
            
                // Stacks can be broken apart if the target is empty
                // Note: In these cases, the player *does not move*
                if target_2.is_some() && (any == MediumSmall || any == LargeSmall || any == LargeMedium) {
                    // Pushing into empty space
                    if let Empty = target_2.unwrap() {
                        let (a, b) = match any {
                            MediumSmall => (Medium, Small),
                            LargeSmall => (Large, Small),
                            LargeMedium => (Large, Medium),
                            _ => panic!("Invalid stack size"),
                        };

                        let mut new_state = self.clone();
                        new_state.set(next, Snowman(a));
                        new_state.set(next_2, Snowman(b));
                        states.push((1, new_state));
                        continue;
                    }

                    // Pushing onto snow
                    if let Snow = target_2.unwrap() {
                        let (a, b) = match any {
                            MediumSmall => (Medium, Medium),
                            LargeSmall => (Large, Medium),
                            LargeMedium => (Large, Large),
                            _ => panic!("Invalid stack size"),
                        };

                        let mut new_state = self.clone();        
                        new_state.set(next, Snowman(a));
                        new_state.set(next_2, Snowman(b));
                        states.push((1, new_state));
                        continue;
                    }
                }
            }
        
            // If we made it this far, the state is invalid for some other reason
            // println!("Invalid state: {:?}", self);
            // println!("Moving: {:?}", delta);
        }

        // If we have states, return them
        if !states.is_empty() {
            Some(states)
        } else {
            None
        }
    }

    fn heuristic(&self, _global: &G) -> i64 {
        use Location::*;
        use Stack::*;

        let mut score = 0;

        // // Add the sum of distances between each pair of incomplete snowmen
        // for y in 0..self.height {
        //     for x in 0..self.width {
        //         match self.get(Point(x as i8, y as i8)) {
        //             Some(Snowman(LargeMediumSmall)) => continue,
        //             Some(Snowman(_)) => {
        //                 for y2 in 0..self.height {
        //                     for x2 in 0..self.width {
        //                         match self.get(Point(x2 as i8, y2 as i8)) {
        //                             Some(Snowman(LargeMediumSmall)) => continue,
        //                             Some(Snowman(_)) => {
        //                                 let dx = x as i64 - x2 as i64;
        //                                 let dy = y as i64 - y2 as i64;
        //                                 score += dx.abs() + dy.abs();
        //                             },
        //                             _ => continue,
        //                         }
        //                     }
        //                 }
        //             },
        //             _ => continue,
        //         }
        //     }
        // }

        // Add the distance from the player to the nearest incomplete snowman
        let mut distance = (self.width as i8 + self.height as i8) as i64;
        for y in 0..self.height {
            for x in 0..self.width {
                match self.get(Point(x as i8, y as i8)) {
                    Some(Snowman(LargeMediumSmall)) => continue,
                    Some(Snowman(_)) => {
                        let dx = (self.player.0 - (x as i8)).abs() as i64;
                        let dy = (self.player.1 - (y as i8)).abs() as i64;
                        distance = distance.min(dx + dy);
                    },
                    _ => continue,
                }
            }
        }
        score += distance;
        
        // // Add various points per incomplete snowman
        // for y in 0..self.height {
        //     for x in 0..self.width {
        //         match self.get(Point(x as i8, y as i8)) {
        //             Some(Snowman(LargeMediumSmall)) => continue,
        //             Some(Snowman(any)) => score += 10 * (any as i64),
        //             _ => continue,
        //         }
        //     }
        // }

        return score;
    }

    fn display(&self, _global: &G) {
        println!("{}", self);
    }
}

fn main() {
    env_logger::init();

    let input = io::read_to_string(io::stdin()).unwrap();
    let initial_state = Map::from(input.as_str());

    println!("initial: {}", initial_state);

    let mut solver = Solver::new((), initial_state.clone());

    while let Some(state) = solver.next() {
        if solver.states_checked() % 10000 != 0 {
            continue;
        }

        println!("===== ===== ===== ===== =====");
        println!("state: {}", state);
        println!(
            "{} states checked, {} in queue, {} seconds, heuristic: {}",
            solver.states_checked(),
            solver.in_queue(),
            solver.time_spent(),
            state.heuristic(&()),
        );

        // if solver.states_checked() > 100 {
        //     break;
        // }
    }

    let solution = solver.get_solution();
    if solution.is_none() {
        println!("no solution found");
        return;
    }
    let solution = solution.unwrap();

    println!("solution: {}", solution);

    let mut steps = String::new();
    let mut player = initial_state.player;
    for step in solver.path(initial_state, solution).unwrap() {
        let delta = step.player - player;
        match delta {
            Point(0, -1) => steps += "N",
            Point(0, 1) => steps += "S",
            Point(1, 0) => steps += "E",
            Point(-1, 0) => steps += "W",
            Point(0, 0) => steps += "X",
            _ => panic!("Invalid move delta: {:?} to {:?} is {:?}", player, step.player, delta),
        }
        player = step.player;
    }
    println!("path: {}", steps);

    println!(
        "{} states, {} seconds",
        solver.states_checked(),
        solver.time_spent()
    );

}


// TODO: Rebecca