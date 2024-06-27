use std::fmt;
use std::io;

use solver::{Solver, State, Point};

const MOVE_EMPTY_COST: i64 = 1;
const MOVE_SNOW_COST: i64 = 1;
const TELEPORT_COST: i64 = 10;
const PUSH_EMPTY_COST: i64 = 2;
const PUSH_SNOW_COST: i64 = 4;
const STACK_COST: i64 = 10;
const UNSTACK_EMPTY_COST: i64 = 4;
const UNSTACK_SNOW_COST: i64 = 4;

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
    Teleporter,
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct Map {
    width: u8,
    height: u8,
    player: Point,
    data: Vec<Location>,
    targets: Option<Vec<Point>>,
    teleporters: Vec<Point>,
}

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
enum Step {
    North,
    South,
    East,
    West,
}

impl Map {
    fn get(&self, p: Point) -> Option<Location> {
        self.get_xy(p.x as i8, p.y as i8)
    }

    fn set(&mut self, p: Point, l: Location) {
        self.data[(p.x * (self.width as isize) + p.y) as usize] = l;
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
        let mut player = Point{x: 0, y: 0};

        let mut targets = Vec::new();
        let mut teleporters = Vec::new();

        for (y, line) in lines.iter().enumerate() {
            for (x, c) in line.chars().enumerate() {
                match c {
                    '-' => data.push(Empty),
                    '=' => {
                        data.push(Empty);
                        targets.push(Point{x: x as isize, y: y as isize});
                    }
                    '*' => data.push(Snow),
                    '+' => {
                        data.push(Snow);
                        targets.push(Point{x: x as isize, y: y as isize});
                    }
                    'X' => data.push(Wall),
                    '#' => {
                        data.push(Empty);
                        player = Point{x: x as isize, y: y as isize};
                    }
                    'T' => {
                        data.push(Teleporter);
                        teleporters.push(Point{x: x as isize, y: y as isize});
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
            targets: if targets.is_empty() {
                None
            } else {
                Some(targets)
            },
            teleporters,
        }
    }
}

impl fmt::Display for Map {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        use Location::*;

        let mut s: String = String::new();

        s += "map:\n";
        for y in 0..(self.height as isize) {
            s += "  ";
            for x in 0..(self.width as isize) {
                if x == self.player.x && y == self.player.y {
                    s += "#";
                } else {
                    match self.get_xy(x as i8, y as i8).unwrap() {
                        Empty => s += "-",
                        Snow => s += "*",
                        Wall => s += "X",
                        Teleporter => s += "T",
                        Snowman(stack) => s += (stack as u8).to_string().as_str(),
                    }
                }
            }
            s += "\n";
        }

        write!(f, "{}", s)
    }
}

impl<G> State<G, Step> for Map {
    fn is_valid(&self, g: &G) -> bool {
        use Location::*;
        use Stack::*;

        // If we're solved, we're valid
        if self.is_solved(g) {
            return true;
        }

        // Check that we have or can build the proper number of snowmen
        let mut small_balls = 0;
        let mut medium_balls = 0;
        let mut large_balls = 0;
        let mut snow_count = 0;

        for y in 0..self.height {
            for x in 0..self.width {
                match self.get(Point{x: x as isize, y: y as isize}) {
                    Some(Snowman(any)) => {
                        small_balls += if (any as u8) & 1 != 0 { 1 } else { 0 };
                        medium_balls += if (any as u8) & 2 != 0 { 1 } else { 0 };
                        large_balls += if (any as u8) & 4 != 0 { 1 } else { 0 };
                    }
                    Some(Snow) => snow_count += 1,
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

        // We have to have enough snow left to transform so that we have N L/M/S
        // Assume we have the correct number of balls.
        // lol
        let target_snowmen = (small_balls + medium_balls + large_balls) / 3;

        // Not enough remaining heads
        if small_balls < target_snowmen {
            return false;
        }

        // Not enough heads to match bellies (even if we used all remaining snow)
        if medium_balls + small_balls.min(snow_count) < target_snowmen {
            return false;
        }

        // Not enough butts, even if we use bellies and heads with remaining snow
        // This is imperfect, but should help
        if large_balls + medium_balls.min(snow_count) + small_balls.min(snow_count) < target_snowmen
        {
            return false;
        }

        // We have too many bottoms; can't downgrade them
        if large_balls > target_snowmen {
            return false;
        }

        // This overestimates
        // Technically a snowman is 'moveable' if you can either stack them
        // Or the target snowman is itself moveable
        // Which we can calculate, but I haven't yet
        let is_moveable = |l: Location| match l {
            Empty => true,
            Snow => true,
            Snowman(_) => true,
            _ => false,
        };

        // Any non-large in a corner is trapped (and thus invalid)
        // Larges in a non-target corner are also bad
        for y in 0..self.height {
            for x in 0..self.width {
                let (is_snowman, is_large) = match self.get(Point{x: x as isize, y: y as isize}) {
                    Some(Snowman(LargeMediumSmall | LargeMedium | Large)) => (true, true),
                    Some(Snowman(_)) => (true, false),
                    _ => (false, false),
                };
                if !is_snowman {
                    continue;
                }
                if is_large && self.targets.is_some() {
                    if self
                        .targets
                        .as_ref()
                        .unwrap()
                        .contains(&Point{x: x as isize, y: y as isize})
                    {
                        continue;
                    }
                }

                let north = self.get(Point{x: x as isize, y: y as isize - 1}).unwrap_or(Wall);
                let south = self.get(Point{x: x as isize, y: y as isize + 1}).unwrap_or(Wall);
                let east = self.get(Point{x: x as isize + 1, y: y as isize}).unwrap_or(Wall);
                let west = self.get(Point{x: x as isize - 1, y: y as isize}).unwrap_or(Wall);

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

        // If there are targets, any completed snowman not on a target is invalid
        if let Some(targets) = &self.targets {
            for y in 0..self.height {
                for x in 0..self.width {
                    if let Some(Snowman(LargeMediumSmall)) = self.get(Point{x: x as isize, y: y as isize}) {
                        if !targets.contains(&Point{x: x as isize, y: y as isize}) {
                            return false;
                        }
                    }
                }
            }
        }

        // Otherwise, assume valid
        return true;
    }

    fn is_solved(&self, _: &G) -> bool {
        use Location::*;
        use Stack::*;

        // All snowmen are completely formed
        for y in 0..self.height {
            for x in 0..self.width {
                match self.get(Point{x: x as isize, y: y as isize}) {
                    Some(Snowman(LargeMediumSmall)) => continue,
                    Some(Snowman(_)) => return false,
                    _ => continue,
                }
            }
        }

        // If there are targets, any completed snowman not on a target is invalid
        // TODO: This probably shouldn't be in both is_valid and is_solved
        if let Some(targets) = &self.targets {
            for y in 0..self.height {
                for x in 0..self.width {
                    if let Some(Snowman(LargeMediumSmall)) = self.get(Point{x: x as isize, y: y as isize}) {
                        if !targets.contains(&Point{x: x as isize, y: y as isize}) {
                            return false;
                        }
                    }
                }
            }
        }

        return true;
    }

    fn next_states(&self, _: &G) -> Option<Vec<(i64, Step, Map)>> {
        use Location::*;
        use Stack::*;
        use Step::*;

        let mut states = Vec::new();

        let moves = [
            (North, Point{x: 0, y: -1}),
            (South, Point{x: 0, y: 1}),
            (East, Point{x: 1, y: 0}),
            (West, Point{x: -1, y: 0}),
        ];

        for (step, delta) in moves.into_iter() {
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

            // The target is a teleporter, queue all others
            if self.teleporters.contains(&next) {
                for teleporter in self.teleporters.iter() {
                    if *teleporter == next {
                        continue;
                    }

                    let mut new_state = self.clone();
                    new_state.player = *teleporter;
                    states.push((TELEPORT_COST, step, new_state));
                }

                continue;
            }

            // The target is empty, just move into it
            if let Empty = target.unwrap() {
                let mut new_state = self.clone();
                new_state.player = next;
                states.push((MOVE_EMPTY_COST, step, new_state));
                continue;
            }

            // Likewise onto snow
            if let Snow = target.unwrap() {
                let mut new_state = self.clone();
                new_state.player = next;
                states.push((MOVE_SNOW_COST, step, new_state));
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
                        states.push((PUSH_EMPTY_COST, step, new_state));
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
                        states.push((PUSH_SNOW_COST, step, new_state));
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
                            states.push((STACK_COST, step, new_state));
                            continue;
                        }
                    }
                }

                // Stacks can be broken apart if the target is empty
                // Note: In these cases, the player *does not move*
                if target_2.is_some()
                    && (any == MediumSmall || any == LargeSmall || any == LargeMedium)
                {
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
                        states.push((UNSTACK_EMPTY_COST, step, new_state));
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
                        states.push((UNSTACK_SNOW_COST, step, new_state));
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
        //         match self.get(Point{x: x as isize, y: y as isize}) {
        //             Some(Snowman(LargeMediumSmall)) => continue,
        //             Some(Snowman(_)) => {
        //                 for y2 in 0..self.height {
        //                     for x2 in 0..self.width {
        //                         match self.get(Point{x: x2 as i8, y: y2 as i8}) {
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
                match self.get(Point{x: x as isize, y: y as isize}) {
                    Some(Snowman(LargeMediumSmall)) => continue,
                    Some(Snowman(_)) => {
                        let dx = (self.player.x - (x as isize)).abs() as i64;
                        let dy = (self.player.y - (y as isize)).abs() as i64;
                        distance = distance.min(dx + dy);
                    }
                    _ => continue,
                }
            }
        }
        score += distance;

        // Add the distance from each small/medium to the nearest available medium/large
        for y in 0..self.height {
            for x in 0..self.width {
                match self.get(Point{x: x as isize, y: y as isize}) {
                    Some(Snowman(Small)) => {
                        let mut distance = (self.width as i8 + self.height as i8) as i64;
                        for y2 in 0..self.height {
                            for x2 in 0..self.width {
                                match self.get(Point{x: x2 as isize, y: y2 as isize}) {
                                    Some(Snowman(Medium)) | Some(Snowman(LargeMedium)) => {
                                        let dx = (x as i8 - x2 as i8).abs() as i64;
                                        let dy = (y as i8 - y2 as i8).abs() as i64;
                                        distance = distance.min(dx + dy);
                                    }
                                    _ => continue,
                                }
                            }
                        }
                        score += distance;
                    }
                    Some(Snowman(Medium)) => {
                        let mut distance = (self.width as i8 + self.height as i8) as i64;
                        for y2 in 0..self.height {
                            for x2 in 0..self.width {
                                match self.get(Point{x: x2 as isize, y: y2 as isize}) {
                                    Some(Snowman(Large)) => {
                                        let dx = (x as i8 - x2 as i8).abs() as i64;
                                        let dy = (y as i8 - y2 as i8).abs() as i64;
                                        distance = distance.min(dx + dy);
                                    }
                                    _ => continue,
                                }
                            }
                        }
                        score += distance;
                    }
                    _ => continue,
                }
            }
        }

        // // Add various points per incomplete snowman
        // for y in 0..self.height {
        //     for x in 0..self.width {
        //         match self.get(Point{x: x as isize, y: y as isize}) {
        //             Some(Snowman(LargeMediumSmall)) => continue,
        //             Some(Snowman(any)) => score += 10 * (any as i64),
        //             _ => continue,
        //         }
        //     }
        // }

        // // Add points for snow
        // for y in 0..self.height {
        //     for x in 0..self.width {
        //         match self.get(Point{x: x as isize, y: y as isize}) {
        //             Some(Snow) => score += 1,
        //             _ => continue,
        //         }
        //     }
        // }

        return score;
    }

    fn stringify(&self, _global: &G) -> String {
        format!("{}", self)
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
        log::info!("{solver}, state:\n{}", state.stringify(&()));
    }

    let solution = solver.get_solution();
    if solution.is_none() {
        println!("no solution found");
        return;
    }
    let solution = solution.unwrap();

    println!("solution: {}", solution);

    let mut steps = String::new();
    for step in solver.path(&initial_state, &solution).unwrap() {
        match step {
            Step::North => steps += "W",
            Step::South => steps += "S",
            Step::East => steps += "D",
            Step::West => steps += "A",
        }
    }
    println!("path: {}", steps);

    println!(
        "{} states, {} seconds",
        solver.states_checked(),
        solver.time_spent()
    );
}

// TODO: Rebecca
