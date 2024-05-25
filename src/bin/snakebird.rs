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

        Ok((
            Global {
                width,
                height,
                tiles,
            },
            local,
        ))
    }

    // Get the tile at a point
    fn tile(&self, point: Point) -> Tile {
        self.tiles.get(&point).copied().unwrap_or(Tile::Empty)
    }

    // Check if a point is in bounds
    fn in_bounds(&self, point: Point) -> bool {
        point.x >= 0
            && point.x < self.width as isize
            && point.y >= 0
            && point.y < self.height as isize
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

        // Cannot move out of bounds
        if !global.in_bounds(new_head) {
            return false;
        }

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
        // If any snake is pushed into a wall or out of bounds, the whole move is invalid
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

        // No snake pushing points can be out of bounds
        if snake_pushing_points
            .iter()
            .skip(1)
            .any(|p| !global.in_bounds(*p))
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


    fn is_supported_by_wall(&self, global: &Global, index: usize) -> bool {
        self.snakes[index].points.iter().any(|point| {
            global.tile(*point + Direction::Down) == Tile::Wall
        })
    }

    fn is_supported_by(&self, _: &Global, index_1: usize, index_2: usize) -> bool {
        self.snakes[index_1].points.iter().any(|point| {
            self.snakes[index_2].points.contains(&(*point + Direction::Down))
        })
    }

    fn try_gravity(&mut self, global: &Global) -> bool {
        // Apply gravity to all snakes
        'still_falling: loop {
            // Check if any snake can fall
            for (falling_snake_index, _) in self.snakes.iter().enumerate() {
                let can_fall = &self.snakes[falling_snake_index].points.iter().all(|point| {
                    let below = *point + Direction::Down;

                    // Apparently we can walk across fruit?
                    if self.fruit.contains(&below) {
                        return false;
                    }

                    // Otherwise, don't fall through walls
                    if self.is_supported_by_wall(global, falling_snake_index) {
                        return false;
                    }

                    // And don't fall on *other* snakes (we can't support ourselves)
                    for (other_snake_index, _) in self.snakes.iter().enumerate() {
                        if falling_snake_index == other_snake_index {
                            continue;
                        }

                        if self.is_supported_by(global, falling_snake_index, other_snake_index) {
                            // Deal with mutually supported snakes
                            // TODO: This is hacky..
                            if self.is_supported_by(global, other_snake_index, falling_snake_index) && !self.is_supported_by_wall(global, other_snake_index) {
                                return true;
                            }

                            // Otherwise, supported by a supported snake, cannot fall
                            return false;
                        }
                    }

                    true
                });
                if !can_fall {
                    continue;
                }

                // If we can fall and we're on the bottom, bad things happen
                // By that, I mean you lose
                if self.snakes[falling_snake_index]
                    .points
                    .iter()
                    .any(|point| !global.in_bounds(*point))
                {
                    return false;
                }

                // If we made it here, the snake is falling move down all points
                for point in self.snakes[falling_snake_index].points.iter_mut() {
                    point.y += 1;
                }

                // Also, spikes
                if self.snakes[falling_snake_index]
                    .points
                    .iter()
                    .any(|point| global.tile(*point) == Tile::Spike)
                {
                    return false;
                }

                // If the snake's head is somehow on the exit, it ... exits? 
                if global.tile(*self.snakes[falling_snake_index].points.first().unwrap()) == Tile::Exit {
                    self.snakes.remove(falling_snake_index);
                }

                continue 'still_falling;
            }

            // If we made it here, no snake can fall anymore
            break;
        }

        // No snakes fell out of the world
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

        println!("{}", local.stringify(&global));

        // Hang off platform
        assert!(local.try_move(&global, 0, Direction::Right));

        println!("{}", local.stringify(&global));

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
        // All snakes must be in bounds
        self.snakes
            .iter()
            .all(|snake| snake.points.iter().all(|point| global.in_bounds(*point)))
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

fn main() {
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
        return;
    }

    // Otherwise, read from stdin and solve

    let mut solver = Solver::new(global.clone(), local.clone());
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
        return;
    }
    let solution = solution.unwrap();

    log::info!(
        "Solved after {} states in {} seconds:\n{}",
        solver.states_checked(),
        solver.time_spent(),
        solution.stringify(&global),
    );

    let mut last_moved_snake = '\0';

    let mut path = String::new();
    for step in solver.path(&local, &solution).unwrap().iter() {
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
    println!("{path}");
}

// Tests (todo)
/*
cargo --bin snakebird
for i in $(seq -f "%02g" 1 12)
do
  echo -n "$i\t"
  cat data/snakebird/primer/$i.txt | ./target/debug/snakebird
done

01	0→→→→→→↑→→→→↑
02	0→→→↑→→↑→→→↑←←↑←←↑↑←←←←↑↑↑←←
03	0→→→↓↓←←←←←←↑↑→↑→→→→↓→→↓→→
04	0↑→→→→↑←↑←←←←↑↑↑→↑
05	0→→→↓→↑←←↑←
06	0→→↑→↑→↑↑↑←←↓↓←↑↑↑←←←↓↓←↓←←↓↓→→→↓→↑↑
07	0↑←←↑←←↑↑→↑↑→→↓↓→↑↑←↑↑←←↑↑→↑←←←↓→→→↑↑↑→↑↑↑←↑
08	0→→→→→→→↑←←←←←↑←←↑↑→→↑→→→→→↓→→→↑→↑↑↑←←↑↑←↑↑←←←↓↓←←↑←↑←←
09	0←←↑←←↑↑←←←←←↑↑→↓↓←←←←←↑↑←←↑→↑→↑
10	0→→→↑→→→→↑←←←←↑←←↑↑→↑→→
11	0↑a←←←←←↑←0↑←a←←0←↑
12	0→→→→a↑0→a↑0→↑
13	a→→↑↑→→→↓↓←↑←0→a←0↑a←←0↑←←
14	0→→→↑←←←←a↑0←←a←0←←
15	0↑a←←↑←←←←0←←a↑0←a←0←↑←
16  0→→↑→a↑↑↑0↑↑↑←←↑→→→→→↑→→a→→→↑→→→0↑
17  0→→→→→→↑↑←↑→→→→→→↑←←↑←←←←↑→↑↑←←←←←←↑←
18  0↑←↑←←←←↑↑↑→→→→→↓←↓↓←←←↓→↓↓→→↑→→↓↓→→→→↑↑↑
19  

x1  0→→→↓↓←←↑←←↑←←↓↓→→↓↓↓→→↑↑↑→→↑↑←←↑←↑
x2  A↑{→A←a→→A↑a↑{↑a→0↑A→→↑←←a→A↓0→a→{↓→a→→0↑
*/
