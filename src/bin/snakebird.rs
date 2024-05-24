use anyhow::{anyhow, Result};
use fxhash::FxHashMap;
use std::{io::{BufRead, Read}, ops::Add};

use solver::{Solver, State};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Tile {
    Empty,
    Wall,
    Exit,
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
            vec![], // a..=z
            vec![], // A..=Z
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
        point.x >= 0 && point.x < self.width as isize && point.y >= 0 && point.y < self.height as isize
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

impl Add<&Direction> for &Point {
    type Output = Point;

    fn add(self, direction: &Direction) -> Point {
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

impl State<Global, Step> for Local {
    fn is_valid(&self, global: &Global) -> bool {
        // All snakes must be in bounds
        self.snakes.iter().all(|snake| {
            snake.points.iter().all(|point| { global.in_bounds(*point) })
        })
    }

    fn is_solved(&self, _map: &Global) -> bool {
        // All snakes exited and fruit eaten
        self.snakes.is_empty() && self.fruit.is_empty()
    }

    fn next_states(&self, global: &Global) -> Option<Vec<(i64, Step, Local)>> {
        let mut next_states = Vec::new();

        // Try to move each snake in each direction 
        for (snake_index, _) in self.snakes.iter().enumerate() {
            'one_direction: for direction in [
                Direction::Up,
                Direction::Down,
                Direction::Left,
                Direction::Right,
            ]
            .iter()
            {
                // Figure out where the moving snake's head will move to
                let step = Step {
                    snake: snake_index,
                    direction: *direction,
                };
                let mut new_local = self.clone();
                let head = new_local.snakes[snake_index].points.first().unwrap();
                let new_head = head + direction;

                // Cannot move out of bounds
                if !global.in_bounds(new_head) {
                    continue 'one_direction;
                }

                // Cannot move into a wall
                if global.tile(new_head) == Tile::Wall {
                    continue 'one_direction;
                }

                // Cannot move into yourself
                if new_local.snakes[snake_index].points.contains(&new_head) {
                    continue 'one_direction;
                }

                // Attempt to push any snakes in the way
                let mut snake_pushing_indexes = Vec::new();
                let mut snake_pushing_points = vec![new_head];

                'daisies_remain: loop {
                    for (other_snake_index, _) in new_local.snakes.iter().enumerate() {
                        // TODO: What if we have a weird loop where we're pushing ourselves? 
                        if snake_index == other_snake_index {
                            continue;
                        }

                        // Already pushing this snake
                        if snake_pushing_indexes.contains(&other_snake_index) {
                            continue;
                        }

                        // If any pushing point is in the new snake, it's getting pushed too
                        if snake_pushing_points.iter().any(|p| new_local.snakes[other_snake_index].points.contains(&p)) {
                            snake_pushing_indexes.push(other_snake_index);
                            
                            snake_pushing_points.extend(
                                new_local.snakes[other_snake_index].points.iter().map(|p| p + direction)
                            );

                            continue 'daisies_remain;
                        }
                    }

                    break;
                }

                // No snake pushing points can hit anything (other than the original head)
                // TODO: Can you push a snake into the exit? 
                if snake_pushing_points.iter().skip(1).any(|p| global.tile(*p) != Tile::Empty) {
                    continue 'one_direction;
                }

                // No snake pushing points can be out of bounds
                if snake_pushing_points.iter().skip(1).any(|p| !global.in_bounds(*p)) {
                    continue 'one_direction;
                }

                // If we have snakes to push, move them all
                for snake_index in snake_pushing_indexes.iter() {
                    for point in new_local.snakes[*snake_index].points.iter_mut() {
                        *point = &*point + direction;
                    }
                }

                // If the currently moving snake moves onto the exit, check if it can leave
                // Otherwise, we're going to move as normal; potentially eating fruit
                if global.tile(new_head) == Tile::Exit {
                    // Check for exit

                    // Cannot exit (or even move onto exit) until all fruit is eaten
                    if !new_local.fruit.is_empty() {
                        continue 'one_direction;
                    }

                    // Otherwise, snake literally exits
                    // Bye bye snake
                    new_local.snakes.remove(snake_index);
                } else {
                    // Update the snake
                    new_local.snakes[snake_index].points.insert(0, new_head);

                    // Potentially eat fruit; don't remove tail if fruit was eaten
                    if let Some(fruit_index) =
                        new_local.fruit.iter().position(|fruit| fruit == &new_head)
                    {
                        new_local.fruit.remove(fruit_index);
                    } else {
                        new_local.snakes[snake_index].points.pop();
                    }
                }

                // Apply gravity to all snakes
                'still_falling: loop {
                    // Check if any snake can fall
                    for (snake_index, _) in new_local.snakes.iter().enumerate() {
                        let can_fall = &new_local.snakes[snake_index].points.iter().all(|point| {
                            let below = Point {
                                x: point.x,
                                y: point.y + 1,
                            };

                            // Apparently we can walk across fruit?
                            if new_local.fruit.contains(&below) {
                                return false;
                            }

                            // Otherwise, don't fall through walls
                            if global.tile(below) == Tile::Wall {
                                return false;
                            }

                            // And don't fall on *other* snakes (we can't support ourselves)
                            for (other_snake_index, other_snake) in self.snakes.iter().enumerate() {
                                if snake_index != other_snake_index
                                    && other_snake.points.contains(&below)
                                {
                                    return false;
                                }
                            }

                            true
                        });
                        if !can_fall {
                            break;
                        }

                        // If we can fall and we're on the bottom, bad things happen
                        // By that, I mean you lose
                        if new_local.snakes[snake_index]
                            .points
                            .iter()
                            .any(|point| !global.in_bounds(*point))
                        {
                            continue 'one_direction;
                        }

                        // If we made it here, the snake is falling move down all points
                        for point in new_local.snakes[snake_index].points.iter_mut() {
                            point.y += 1;
                        }

                        continue 'still_falling;
                    }

                    // If we made it here, no snake can fall anymore
                    break;
                }

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
                };

                for (i, snake) in self.snakes.iter().enumerate() {
                    for (j, pt) in snake.points.iter().enumerate() {
                        if pt.x == x && pt.y == y {
                            c = match i {
                                0 => ('0' as u8 + j as u8) as char,
                                1 => ('a' as u8 + j as u8) as char,
                                2 => ('A' as u8 + j as u8) as char,
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
// 01 0→→→→→→↑→→→→↑
// 02 0→→→↑→→↑→→→↑←←↑←←↑↑←←←←↑↑↑←←
// 03 0→→→↓↓←←←←←←↑↑→↑→→→→↓→→↓→→
// 04 0↑→→→→↑←↑←←←←↑↑↑→↑
// 05 0→→→↓→↑←←←↑
// 06 0→→↑→↑→↑←←←↑↑→→↑←←←←←↓←↓←↓←↓↓→→→↓→↑↑ // 5 seconds!
// 07 0↑←←↑←←↑↑→↑↑→→↓↓→↑↑←↑↑←←↑↑→↑←←←↓→→→↑↑↑→↑↑↑←↑
// 08 0→→→→→→→↑←←←←←↑←←↑↑→→↑→→→→→↓→→→↑→↑↑↑←←↑↑←↑↑←←←↓↓←←↑←↑←←
// 09 0←←↑←←↑↑←←←←←↑↑→↓↓←←←←←↑↑←←↑→↑→↑
// 10 0→→→↑→→→→↑←←←←↑←←↑↑→↑→→
// 11 a←0←a←←←←←0←a↑←↑0↑↑a←0←↑←
// 47 0→→→→↑→↑→↑←↓←←↓←↑←↑←
