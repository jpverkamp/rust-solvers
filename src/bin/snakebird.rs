use std::{collections::HashSet, io::{BufRead, Read}, ops::{Range, RangeInclusive}};
use anyhow::{Result, anyhow};

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
    tiles: Vec<Tile>,
}

impl Global {
    // Read a global + local from a Readable
    fn read<R: Read + BufRead>(reader: &mut R) -> Result<(Global, Local)> {
        let mut tiles = Vec::new();
        let mut width = 0;
        let mut height = 0;

        let mut snakes = vec![
            vec![], // 0..=9
            vec![], // a..=z
            vec![], // A..=Z
        ];

        let mut fruit = Vec::new();

        for (y, line) in reader.lines().enumerate() {
            let line = line?;

            if width == 0 {
                width = line.len();
            } else if line.len() != width {
                return Err(anyhow!("Inconsistent line length"));
            }

            height += 1;

            for (x, c) in line.chars().enumerate() {
                // Default to empty to put tiles under the snakes
                let mut tile = Tile::Empty;

                match c {
                    '-' => tile = Tile::Empty,
                    '#' => tile = Tile::Wall,
                    '=' => tile = Tile::Exit,
                    '+' => fruit.push(Point { x, y }),
                    '0'..='9' => snakes[0].push((c, Point { x, y })),
                    'a'..='z' => snakes[1].push((c, Point { x, y })),
                    'A'..='Z' => snakes[2].push((c, Point { x, y })),
                    _ => return Err(anyhow!("Invalid character {c} in map")),
                }

                tiles.push(tile);
            }
        }

        dbg!(width * height, tiles.len());

        // Remove empty snakes, sort points by character, convert to snake struct
        snakes.retain(|snake| !snake.is_empty());
        let snakes = snakes.iter().map(|points| {
            let mut points = points.clone();
            points.sort();
            Snake { 
                head: points.first().unwrap().0,
                points: points.iter().map(|(_, point)| *point).collect::<Vec<_>>(),
            }
        }).collect::<Vec<_>>();
        let local = Local { snakes, fruit };

        Ok((Global { width, height, tiles }, local))
    }

    // Get the tile at a point
    fn tile(&self, point: Point) -> Tile {
        self.tiles[point.y * self.width + point.x]
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Point {
    x: usize,
    y: usize,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Step {
    snake: usize,
    direction: Direction,
}

impl State<Global, Step> for Local {
    fn is_valid(&self, _global: &Global) -> bool {
        true
    }

    fn is_solved(&self, _map: &Global) -> bool {
        // All snakes exited and fruit eaten
        self.snakes.is_empty() && self.fruit.is_empty()
    }

    fn next_states(&self, global: &Global) -> Option<Vec<(i64, Step, Local)>> {
        let mut next_states = Vec::new();

        let is_occupied = |point: Point| {
            global.tile(point) == Tile::Wall || self.snakes.iter().any(|snake| snake.points.contains(&point))
        };

        for snake_index in 0..=self.snakes.len() {
            'one_direction: for direction in [Direction::Up, Direction::Down, Direction::Left, Direction::Right].iter() {
                let step = Step { snake: snake_index, direction: *direction };
                let mut new_local = self.clone();

                if let Some(snake) = new_local.snakes.get(snake_index) {
                    let head = snake.points.first().unwrap();

                    if head.x == 0 && direction == &Direction::Left {
                        continue 'one_direction;
                    } else if head.y == 0 && direction == &Direction::Up {
                        continue 'one_direction;
                    } else if head.x == global.width - 1 && direction == &Direction::Right {
                        continue 'one_direction;
                    } else if head.y == global.height - 1 && direction == &Direction::Down {
                        continue 'one_direction;
                    }

                    let new_head = match direction {
                        Direction::Up => Point { x: head.x, y: head.y - 1 },
                        Direction::Down => Point { x: head.x, y: head.y + 1 },
                        Direction::Left => Point { x: head.x - 1, y: head.y },
                        Direction::Right => Point { x: head.x + 1, y: head.y },
                    };

                    // Cannot move into a wall or another snake
                    if is_occupied(new_head) {
                        continue 'one_direction;
                    }

                    // Check for exit
                    if global.tile(new_head) == Tile::Exit {
                        // Cannot exit until all fruit is eaten
                        if !new_local.fruit.is_empty() {
                            continue 'one_direction;
                        }

                        // Otherwise, snake literally exits
                        // Bye bye snake
                        new_local.snakes.remove(snake_index);
                    } else {
                        // Update the snake
                        if let Some(snake) = new_local.snakes.get_mut(snake_index) {
                            snake.points.insert(0, new_head);

                            // Potentially eat fruit; don't remove tail if fruit was eaten
                            if let Some(fruit_index) = new_local.fruit.iter().position(|fruit| fruit == &new_head) {
                                new_local.fruit.remove(fruit_index);
                            } else {
                                snake.points.pop();
                            }
                        }

                        // Fall
                        loop {
                            let on_bottom = new_local.snakes.iter().any(|snake| {
                                let head = snake.points.first().unwrap();
                                head.y == global.height - 1
                            });
                            if on_bottom {
                                todo!("handle snakes falling out of the world");
                            }

                            let can_fall = new_local.snakes.iter().all(|snake| {
                                let head = snake.points.first().unwrap();
                                let below = Point { x: head.x, y: head.y + 1 };
                                !is_occupied(below)
                            });

                            if !can_fall {
                                break;
                            }

                            for snake in new_local.snakes.iter_mut() {
                                for point in snake.points.iter_mut() {
                                    point.y += 1;
                                }
                            }
                        }
                    }

                    next_states.push((1, step, new_local));
                }
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
            for x in 0..global.width {
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
