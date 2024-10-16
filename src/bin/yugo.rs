use std::{io::Read, ops::AddAssign};

use solver::{Direction, Point, Solver, State};

#[derive(Debug, Clone)]
struct Global {
    width: isize,
    height: isize,
    walls: Vec<bool>,
}

impl Global {
    fn new(width: isize, height: isize, walls: Vec<Point>) -> Self {
        let walls = (0..width * height)
            .map(|i| {
                let x = i % width;
                let y = i / width;
                walls
                    .iter()
                    .any(|&Point { x: wx, y: wy }| x == wx && y == wy)
            })
            .collect();

        Global {
            width,
            height,
            walls,
        }
    }

    fn is_wall(&self, pt: Point) -> bool {
        if pt.x < 0 || pt.y < 0 || pt.x >= self.width || pt.y >= self.height {
            return true;
        }

        let index = (pt.y * self.width + pt.x) as usize;
        self.walls[index]
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Color {
    Red,
    Green,
    Blue,
    Yellow,
}

impl Color {
    fn all() -> Vec<Color> {
        vec![Color::Red, Color::Green, Color::Blue, Color::Yellow]
    }
}

impl TryFrom<char> for Color {
    type Error = ();

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value.to_ascii_uppercase() {
            'R' => Ok(Color::Red),
            'G' => Ok(Color::Green),
            'B' => Ok(Color::Blue),
            'Y' => Ok(Color::Yellow),
            _ => Err(()),
        }
    }
}

impl From<Color> for char {
    fn from(value: Color) -> Self {
        match value {
            Color::Red => 'R',
            Color::Green => 'G',
            Color::Blue => 'B',
            Color::Yellow => 'Y',
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Block {
    color: Color,
    points: Vec<Point>,
}

impl AddAssign<Direction> for Block {
    fn add_assign(&mut self, direction: Direction) {
        for point in self.points.iter_mut() {
            *point = *point + direction.into();
        }
    }
}

impl Block {
    fn overlaps(&self, other: &Block) -> bool {
        self.points.iter().any(|point| other.points.contains(point))
    }

    fn color_matches(&self, other: &Block) -> bool {
        self.color == other.color
    }

    fn adjacent_to(&self, other: &Block) -> bool {
        self.points.iter().any(|point| {
            other
                .points
                .iter()
                .any(|other_point| point.manhattan_distance(*other_point) == 1)
        })
    }

    fn can_merge(&self, other: &Block) -> bool {
        self.color_matches(other) && self.adjacent_to(other)
    }

    fn pushes(&self, other_block: &Block, direction: Direction) -> bool {
        self.points.iter().any(|point| {
            other_block
                .points
                .iter()
                .any(|other_point| *point + direction.into() == *other_point)
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Local {
    blocks: Vec<Block>,
}

impl Local {
    fn new(blocks: Vec<Block>) -> Self {
        Local { blocks }
    }

    fn is_empty(&self, point: Point) -> bool {
        self.blocks
            .iter()
            .all(|block| !block.points.contains(&point))
    }

    // Returns None if cannot move
    // Returns a list of block IDs that need to move otherwise
    fn can_move_block(
        &self,
        global: &Global,
        block_id: usize,
        direction: Direction,
    ) -> Option<Vec<usize>> {
        let mut moving_blocks = vec![block_id];

        // Cannot move into walls
        if self.blocks[block_id]
            .points
            .iter()
            .any(|point| global.is_wall(*point + direction.into()))
        {
            return None;
        }

        // If we move into a block, try to move that block the same direction first
        for (other_block_id, other_block) in self.blocks.iter().enumerate() {
            if block_id == other_block_id {
                continue;
            }

            if self.blocks[block_id].pushes(other_block, direction) {
                // Check if we can move the other block
                if let Some(recursive_moving_blocks) =
                    self.can_move_block(global, other_block_id, direction)
                {
                    // TODO: This is probably wrong, I think we just don't want the ID twice
                    if recursive_moving_blocks.contains(&block_id) {
                        return None;
                    }

                    moving_blocks.extend(recursive_moving_blocks);
                } else {
                    // Recursive move failed, so we can't move
                    return None;
                }
            }
        }

        // If we've made it this far, the move is valid
        Some(moving_blocks)
    }

    fn move_block(&mut self, block_id: usize, direction: Direction) {
        self.blocks[block_id] += direction.into();
    }
    
    fn merge(&mut self) {
        loop {
            let mut to_merge = None;

            'find_a_friend: for block_id in 0..self.blocks.len() {
                for other_block_id in (block_id + 1)..self.blocks.len() {
                    if self.blocks[block_id]
                        .can_merge(&self.blocks[other_block_id])
                    {
                        to_merge = Some((block_id, other_block_id));
                        break 'find_a_friend;
                    }
                }
            }

            // If we found a pair that can merge do it
            // Because of the order of the loops, block_id < other_block_id
            if let Some((block_id, other_block_id)) = to_merge {
                let mut merged_block = self.blocks[block_id].clone();
                merged_block
                    .points
                    .extend(self.blocks[other_block_id].points.iter());
                self.blocks.remove(other_block_id);
                self.blocks[block_id] = merged_block;
            } else {
                break;
            }
        }

    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Step {
    color: Color,
    index_of_that_color: usize,
    point: Point,
    direction: Direction,
}

impl State<Global, Step> for Local {
    fn is_valid(&self, _global: &Global) -> bool {
        true
    }

    fn is_solved(&self, _map: &Global) -> bool {
        // No more than one of each block of each color
        Color::all().iter().all(|color| {
            self.blocks
                .iter()
                .filter(|block| block.color == *color)
                .count()
                <= 1
        })
    }

    fn next_states(&self, global: &Global) -> Option<Vec<(i64, Step, Local)>> {
        let mut next_states = vec![];

        for block_id in 0..self.blocks.len() {
            // The 'nth' blue block for example
            let index_of_that_color = 1 + self.blocks
                .iter()
                .filter(|block| block.color == self.blocks[block_id].color)
                .position(|block| block.points == self.blocks[block_id].points)
                .unwrap();

            for direction in [Direction::Left, Direction::Right] {
                if let Some(moving_blocks) = self.can_move_block(global, block_id, direction) {
                    // Move first
                    let mut next_local = self.clone();
                    for moving_block in moving_blocks {
                        next_local.move_block(moving_block, direction);
                    }

                    // // TODO: Temporary sanity check, no point appears more than once
                    // assert!(next_local.blocks.iter().all(|block| {
                    //     block.points.iter().all(|point| {
                    //         next_local
                    //             .blocks
                    //             .iter()
                    //             .filter(|other_block| other_block.points.contains(point))
                    //             .count()
                    //             == 1
                    //     })
                    // }));

                    // Apply gravity
                    'still_falling: loop {
                        for block_id in 0..next_local.blocks.len() {
                            if let Some(falling_blocks) =
                                next_local.can_move_block(global, block_id, Direction::Down)
                            {
                                // println!("Falling blocks: {:?}", falling_blocks);
                                for falling_block in falling_blocks {
                                    next_local.move_block(falling_block, Direction::Down);
                                }
                                continue 'still_falling;
                            }
                        }
                        break;
                    }

                    // // TODO: Temporary sanity check, no point appears more than once
                    // assert!(next_local.blocks.iter().all(|block| {
                    //     block.points.iter().all(|point| {
                    //         next_local
                    //             .blocks
                    //             .iter()
                    //             .filter(|other_block| other_block.points.contains(point))
                    //             .count()
                    //             == 1
                    //     })
                    // }));

                    next_local.merge();

                    // // TODO: Temporary sanity check, no point appears more than once
                    // assert!(next_local.blocks.iter().all(|block| {
                    //     block.points.iter().all(|point| {
                    //         next_local
                    //             .blocks
                    //             .iter()
                    //             .filter(|other_block| other_block.points.contains(point))
                    //             .count()
                    //             == 1
                    //     })
                    // }));

                    let step = Step {
                        color: self.blocks[block_id].color,
                        index_of_that_color,
                        point: self.blocks[block_id].points[0],
                        direction,
                    };

                    next_states.push((1, step, next_local));
                }
            }
        }

        if next_states.is_empty() {
            None
        } else {
            Some(next_states)
        }
    }

    fn heuristic(&self, _global: &Global) -> i64 {
        0
    }

    fn stringify(&self, global: &Global) -> String {
        let mut chars = vec![vec![' '; global.width as usize]; global.height as usize];

        for x in 0..global.width {
            for y in 0..global.height {
                if global.is_wall(Point { x, y }) {
                    chars[y as usize][x as usize] = '█';
                }
            }
        }

        for block in &self.blocks {
            for Point { x, y } in &block.points {
                chars[*y as usize][*x as usize] = char::from(block.color);
            }
        }

        let mut result = String::new();
        for row in chars {
            result.push_str(&row.iter().collect::<String>());
            result.push('\n');
        }
        result
    }
}

fn load(input: &str) -> (Global, Local) {
    let mut walls = vec![];
    let mut blocks = vec![];
    let mut width = 0;
    let mut height = 0;

    for (y, line) in input.lines().enumerate() {
        height = height.max(1 + y as isize);

        for (x, c) in line.chars().enumerate() {
            width = width.max(1 + x as isize);

            if let Ok(c) = Color::try_from(c) {
                blocks.push(Block {
                    color: c,
                    points: vec![Point {
                        x: x as isize,
                        y: y as isize,
                    }],
                });
            } else if c == '#' {
                walls.push(Point {
                    x: x as isize,
                    y: y as isize,
                });
            }
        }
    }

    let global = Global::new(width, height, walls);
    let mut local = Local::new(blocks);
    local.merge();

    (global, local)
}

fn main() {
    env_logger::init();

    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    let (global, local) = load(&input);
    log::debug!("Initial state:\n{}", local.stringify(&global));

    let mut solver = Solver::new(global.clone(), local.clone());

    while let Some(state) = solver.next() {
        if solver.states_checked() % 10000 == 0 {
            log::debug!("{}", state.stringify(&global));
            log::debug!("{solver}");
        }
    }
    let solution = solver.get_solution();

    if let Some(solution) = solution {
        println!("{}", solver.stringify(&solution));

        let path = solver.path(&local, &solution).unwrap();
        for Step { color, index_of_that_color, direction, .. } in path {
            println!("{color:?} {index_of_that_color} {direction:?}");
        }
    } else {
        println!("No solution found");
    }
}
