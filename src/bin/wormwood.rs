use std::io::Read;

use itertools::{Either, Itertools};
use solver::{Direction, Point, Solver, State};

#[derive(Debug, Clone, Default)]
struct Global {
    // Map settings
    width: isize,
    height: isize,
    cells: Vec<bool>,
}

impl From<&str> for Global {
    fn from(input: &str) -> Self {
        let mut global = Global::default();

        for line in input.lines() {
            let chars = line.chars();
            let mut width = 0;
            global.height += 1;

            for c in chars {
                width += 1;
                if c == '#' {
                    global.cells.push(true);
                } else {
                    global.cells.push(false);
                }
            }

            if global.width == 0 {
                global.width = width;
            } else if global.width != width {
                panic!("Map width mismatch: {} != {}", global.width, width);
            }
        }

        global
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Local {
    worm: Vec<Point>,
    blocks: Vec<Vec<Point>>,
}

impl Local {
    #[tracing::instrument(skip(self, global))]
    fn step(&self, d: Direction, global: &Global) -> Result<Self, String> {
        // The new head of the worm cannot go more than worm length out of bounds
        // And cannot eat through the floor at all
        let new_head = self.worm[0] + d.into();
        let worm_len = self.worm.len() as isize;

        if new_head.x < -worm_len
            || new_head.x >= global.width + worm_len
            || new_head.y < -worm_len
            || new_head.y >= global.height
        {
            return Err("Out of bounds".to_string());
        }

        // Cannot double back
        if self.worm.contains(&new_head) {
            return Err("Cannot double back".to_string());
        }

        // Otherwise, update state:
        // Worms moves, active cells are eaten, gravity is applied
        let mut new_blocks = self.blocks.clone();
        let mut new_worm = vec![new_head];

        // Move the worm
        new_worm.extend(self.worm.iter().take(self.worm.len() - 1).cloned());

        // Eat a chunk from a block (if possible)
        let mut eat_index = None;
        for (i, b) in new_blocks.iter_mut().enumerate() {
            if b.contains(&new_head) {
                eat_index = Some(i);
                b.retain(|p| p != &new_head);
                break;
            }
        }

        // If we did, potentially split that block
        if let Some(i) = eat_index {
            let b = &new_blocks[i];

            if b.is_empty() {
                // If it was empty, just remove it
                new_blocks.remove(i);
            } else {
                // If you can't reach all points from any remaining point, split the block
                let mut visited = vec![false; b.len()];
                let mut stack = vec![b[0]];
                while let Some(p) = stack.pop() {
                    if let Some(i) = b.iter().position(|x| x == &p) {
                        if !visited[i] {
                            visited[i] = true;
                            for n in p.neighbors() {
                                if n != new_head
                                    && b.contains(&n)
                                    && !visited[b.iter().position(|x| x == &n).unwrap()]
                                {
                                    stack.push(n);
                                }
                            }
                        }
                    }
                }

                if visited.iter().any(|&x| !x) {
                    let (b1, b2) = b
                        .iter()
                        .enumerate()
                        .filter(|(_, p)| new_head != **p)
                        .partition_map(|(i, p)| {
                            if visited[i] {
                                Either::Left(*p)
                            } else {
                                Either::Right(*p)
                            }
                        });

                    new_blocks.remove(i);
                    new_blocks.push(b1);
                    new_blocks.push(b2);
                }
            }
        }

        'falling: loop {
            // Apply gravity to worm
            'worm_falling: {
                // If any point is directly above a block or the ground, supported
                if new_worm.iter().any(|p| {
                    let down = *p + Direction::Down.into();
                    down.y >= global.height || new_blocks.iter().any(|b| b.contains(&down))
                }) {
                    tracing::debug!("Worm is supported directly");
                    break 'worm_falling;
                }

                // If two adjacent segments are supported on the sides, supported
                // TODO: Is this the actual condition?
                if new_worm.windows(2).any(|w| {
                    let left0 = w[0] + Direction::Left.into();
                    let right0 = w[0] + Direction::Right.into();
                    let left1 = w[1] + Direction::Left.into();
                    let right1 = w[1] + Direction::Right.into();

                    (new_blocks.iter().any(|b| b.contains(&left0))
                        || new_blocks.iter().any(|b| b.contains(&right0)))
                        && (new_blocks.iter().any(|b| b.contains(&left1))
                            || new_blocks.iter().any(|b| b.contains(&right1)))
                }) {
                    tracing::debug!("Worm is supported on the sides");
                    break 'worm_falling;
                }

                // If we passed all other conditions, still falling, update worm and continue loop
                tracing::debug!("Worm is falling");
                new_worm
                    .iter_mut()
                    .for_each(|p| *p = *p + Direction::Down.into());

                continue 'falling;
            }

            // Apply gravity to blocks
            for (i, b) in new_blocks.iter().enumerate() {
                'block_falling: {
                    // Supported by the ground
                    if b.iter().any(|p| p.y >= global.height - 1) {
                        tracing::debug!("Block {i} is supported by the ground");
                        break 'block_falling;
                    }

                    // Supported by the worm
                    if b.iter()
                        .any(|p| new_worm.contains(&(*p + Direction::Down.into())))
                    {
                        tracing::debug!("Block {i} is supported by the worm");
                        break 'block_falling;
                    }

                    // Supported by another block
                    if let Some((j, _)) = new_blocks.iter().enumerate().find(|(j, b2)| {
                        i != *j
                            && b.iter().any(|p| {
                                let down = *p + Direction::Down.into();
                                b2.contains(&down)
                            })
                    }) {
                        tracing::debug!("Block {i} is supported by block {j}");
                        break 'block_falling;
                    }

                    // If the worm was solely supported by this block by the sides, it falls too
                    // If it was supported under by this block this will be handled next loop
                    if new_worm.windows(2).any(|w| {
                        let left0 = w[0] + Direction::Left.into();
                        let right0 = w[0] + Direction::Right.into();
                        let left1 = w[1] + Direction::Left.into();
                        let right1 = w[1] + Direction::Right.into();

                        (b.contains(&left0) || b.contains(&right0))
                            && (b.contains(&left1) || b.contains(&right1))
                    }) {
                        tracing::debug!("Worm might be slide falling due to block {i}");

                        if new_worm.iter().any(|p| {
                            let down = *p + Direction::Down.into();
                            self.blocks.iter().any(|b| b.contains(&down))
                        }) {
                            // Another block is directly supporting the worm
                            tracing::debug!("^ Just kidding, worm is directly supported");
                        } else if false {
                            // TODO: Another block is side supporting the worm
                        } else {
                            // Otherwise, fall with the block
                            tracing::debug!("^ Yes, it's falling");
                            new_worm
                                .iter_mut()
                                .for_each(|p| *p = *p + Direction::Down.into());
                        }
                    }

                    // If we made it to this point, the block is falling
                    // Update it and continue falling
                    // This may have allowed the worm to fall
                    tracing::debug!("Block {i} ({b:?}) is falling");
                    for p in new_blocks[i].iter_mut() {
                        *p = *p + Direction::Down.into();
                    }

                    continue 'falling;
                }
            }

            // If we made it out of both loops, we're done with gravity
            break 'falling;
        }

        // If we made it all the way here, we have a valid state
        Ok(Local {
            worm: new_worm,
            blocks: new_blocks,
        })
    }
}

impl std::fmt::Display for Local {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut result = String::new();
        result.push_str("Worm{");
        for (i, p) in self.worm.iter().enumerate() {
            if i > 0 {
                result.push(' ');
            }
            result.push_str(&format!("{},{}", p.x, p.y));
        }
        result.push_str("} ");

        result.push_str("Blocks{");
        for b in self.blocks.iter() {
            result.push_str(&format!(
                "[{}]",
                b.iter()
                    .map(|p| format!("{},{}", p.x, p.y))
                    .collect::<Vec<_>>()
                    .join(" ")
            ));
        }

        write!(f, "{}", result)
    }
}

impl Global {
    fn make_local(&self) -> Local {
        // Start at the lower left 3 in length
        let worm = vec![
            Point {
                x: -1,
                y: self.height - 1,
            },
            Point {
                x: -2,
                y: self.height - 1,
            },
            Point {
                x: -3,
                y: self.height - 1,
            },
        ];

        // Create a single initial block
        let mut block = vec![];
        for y in 0..self.height {
            for x in 0..self.width {
                block.push(Point { x, y });
            }
        }
        let blocks = vec![block];

        Local { worm, blocks }
    }
}

impl State<Global, Direction> for Local {
    #[tracing::instrument(skip(global), ret)]
    fn is_valid(&self, global: &Global) -> bool {
        // If we don't have enough blocks left to match the global
        // if self.blocks.iter().map(|b| b.len()).sum::<usize>()
        //     < global.cells.iter().filter(|&&x| x).count()
        // {
        //     tracing::debug!(
        //         "Not enough blocks overall: {} < {}",
        //         self.blocks.iter().map(|b| b.len()).sum::<usize>(),
        //         global.cells.iter().filter(|&&x| x).count()
        //     );
        //     return false;
        // }

        // Blocks can only fall straight down, if any column doesn't have enough
        // This is measured top down since chunks can never go up or side to side
        for x in 0..global.width {
            let mut chunk_count = 0;
            let mut target_count = 0;

            for y in 0..global.height {
                if self.blocks.iter().any(|b| b.contains(&Point { x, y })) {
                    chunk_count += 1;
                }

                if global.cells[(y * global.width + x) as usize] {
                    target_count += 1;
                }

                if chunk_count < target_count {
                    tracing::debug!(
                        "Not enough blocks in column {}: {} < {}",
                        x,
                        chunk_count,
                        target_count
                    );
                    return false;
                }
            }
        }

        true
    }

    #[tracing::instrument(skip(global), ret)]
    fn is_solved(&self, global: &Global) -> bool {
        let mut local_cells = vec![false; global.width as usize * global.height as usize];

        for b in self.blocks.iter() {
            for p in b.iter() {
                let index = (p.y * global.width + p.x) as usize;
                local_cells[index] = true;
            }
        }

        global.cells == local_cells
    }

    #[tracing::instrument(skip(self, global), fields(self = %self))]
    fn next_states(&self, global: &Global) -> Option<Vec<(i64, Direction, Local)>> {
        let mut next_states = vec![];

        for d in Direction::all() {
            match self.step(d, global) {
                Ok(new_state) => next_states.push((1, d, new_state)),
                Err(e) => {
                    tracing::debug!("Invalid step in direction {:?}: {}", d, e);
                }
            }
        }

        // If we have any new states, return them
        if !next_states.is_empty() {
            Some(next_states)
        } else {
            None
        }
    }

    fn heuristic(&self, global: &Global) -> i64 {
        // The number of chunks not in the correct position
        let mut local_cells = vec![false; global.width as usize * global.height as usize];
        for b in self.blocks.iter() {
            for p in b.iter() {
                let index = (p.y * global.width + p.x) as usize;
                local_cells[index] = true;
            }
        }

        let mut heuristic = 0;
        for (i, &cell) in local_cells.iter().enumerate() {
            if cell != global.cells[i] {
                heuristic += 1;
            }
        }

        heuristic
    }

    fn stringify(&self, global: &Global) -> String {
        let worm_len = self.worm.len() as isize;
        let print_width = global.width + 2 * worm_len;
        let print_height = global.height + worm_len;

        let mut chars = vec![vec!['.'; print_width as usize]; print_height as usize];

        for (i, b) in self.blocks.iter().enumerate() {
            for p in b.iter() {
                chars[(p.y + 3) as usize][(p.x + 3) as usize] = (b'0' + (i as u8)) as char;
            }
        }

        for worm_p in self.worm.iter() {
            chars[(worm_p.y + 3) as usize][(worm_p.x + 3) as usize] = '*';
        }

        chars[(self.worm[0].y + 3) as usize][(self.worm[0].x + 3) as usize] = '☺';

        let mut result = String::new();
        result.push_str(&format!("{}x{}\n", global.width, global.height));
        for row in chars {
            result.push_str(&row.iter().collect::<String>());
            result.push('\n');
        }
        result
    }
}

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    let global = Global::from(input.as_str());
    let local = global.make_local();

    let tracing_enabled = std::env::var("WORMWOOD_TRACE").is_ok();

    if tracing_enabled {
        tracing_subscriber::fmt()
            .without_time()
            .with_max_level(tracing::Level::DEBUG)
            .init();
    } else {
        env_logger::init();
    }

    log::info!("Initial state:\n{}", local.stringify(&global));

    // If we have args, run each of those as a solution
    if std::env::args().len() > 1 {
        for sequence in std::env::args().skip(1) {
            let mut local = global.make_local();

            for d in sequence.chars() {
                let d = match d {
                    'U' | 'u' => Direction::Up,
                    'D' | 'd' => Direction::Down,
                    'L' | 'l' => Direction::Left,
                    'R' | 'r' => Direction::Right,
                    ' ' => continue,
                    _ => panic!("Invalid direction: {}", d),
                };

                match local.step(d, &global) {
                    Ok(new_state) => {
                        log::debug!("Step: {d:?}");
                        log::debug!("{}", new_state.stringify(&global));
                        log::debug!("{new_state}");
                        local = new_state;
                    }
                    Err(e) => {
                        log::error!("Invalid step in direction {:?}: {}", d, e);
                        break;
                    }
                }
            }
        }

        return;
    }

    // Otherwise, run the solver
    let mut solver = Solver::new(global.clone(), local.clone());

    while let Some(state) = solver.next() {
        if solver.states_checked() % 100_000 == 0 {
            log::debug!("\n{}", state.stringify(&global));
            log::debug!("{solver}");
        }
    }
    let solution = solver.get_solution();

    if let Some(solution) = solution {
        log::info!("{}", solver.stringify(&solution));

        let path = solver.path(&local, &solution).unwrap();
        log::info!("Path: {:?}", path);

        let sequence = path
            .iter()
            .map(|d| match d {
                Direction::Up => 'U',
                Direction::Down => 'D',
                Direction::Left => 'L',
                Direction::Right => 'R',
            })
            .collect::<String>();

        let sequence = sequence
            .chars()
            .chunks(5)
            .into_iter()
            .map(|chunk| chunk.collect::<String>())
            .collect::<Vec<_>>()
            .join(" ");

        println!("Sequence: {}", sequence);
    } else {
        println!("No solution found");
        std::process::exit(1);
    }
}
