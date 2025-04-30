use direction::Direction;
use point::Point;
use solver::State;

use crate::{local::Local, Global};

impl State<Global, Direction> for Local {
    #[tracing::instrument(skip(global), ret)]
    fn is_valid(&self, global: &Global) -> bool {
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
            let _p = Point {
                x: (i % global.width as usize) as isize,
                y: (i / global.width as usize) as isize,
            };
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
