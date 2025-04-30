use direction::Direction;
use itertools::{Either, Itertools};

use crate::{local::Local, Global};

impl Local {
    #[tracing::instrument(skip(self, global))]
    pub(crate) fn step(&self, d: Direction, global: &Global) -> Result<Self, String> {
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
