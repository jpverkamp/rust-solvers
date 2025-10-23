use direction::Direction;
use point::Point;
use solver::State;

use crate::{
    Global,
    global::{Critter, CritterKind, Tile},
    local::Local,
};

type Step = (usize, Direction);

impl State<Global, (usize, Direction)> for Local {
    #[tracing::instrument(skip(_global), ret)]
    fn is_valid(&self, _global: &Global) -> bool {
        // TODO
        true
    }

    #[tracing::instrument(skip(global), ret)]
    fn is_solved(&self, global: &Global) -> bool {
        // All nests have a matching penguin on them
        for x in 0..global.width {
            for y in 0..global.height {
                if let Tile::Nest(nest_color) = global.tile_at((x, y).into())
                    && !self.critters.iter().any(|c| {
                        c.kind == CritterKind::Penguin
                            && c.location == (x, y).into()
                            && c.color == nest_color
                    })
                {
                    return false;
                }
            }
        }

        // There can't be any seals left (they all have to leave the level)
        if self.critters.iter().any(|c| c.kind == CritterKind::Seal) {
            return false;
        }

        true
    }

    #[tracing::instrument(skip(global))]
    fn next_states(&self, global: &Global) -> Option<Vec<(i64, Step, Local)>> {
        let mut next_states = vec![];

        for (index, _critter) in self.critters.iter().enumerate() {
            for d in Direction::all() {
                if let Some(move_to) = self.try_move(global, index, d) {
                    let mut new_local = self.clone();
                    
                    if global.tile_at(move_to) == Tile::Water {
                        // The critter went swimming
                        new_local.critters.remove(index);
                    } else {
                        new_local.critters[index].location = move_to;
                    }
                    
                    next_states.push((1, (index, d), new_local))
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

    fn heuristic(&self, _global: &Global) -> i64 {
        // TODO
        0
    }

    fn stringify(&self, global: &Global) -> String {
        let mut result = String::new();

        let mut critters_to_print = vec![];

        for row in 0..global.height {
            // The horizontal walls
            for col in 0..global.width {
                let p: Point = (col, row).into();

                result.push(' ');
                result.push(global.wall_at(p, Direction::Up).as_horizontal_char());
            }
            result.push('\n');

            // The vertical walls and tiles
            for col in 0..global.width {
                let p: Point = (col, row).into();

                result.push(global.wall_at(p, Direction::Left).as_vertical_char());

                if let Some(critter) = self
                    .critters
                    .iter()
                    .find(|&c| c.location.x == col as isize && c.location.y == row as isize)
                {
                    let c = Critter::index_char(critters_to_print.len());
                    result.push(c);
                    critters_to_print.push((c, critter));
                } else {
                    result.push(global.tile_at(p).into())
                }
            }

            // The last vertical wall
            result.push(
                global
                    .wall_at((global.width - 1, row).into(), Direction::Right)
                    .as_vertical_char(),
            );
            result.push('\n');
        }

        // The last row of horizontal walls!
        for col in 0..global.width {
            result.push(' ');
            result.push(
                global
                    .wall_at((col, global.height - 1).into(), Direction::Down)
                    .as_horizontal_char(),
            );
        }
        result.push('\n');
        result.push('\n');
        for (c, critter) in critters_to_print {
            let Critter { kind, color, .. } = critter;
            result.push_str(format!("{c}: {color:?} {kind:?}\n").as_str());
        }

        result
    }
}
