use direction::Direction;
use point::Point;
use solver::State;

type Global = ();

#[derive(Debug, Clone, Copy)]
pub(crate) enum Step {
    SwitchCritter {
        critter: Critter,
    },
    Move {
        direction: Direction,
        new_critter: Option<Critter>,
    },
}

use crate::model::critter::{Critter, CritterKind};
use crate::model::map::Map;
use crate::model::thing::{Thing, ThingKind};
use crate::model::tile::Tile;
use crate::model::wall::WallKind;

impl Map {
    // Try to move the active critter in the given direction
    // Returns the point the critter moves to (if it moves) + if the critter changed
    #[tracing::instrument(skip(self), ret)]
    pub(crate) fn try_move(&self, direction: Direction) -> Option<(Map, bool)> {
        if self.critters.is_empty() {
            return None;
        }

        let mut pt = self.critters[self.active_critter].location;
        let mut moved = false;

        // We will throw this away if it's invalid, but this is necessary to update cracked walls/floors
        let mut new_map = self.clone();

        loop {
            match new_map.tile_at(pt) {
                Tile::Water => {
                    // If we're on water, stop moving
                    tracing::debug!("{pt:?} stopped at water");
                    break;
                }
                Tile::CrackedFloor => {
                    // Cracked tiles turn into water
                    // But we're allowed to continue (will stop if we hit it again)
                    tracing::debug!("{pt:?} broke the floor");
                    new_map.break_floor(pt);
                }
                Tile::Nest(_) | Tile::Floor => {
                    // Everything else just keep on sliding
                }
            }

            // Standing on a thing, pick it up
            // TODO: Only if not carrying something, is this correct?
            if new_map.critters[new_map.active_critter].carrying.is_none()
                && let Some(index) = self.things.iter().position(|t| t.location == pt)
            {
                let thing = new_map.things.remove(index);
                tracing::debug!("picked up {thing:?}");
                new_map.critters[new_map.active_critter].carrying = Some(thing.kind);
                moved = true;
            }

            match new_map.wall_at(pt, direction) {
                WallKind::Empty => {}
                WallKind::Solid => {
                    // Bumped into a wall
                    tracing::debug!("{pt:?} stopped at wall");
                    break;
                }
                WallKind::Cracked => {
                    // Bumped into a cracked wall, break it
                    // This counts as moving even even though we stopped
                    tracing::debug!("{pt:?} stopped at cracked wall, breaking it");
                    new_map.break_wall(pt, direction);
                    moved = true;
                    break;
                }
            }

            // Bumped into any other critter
            if self
                .critters
                .iter()
                .any(|c| c.location == pt + direction.into())
            {
                tracing::debug!("{pt:?} stopped at critter");
                break;
            }

            pt = pt + direction.into();
            tracing::debug!("moved to {pt:?}");
            moved = true;
        }

        // If, at the end of moving, the critter is carrying a spring, they bounce backwards one
        // TODO: Handle bouncing backwards over a wall
        // TODO: Is this handling of water correct?
        if new_map.critters[new_map.active_critter].carrying == Some(ThingKind::Spring)
            && new_map.tile_at(pt) != Tile::Water
        {
            tracing::debug!("bouncing backwards");
            pt = pt - direction.into();
            moved = true;
        }

        // If we didn't move, this is invalid location
        // TODO: Handle bouncing etc
        if !moved {
            return None;
        }

        // If the critter is on water, remove it and choose a new active critter
        if self.tile_at(pt) == Tile::Water {
            new_map.critters.remove(self.active_critter);
            if new_map.active_critter >= new_map.critters.len() {
                new_map.active_critter = 0;
            }

            return Some((new_map, true));
        }

        // Otherwise, the critter just moved
        new_map.critters[self.active_critter].location = pt;
        Some((new_map, false))
    }
}

impl State<Global, Step> for Map {
    #[tracing::instrument(skip(_global), ret)]
    fn is_valid(&self, _global: &Global) -> bool {
        // TODO
        true
    }

    #[tracing::instrument(ret)]
    fn is_solved(&self, _: &Global) -> bool {
        // All nests have a matching penguin on them
        for x in 0..self.width {
            for y in 0..self.height {
                if let Tile::Nest(nest_color) = self.tile_at((x, y).into())
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

    #[tracing::instrument()]
    fn next_states(&self, _: &Global) -> Option<Vec<(i64, Step, Map)>> {
        let mut next_states = vec![];

        // Try moving the active critter
        for d in Direction::all() {
            if let Some((new_map, new_critter)) = self.try_move(d) {
                // If the step resulted in a critter switch, record that in the step
                let step = Step::Move {
                    direction: d,
                    new_critter: if new_critter && !new_map.critters.is_empty() {
                        Some(new_map.critters[new_map.active_critter])
                    } else {
                        None
                    },
                };
                next_states.push((1, step, new_map))
            }
        }

        // Try switching to each other critter
        for i in 0..self.critters.len() {
            if i == self.active_critter {
                continue;
            }

            let mut new_map = self.clone();
            new_map.active_critter = i;
            let step = Step::SwitchCritter {
                critter: new_map.critters[i],
            };
            next_states.push((0, step, new_map));
        }

        // If we have any new states, return them
        if !next_states.is_empty() {
            Some(next_states)
        } else {
            None
        }
    }

    fn heuristic(&self, _: &Global) -> i64 {
        // TODO
        0
    }

    fn stringify(&self, _: &Global) -> String {
        let mut result = String::new();
        let mut critters_to_print = vec![];
        let mut things_to_print = vec![];
        let mut index = 0;

        fn index_char(index: usize) -> char {
            if index < 10 {
                (b'0' + (index as u8)) as char
            } else if index < 10 + 26 {
                (b'a' + ((index - 10) as u8)) as char
            } else if index < 10 + 26 + 26 {
                (b'A' + ((index - 10 - 26) as u8)) as char
            } else {
                unimplemented!("Too many critters!")
            }
        }

        for row in 0..self.height {
            // The horizontal walls
            for col in 0..self.width {
                let p: Point = (col, row).into();

                result.push(' ');
                result.push(self.wall_at(p, Direction::Up).as_horizontal_char());
            }
            result.push('\n');

            // The vertical walls and tiles
            for col in 0..self.width {
                let p: Point = (col, row).into();

                result.push(self.wall_at(p, Direction::Left).as_vertical_char());

                if let Some(critter) = self
                    .critters
                    .iter()
                    .find(|&c| c.location.x == col as isize && c.location.y == row as isize)
                {
                    let c = index_char(index);
                    index += 1;

                    result.push(c);
                    critters_to_print.push((c, critter));
                } else if let Some(thing) = self
                    .things
                    .iter()
                    .find(|&t| t.location.x == col as isize && t.location.y == row as isize)
                {
                    // TODO: Combine critters and things?

                    let c = index_char(index);
                    index += 1;

                    result.push(c);
                    things_to_print.push((c, thing));
                } else {
                    result.push(self.tile_at(p).into())
                }
            }

            // The last vertical wall
            result.push(
                self.wall_at((self.width - 1, row).into(), Direction::Right)
                    .as_vertical_char(),
            );
            result.push('\n');
        }

        // The last row of horizontal walls!
        for col in 0..self.width {
            result.push(' ');
            result.push(
                self.wall_at((col, self.height - 1).into(), Direction::Down)
                    .as_horizontal_char(),
            );
        }
        result.push('\n');
        result.push('\n');
        for (c, critter) in critters_to_print {
            let Critter {
                kind,
                color,
                carrying,
                ..
            } = critter;
            result.push_str(
                format!(
                    "{c}: {color:?} {kind:?}{carrying}\n",
                    carrying = if let Some(thing) = carrying {
                        format!(" w/{thing:?}")
                    } else {
                        String::new()
                    }
                )
                .as_str(),
            );
        }

        result.push('\n');
        for (c, thing) in things_to_print {
            let Thing { kind, .. } = thing;
            result.push_str(format!("{c}: {kind:?}\n").as_str());
        }

        result
    }
}
