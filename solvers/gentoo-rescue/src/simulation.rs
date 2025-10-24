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

use crate::model::color::Color;
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

        // We will throw this away if it's invalid, but this is necessary to update cracked walls/floors
        let mut new_map = self.clone();
        while new_map.try_move_one(direction, 0) {
            // Keep on moving
            // It feels weird to have an empty loop
        }

        // If nothing changed, this is invalid location
        if self == &new_map {
            return None;
        }

        // If the critter is on water, remove it and choose a new active critter
        if new_map.tile_at(new_map.critters[new_map.active_critter].location) == Tile::Water {
            tracing::info!("critter ESCAPED into the water");
            new_map.critters.remove(self.active_critter);
            if new_map.active_critter >= new_map.critters.len() {
                new_map.active_critter = 0;
            }

            return Some((new_map, true));
        }

        // Otherwise, the critter just moved
        Some((new_map, false))
    }

    // Internal function to move a single tile in a direction, looped to slide or used once to bounce
    // Modifies the map in place
    // Returns if we should continue moving
    #[tracing::instrument(skip(self), ret, fields(pt = ?self.critters[self.active_critter].location))]
    fn try_move_one(&mut self, direction: Direction, depth: usize) -> bool {
        // If we're stuck in a bouncing loop, launch off the map 
        // TODO: Do we have to actually stop at a specific point or just 'off'?
        // TODO: Magick constants!
        if depth > 10 { 
            self.critters[self.active_critter].location = Point { x: -10, y: -10 };
            return false;
        }

        let me = self.critters[self.active_critter];

        match self.tile_at(me.location) {
            Tile::Water => {
                // If we're on water, don't move
                tracing::debug!("stopped at water");
                return false;
            }
            Tile::CrackedFloor => {
                // Cracked tiles turn into water
                // But we're allowed to continue (will stop if we hit it again)
                tracing::debug!("broke the floor");
                self.break_floor(me.location);
            }
            Tile::Nest(_) | Tile::Floor => {
                // Everything else just keep on sliding
            },
            Tile::Wall => {
                unreachable!("Wall tiles should always be surrounded, so this should be impossible");
            }
        }

        // Standing on a thing, pick it up
        // TODO: Only if not carrying something, is this correct?
        if self.critters[self.active_critter].carrying.is_none()
            && let Some(index) = self.things.iter().position(|t| t.location == me.location)
        {
            let thing = self.things.remove(index);
            tracing::debug!("picked up {thing:?}");
            self.critters[self.active_critter].carrying = Some(thing.kind);
        }

        let wall = self.wall_at(me.location, direction);
        match wall {
            WallKind::Empty => {}
            WallKind::Solid | WallKind::Cracked => {
                if wall == WallKind::Cracked {
                    tracing::debug!("hit a cracked wall, breaking it");
                    self.break_wall(me.location, direction);
                }

                if self.critters[self.active_critter].carrying == Some(ThingKind::Spring) {
                    tracing::debug!("bounced off a wall");
                    self.try_move_one(direction.flip(), depth + 1);
                } else {
                    tracing::debug!("hit wall");    
                }

                // Either way, don't keep moving
                return false;
            }
        }

        // Bumped into any other critter
        // TODO: Do we bounce off critters? 
        if let Some(other_critter) = self
            .critters
            .iter()
            .position(|c| c.location == me.location + direction.into())
        {
            match self.critters[self.active_critter].carrying {
                Some(ThingKind::Spring) => {
                    tracing::debug!("bounced off another critter");
                    self.try_move_one(direction.flip(), depth + 1);
                },
                Some(ThingKind::Hammer) => {
                    tracing::debug!("hammered off another critter");
                    let my_index = self.active_critter;
                    let my_position = self.critters[my_index].location;

                    // The other critter gets bumped out of our way
                    // TODO: Handle recursion that moves the original critter out of the way
                    self.active_critter = other_critter;
                    match self.try_move(direction) {
                        Some((mut new_map, false)) => {
                            // The critter being bumped did not leave the map
                            std::mem::swap(self, &mut new_map);
                            self.active_critter = my_index;
                        },

                        Some((mut new_map, true)) => {
                            // The critter being bumped left the map
                            // This might have screwed up the active index, so (try to) find it again
                            std::mem::swap(self, &mut new_map);
                            match self.critters.iter().position(|oc| oc.location == my_position) {
                                Some(my_index) => self.active_critter = my_index,
                                None => panic!("Could not find original critter after hammer time")
                            }
                        },
                        None => {
                            self.active_critter = my_index;
                        }
                    }
                    
                    // We take that spot
                    self.try_move_one(direction, depth + 1);
                },
                None => {
                    tracing::debug!("hit another critter");
                },
            }
            return false;
        }

        let dst = me.location + direction.into();
        tracing::debug!("moved to {dst:?}");
        self.critters[self.active_critter].location = dst;
        true
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
        // Gray critters don't move under our control
        for i in 0..self.critters.len() {
            if i == self.active_critter {
                continue;
            }

            if self.critters[i].color == Color::Gray {
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
