use direction::Direction;
use point::Point;
use solver::State;

type Global = ();

#[derive(Debug, Clone, Copy)]
pub(crate) enum Step {
    Move {
        critter_index: usize,
        direction: Direction,
    },
}

use crate::model::color::Color;
use crate::model::critter::CritterKind;
use crate::model::map::Map;
use crate::model::thing::{Thing, ThingKind};
use crate::model::tile::Tile;
use crate::model::toggle::ToggleKind;
use crate::model::wall::WallKind;

impl Map {
    // Reset any local state
    fn reset(&mut self) {
        self.teleport_cooldown = true;
        self.used_teleports.clear();
    }

    // Try to move the active critter in the given direction
    // Returns if the move is possible (and something actually changed)
    #[tracing::instrument(skip(self), ret, fields(critter = %self.critters[critter_index]))]
    pub(crate) fn try_move(
        &mut self,
        critter_index: usize,
        direction: Direction,
        first_call: bool,
    ) -> bool {
        // The first time try_move is called, we have to reset local state
        if first_call {
            self.reset();
        }

        // Do not move escaped critters
        if self.critters[critter_index].escaped() {
            return false;
        }

        // If the current critter is on a dust cloud it cannot move
        if first_call
            && matches!(
                self.tile_at(self.critters[critter_index].location()),
                Tile::Dust | Tile::Nest { dusty: true, .. }
            )
        {
            return false;
        }

        // Take steps until we should not move any more
        let original_state = self.clone(); // TODO: Expensive...
        while self.try_move_one(critter_index, direction, 0, false) {}

        // HACK: Escape any animals that ended up on water or cracked floor after movement
        // I'm not sure how we ended up in this state; it's probably crutches returning early
        for index in 0..self.critters.len() {
            if self.critters[index].escaped() {
                continue;
            }

            match self.tile_at(self.critters[index].location()) {
                Tile::Water | Tile::CrackedFloor => {
                    tracing::debug!(
                        "Critter {:?} ended up on {:?}, escaping",
                        self.critters[index],
                        self.tile_at(self.critters[index].location())
                    );
                    self.critters[index].escape();
                }
                _ => {}
            }
        }

        // If nothing changed, this is invalid location
        if self == &original_state {
            return false;
        }

        true
    }

    // Internal function to move a single tile in a direction, looped to slide or used once to bounce
    // Modifies the map in place
    // Returns if we should continue moving
    #[tracing::instrument(skip(self), ret, fields(critter = %self.critters[critter_index]))]
    fn try_move_one(
        &mut self,
        critter_index: usize,
        direction: Direction,
        depth: usize,
        ignore_water: bool,
    ) -> bool {
        // If we're stuck in a bouncing loop, launch off the map
        // TODO: Magick constants!
        if depth > 10 {
            self.critters[critter_index].escape();
            return false;
        }

        let me = self.critters[critter_index];

        match self.tile_at(me.location()) {
            Tile::Water => {
                // If we're on water, don't move
                if ignore_water {
                    // Picked up a hammer on a cracked floor that broke
                    // Comes up in 025
                    tracing::debug!("on water, but just sliiiiiding on by");
                } else {
                    tracing::debug!("escaping into the water at {:?}", me.location());
                    self.critters[critter_index].escape();
                    return false;
                }
            }
            Tile::CrackedFloor => {
                // Cracked tiles turn into water
                // But we're allowed to continue (will stop if we hit it again)
                tracing::debug!("broke the floor");
                self.break_floor(me.location());
            }
            Tile::Nest { .. } | Tile::Floor | Tile::Dust => {
                // Everything else just keep on sliding
            }
            Tile::Wall => {
                unreachable!(
                    "Wall tiles should always be surrounded, so this should be impossible"
                );
            }
            Tile::Teleport(_) => {
                // Handle below
            }
            Tile::Toggle(_) => {
                // Handle below
            }
        }

        match self.maybe_do_teleport(critter_index, direction) {
            Some(end_movement) => return end_movement,
            None => {
                // Didn't teleport
            }
        }

        let wall = self.wall_at(me.location(), direction);
        match wall {
            WallKind::Empty => {}
            WallKind::Solid => {
                if self.critters[critter_index].carrying() == Some(ThingKind::Spring) {
                    tracing::debug!("bounced off a wall");
                    self.try_move_one(critter_index, direction.flip(), depth + 1, true);
                } else {
                    tracing::debug!("hit wall");
                }

                // Either way, don't keep moving
                return false;
            }
            WallKind::Color(c) => {
                if c == me.color() {
                    // Go right through my own colored walls!
                } else {
                    // Treat every other color as solid
                    if self.critters[critter_index].carrying() == Some(ThingKind::Spring) {
                        tracing::debug!("bounced off a mis-matched colored wall");
                        self.try_move_one(critter_index, direction.flip(), depth + 1, true);
                        self.maybe_do_teleport(critter_index, direction.flip());
                    } else {
                        tracing::debug!("hit colored wall");
                    }

                    // Either way, don't keep moving
                    return false;
                }
            }
            WallKind::Cracked => {
                tracing::debug!("hit a cracked wall, breaking it");
                self.break_wall(me.location(), direction);

                match self.critters[critter_index].carrying() {
                    Some(ThingKind::Spring) => {
                        tracing::debug!("bounced off a wall");
                        self.try_move_one(critter_index, direction.flip(), depth + 1, false);
                        self.maybe_do_teleport(critter_index, direction.flip());
                        return false;
                    }
                    Some(ThingKind::Hammer) => {
                        tracing::debug!("smashed right on through it");
                    }
                    Some(ThingKind::Crutch) | Some(ThingKind::Bomb) | None => {
                        return false;
                    }
                }
            }
        }

        // Bumped into any other critter
        if let Some(other_critter) = self
            .critters
            .iter()
            .position(|c| c.location() == me.location() + direction.into())
        {
            // Behavior based on what we're carrying
            match self.critters[critter_index].carrying() {
                Some(ThingKind::Spring) => {
                    tracing::debug!("bounced off another critter");
                    self.try_move_one(critter_index, direction.flip(), depth + 1, false);
                }
                Some(ThingKind::Hammer) => {
                    tracing::debug!("hammered off another critter");

                    // If the other critter is standing on a teleport, we'll end up going through that first
                    // Which is a crazy edge case, first discovered in 041 Flow
                    // TODO: Handle this better (?)

                    // So:
                    // - Move the other critter out of the way (so we don't hit them again)
                    // - Move me into that location; might be able to teleport
                    // - Put them back (overlapping, it's fine)
                    // - Let them move now

                    let other_location = self.critters[other_critter].location();
                    self.critters[other_critter].move_to(Point { x: -10, y: -10 });

                    self.try_move_one(critter_index, direction, depth + 1, true);

                    self.teleport_cooldown = false;
                    self.maybe_do_teleport(critter_index, direction);

                    self.critters[other_critter].move_to(other_location);
                    if self.try_move(other_critter, direction, false) {
                        // The other could move, all is well
                    } else {
                        // The other couldn't move, remove it
                        self.critters[other_critter].escape();
                    }
                }
                Some(ThingKind::Crutch | ThingKind::Bomb) | None => {
                    tracing::debug!("hit another critter");
                }
            }

            // Behavior based on what they are carrying
            match self.critters[other_critter].carrying() {
                Some(ThingKind::Bomb) => {
                    tracing::debug!("other critter is carrying a bomb, BOOM");
                    self.critters[critter_index].escape();
                    self.critters[other_critter].escape();
                    return false;
                }
                Some(ThingKind::Spring | ThingKind::Hammer | ThingKind::Crutch) | None => {}
            }

            self.maybe_do_teleport(critter_index, direction);
            return false;
        }

        let dst = me.location() + direction.into();
        tracing::debug!(
            "{critter:?} moved to {dst:?}",
            critter = self.critters[critter_index]
        );
        self.critters[critter_index].move_to(dst);

        // Moved onto a thing, pick it up
        // If we were already holding something, chuck our current thing into the water
        if let Some(index) = self.things.iter().position(|t| t.location == dst) {
            let thing = self.things.remove(index);
            tracing::debug!("picked up {thing:?}");
            self.critters[critter_index].pick_up(thing.kind);
        }

        // If we moved onto a toggle, trigger any matching rules
        if let Tile::Toggle(c) = self.tile_at(dst) {
            tracing::debug!("stepped on toggle {c}, checking rules");

            let matching_rules: Vec<_> = self
                .toggle_rules
                .iter()
                .filter(|r| r.key == c)
                .map(|r| r.kind)
                .collect();

            for kind in matching_rules {
                match kind {
                    ToggleKind::Floor(p) => {
                        tracing::debug!("toggle {c} rule: toggling floor at {p:?}");
                        self.toggle_floor(p);
                    }
                    ToggleKind::Wall(p, d) => {
                        tracing::debug!(
                            "toggle {c} rule: toggling wall at {p:?} in direction {d:?}"
                        );
                        self.toggle_wall(p, d);
                    }
                }
            }
        }

        if self.teleport_cooldown {
            tracing::debug!("ending teleport cooldown");
            self.teleport_cooldown = false;
        }

        // If we're carrying a crutch, only move once
        if self.critters[critter_index].carrying() == Some(ThingKind::Crutch) {
            self.maybe_do_teleport(critter_index, direction);
            return false;
        }

        true
    }

    fn maybe_do_teleport(&mut self, critter_index: usize, direction: Direction) -> Option<bool> {
        let me = self.critters[critter_index];
        if let Tile::Teleport(target) = self.tile_at(me.location()) {
            if self.teleport_cooldown {
                tracing::debug!("Cannot teleport, on cooldown");
                return None;
            }

            if self.used_teleports.contains(&(direction, me.location())) {
                tracing::debug!("teleport loop detected, not using teleport LAUNCHING");
                self.critters[critter_index].escape();
                return Some(false);
            }

            match self.critters.iter().position(|c| c.location() == target) {
                Some(other_critter) => {
                    tracing::debug!("teleporting to {target:?}, TELEFRAG");
                    self.critters[other_critter].escape();
                    return Some(true);
                }
                None => {
                    // Teleport there and keep going!
                    tracing::debug!("teleporting to {target:?}");
                    self.used_teleports.push((direction, me.location()));
                    self.teleport_cooldown = true;
                    self.critters[critter_index].move_to(target);
                    return Some(true);
                }
            }
        }

        None
    }
}

impl State<Global, Step> for Map {
    #[tracing::instrument(skip(self), ret)]
    fn is_valid(&self, _: &Global) -> bool {
        // TODO
        true
    }

    #[tracing::instrument(skip(self), ret)]
    fn is_solved(&self, _: &Global) -> bool {
        // All nests have a matching penguin on them
        for x in 0..self.width {
            for y in 0..self.height {
                if let Tile::Nest {
                    color: nest_color, ..
                } = self.tile_at((x, y).into())
                    && !self.critters.iter().any(|c| {
                        c.kind() == CritterKind::Penguin
                            && c.location() == (x, y).into()
                            && c.color() == nest_color
                    })
                {
                    tracing::debug!("Unsolved nest");
                    return false;
                }
            }
        }

        // There can't be any seals left (they all have to leave the level)
        if self
            .critters
            .iter()
            .any(|c| c.kind() == CritterKind::Seal && !c.escaped())
        {
            tracing::debug!("Unsolved seal");
            return false;
        }

        true
    }

    #[tracing::instrument(skip(self))]
    fn next_states(&self, _: &Global) -> Option<Vec<(i64, Step, Map)>> {
        let mut next_states = vec![];

        // Try moving the active critter
        for critter_index in 0..self.critters.len() {
            if self.critters[critter_index].escaped() {
                continue;
            }

            if self.critters[critter_index].color() == Color::Gray {
                continue;
            }

            if self.tile_at(self.critters[critter_index].location()) == Tile::Dust {
                continue;
            }

            for direction in Direction::all() {
                let mut new_map = self.clone();
                if new_map.try_move(critter_index, direction, true) {
                    next_states.push((
                        1,
                        Step::Move {
                            critter_index,
                            direction,
                        },
                        new_map,
                    ))
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
                    .find(|&c| c.location() == (col, row).into())
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
            let kind = critter.kind();
            let color = critter.color();
            let carrying = critter.carrying();
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
