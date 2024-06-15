use anyhow::{anyhow, Result};
use fxhash::FxHashMap;
use std::{
    io::{BufRead, Read},
    ops::{Add, Sub},
};

use solver::{Solver, State};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Tile {
    Empty,
    Wall,
    Exit,
    Spike,
    Portal,
}

#[derive(Debug, Clone)]
struct Global {
    width: usize,
    height: usize,
    tiles: FxHashMap<Point, Tile>,
    solutions: Vec<String>,
}

// Stores the minimum and maximum character for each snake type
// Third is true if the snake is actually a block
const SNAKE_COUNT: usize = 13;
const SNAKE_RANGES: [(char, char, bool); SNAKE_COUNT] = [
    ('0', '9', false),
    ('a', 'm', false),
    ('A', 'M', false),
    ('{', '~', false),
    ('W', 'W', true),
    ('X', 'X', true),
    ('Y', 'Y', true),
    ('Z', 'Z', true),
    ('w', 'w', true),
    ('x', 'x', true),
    ('y', 'y', true),
    ('z', 'z', true),
    ('?', '?', true),
];

impl Global {
    // Read a global + local from a Readable
    fn read<R: Read + BufRead>(reader: &mut R) -> Result<(Global, Local)> {
        let mut width = 0;
        let mut height = 0;

        let mut tiles = FxHashMap::default();
        let mut snakes = vec![vec![]; SNAKE_COUNT];
        let mut fruit = Vec::new();

        for (y, line) in reader.lines().enumerate() {
            let y = y as isize;
            let line = line?;

            if line.is_empty() {
                break;
            }

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
                    '!' => {
                        tiles.insert(p, Tile::Spike);
                    }
                    '@' => {
                        tiles.insert(p, Tile::Portal);
                    }
                    // Fruit is stored in the local state
                    '+' => {
                        fruit.push(p);
                    }
                    // Otherwise, probably a snake
                    _ => {
                        let (snake_type, _) = SNAKE_RANGES
                            .iter()
                            .enumerate()
                            .find(|(_, (min, max, _))| c >= *min && c <= *max)
                            .ok_or_else(|| anyhow!("Invalid character {c} in map"))?;

                        snakes[snake_type].push((c, p));
                    }
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

                let is_statue = SNAKE_RANGES
                    .iter()
                    .find_map(|(min, max, is_block)| {
                        if points[0].0 >= *min && points[0].0 <= *max {
                            Some(*is_block)
                        } else {
                            None
                        }
                    })
                    .expect("Invalid snake type");

                Snake {
                    head: points.first().unwrap().0,
                    points: points.iter().map(|(_, point)| *point).collect::<Vec<_>>(),
                    is_statue,
                }
            })
            .collect::<Vec<_>>();

        let local = Local { snakes, fruit };

        // Any remaining lines are solutions
        let mut solutions = Vec::new();
        for line in reader.lines() {
            let line = line?.trim().to_string();
            if line.is_empty() {
                break;
            }
            solutions.push(line);
        }

        Ok((
            Global {
                width,
                height,
                tiles,
                solutions,
            },
            local,
        ))
    }

    // Get the tile at a point
    fn tile(&self, point: Point) -> Tile {
        self.tiles.get(&point).copied().unwrap_or(Tile::Empty)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Point {
    x: isize,
    y: isize,
}

impl Point {
    fn manhattan_distance(&self, other: Point) -> isize {
        (self.x - other.x).abs() + (self.y - other.y).abs()
    }
}

impl Add<Point> for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        Point {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Sub<Point> for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        Point {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Snake {
    head: char,
    points: Vec<Point>,
    is_statue: bool,
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

impl Add<Direction> for Point {
    type Output = Point;

    fn add(self, direction: Direction) -> Point {
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
    snake_head: char,
    direction: Direction,
}

impl Local {
    fn try_move(&mut self, global: &Global, index: usize, direction: Direction) -> bool {
        // Snakes that are actually statues cannot move
        // TODO: Don't generate this move in the first place?
        if self.snakes[index].is_statue {
            return false;
        }

        let head = self.snakes[index].points.first().unwrap().clone();
        let new_head = head + direction;

        // Cannot move into a wall
        if global.tile(new_head) == Tile::Wall {
            return false;
        }

        // Cannot move into a spike
        if global.tile(new_head) == Tile::Spike {
            return false;
        }

        // Cannot move into yourself
        if self.snakes[index].points.contains(&new_head) {
            return false;
        }

        // If we've eaten all the fruit and we're on the exit, we can leave
        // Otherwise, treat exits as empty space
        if self.fruit.is_empty() && global.tile(new_head) == Tile::Exit {
            self.snakes.remove(index);
        } else {
            // Update the snake
            self.snakes[index].points.insert(0, new_head);

            // Potentially eat fruit; don't remove tail if fruit was eaten
            if let Some(fruit_index) = self.fruit.iter().position(|fruit| fruit == &new_head) {
                self.fruit.remove(fruit_index);
            } else {
                self.snakes[index].points.pop();
            }
        }

        // Apply portals
        // If a portal fails, we just don't go through it; this is still a valid move
        // TODO: Push through portals?
        if global.tile(new_head) == Tile::Portal {
            let _ = self.try_portal(global, index, new_head);
        }

        // Push any other snakes
        // If any snake is pushed into a wall, the whole move is invalid
        if !self.try_push(global, index, direction, new_head) {
            return false;
        }

        // Apply gravity to all snakes; if any snakes fall out fo the world, the whole move is invalid
        if !self.try_gravity(global) {
            return false;
        }

        // If we ended up with snakes overlapping, something funny (and bad!) happened
        if self.snakes.iter().enumerate().any(|(i, snake)| {
            self.snakes.iter().enumerate().any(|(j, other_snake)| {
                i != j
                    && snake
                        .points
                        .iter()
                        .any(|point| other_snake.points.contains(point))
            })
        }) {
            return false;
        }

        // If we made it here, the move was successful
        true
    }

    fn try_push(
        &mut self,
        global: &Global,
        index: usize,
        direction: Direction,
        new_head: Point,
    ) -> bool {
        // Attempt to push any snakes in the way
        let mut pushing_indexes = Vec::new();
        let mut pushing_points = vec![new_head];

        'daisies_remain: loop {
            for (other_index, _) in self.snakes.iter().enumerate() {
                // TODO: What if we have a weird loop where we're pushing ourselves?
                if index == other_index {
                    if pushing_points.iter().any(|p| {
                        *p != self.snakes[other_index].points[0]
                            && self.snakes[other_index].points.contains(&p)
                    }) {
                        return false;
                    } else {
                        continue;
                    }
                }

                // Already pushing this snake
                if pushing_indexes.contains(&other_index) {
                    continue;
                }

                // If any pushing point is in the new snake, it's getting pushed too
                if pushing_points
                    .iter()
                    .any(|p| self.snakes[other_index].points.contains(&p))
                {
                    pushing_indexes.push(other_index);

                    pushing_points.extend(
                        self.snakes[other_index]
                            .points
                            .iter()
                            .map(|p| *p + direction),
                    );

                    continue 'daisies_remain;
                }
            }

            break;
        }

        // No snake pushing points can hit anything (other than the original head)
        // Treat exits as empty space when pushing
        // TODO: Is this correct?
        if pushing_points.iter().skip(1).any(|p| {
            !(global.tile(*p) == Tile::Empty
                || global.tile(*p) == Tile::Exit
                || global.tile(*p) == Tile::Portal)
        }) {
            return false;
        }

        // You cannot be pushed into fruit either
        if pushing_points.iter().any(|p| self.fruit.contains(p)) {
            return false;
        }

        // If we have snakes to push, move them all
        for snake_index in pushing_indexes.iter() {
            for point in self.snakes[*snake_index].points.iter_mut() {
                *point = *point + direction;
            }
        }

        // If any snake ends up on a portal, try to portal it
        for snake_index in pushing_indexes.iter() {
            if let Some(portal_index) = self.snakes[*snake_index]
                .points
                .iter()
                .position(|point| global.tile(*point) == Tile::Portal)
            {
                // If going through the portal fails, that's okay! Just don't use it.
                let _ = self.try_portal(
                    global,
                    *snake_index,
                    self.snakes[*snake_index].points[portal_index],
                );
            }
        }

        true
    }

    fn try_gravity(&mut self, global: &Global) -> bool {
        // Apply gravity to all snakes
        'still_falling: loop {
            // Calculate all snakes directly supported
            let mut supported_indexes = Vec::new();

            'finding_support: loop {
                // If all snakes are supported, we're done with both loops
                if supported_indexes.len() == self.snakes.len() {
                    break 'still_falling;
                }

                for (index, _) in self.snakes.iter().enumerate() {
                    // Skip snakes we've already supported (but only inner loop)
                    if supported_indexes.contains(&index) {
                        continue;
                    }

                    // Supported by a wall
                    if self.snakes[index]
                        .points
                        .iter()
                        .any(|point| global.tile(*point + Direction::Down) == Tile::Wall)
                    {
                        supported_indexes.push(index);
                        continue 'finding_support;
                    }

                    // Supported by fruit? (we can walk across fruit)
                    if self.snakes[index]
                        .points
                        .iter()
                        .any(|point| self.fruit.contains(&(*point + Direction::Down)))
                    {
                        supported_indexes.push(index);
                        continue 'finding_support;
                    }

                    // Statues are supported by spikes
                    if self.snakes[index].is_statue {
                        if self.snakes[index]
                            .points
                            .iter()
                            .any(|point| global.tile(*point + Direction::Down) == Tile::Spike)
                        {
                            supported_indexes.push(index);
                            continue 'finding_support;
                        }
                    }

                    // Supported by another supported snake
                    for other_index in supported_indexes.iter() {
                        if self.snakes[index].points.iter().any(|point| {
                            self.snakes[*other_index]
                                .points
                                .contains(&(*point + Direction::Down))
                        }) {
                            supported_indexes.push(index);
                            continue 'finding_support;
                        }
                    }
                }

                // No more supported snakes found
                break 'finding_support;
            }

            // Otherwise, all non supported snakes fall by one
            for index in 0..self.snakes.len() {
                if supported_indexes.contains(&index) {
                    continue;
                }

                // Any portals that we're covering before moving don't double trigger
                // If a portal is covered before falling and after, it doesn't trigger
                let hidden_portals = self.snakes[index]
                    .points
                    .iter()
                    .filter(|point| global.tile(**point) == Tile::Portal)
                    .map(|point| *point)
                    .collect::<Vec<_>>();

                for point in self.snakes[index].points.iter_mut() {
                    *point = *point + Direction::Down;
                }

                // If a snake falls onto a portal, go through it
                // This applies even if it isn't their head
                // TODO: This currently fails on x5 because of this:
                /*
                 * --@--
                 * 210--
                 * #####
                 */
                // If the snake goes straight up through the portal (assume it has to be vertical)
                // It will trigger the fall when the snake goes up three times (sort of jumping)
                // But the game doesn't actually trigger here
                //
                // Also: If a snake was already covering a portal, it can't trigger by falling
                // TODO: Did this fix the above error?
                if let Some(portal_index) = self.snakes[index].points.iter().position(|point| {
                    global.tile(*point) == Tile::Portal && !hidden_portals.contains(point)
                }) {
                    // If going through the portal fails, that's okay! Just don't use it.
                    let _ = self.try_portal(global, index, self.snakes[index].points[portal_index]);
                }
            }

            // If any snakes fall out of the world, the whole move is invalid
            // Statues are allowed to fall out of the world
            let mut statues_to_remove = Vec::new();

            for index in 0..self.snakes.len() {
                if self.snakes[index]
                    .points
                    .iter()
                    .all(|point| point.y > global.height as isize)
                {
                    if self.snakes[index].is_statue {
                        statues_to_remove.push(index);
                    } else {
                        return false;
                    }
                }
            }

            // Remove statues that fell out of the world
            for index in statues_to_remove.iter().rev() {
                self.snakes.remove(*index);
            }

            // If any snakes fell onto spikes, the whole move is invalid
            if self.snakes.iter().any(|snake| {
                snake
                    .points
                    .iter()
                    .any(|point| global.tile(*point) == Tile::Spike)
            }) {
                return false;
            }

            // If any snakes fell into a valid exit, they're gone
            if self.fruit.is_empty() {
                self.snakes
                    .retain(|snake| global.tile(snake.points[0]) != Tile::Exit);
            }

            // But bodies going through the exit are bad (or going through exit before all fruit is gone)
            // if self
            //     .snakes
            //     .iter()
            //     .any(|snake| global.tile(snake.points[0]) == Tile::Exit)
            // {
            //     return false;
            // }
        }

        true
    }

    fn try_portal(&mut self, global: &Global, index: usize, portal_point: Point) -> bool {
        let other_portal = global
            .tiles
            .iter()
            .find_map(|(point, tile)| {
                if *tile == Tile::Portal && *point != portal_point {
                    Some(*point)
                } else {
                    None
                }
            })
            .expect("There must always be a pair of portals");

        // Verify that the other end of the portal is empty
        if self.snakes[index].points.iter().any(|point| {
            let target = *point - portal_point + other_portal;

            // Jumping onto a wall or spike
            if global.tile(target) == Tile::Wall || global.tile(target) == Tile::Spike {
                return true;
            }

            // Jumping onto fruit
            // TODO: Is this fine?
            if self.fruit.contains(&target) {
                return true;
            }

            // Jumping onto another snake
            if self
                .snakes
                .iter()
                .any(|snake| snake.points.contains(&target))
            {
                return true;
            }

            // This point is fine to portal
            false
        }) {
            return false;
        }

        // Move the snake to the other portal
        self.snakes[index].points.iter_mut().for_each(|point| {
            *point = *point - portal_point + other_portal;
        });

        true
    }
}

#[cfg(test)]
mod local_move_tests {
    use super::*;

    fn test_world_1() -> (Global, Local) {
        let input = "\
----
10--
###-";
        Global::read(&mut std::io::Cursor::new(input)).unwrap()
    }

    fn test_world_2() -> (Global, Local) {
        let input = "\
----
10#-
###-";
        Global::read(&mut std::io::Cursor::new(input)).unwrap()
    }

    fn test_big_world_1() -> (Global, Local) {
        let input = "\
------------
------------
------------
3210--------
######---###";
        Global::read(&mut std::io::Cursor::new(input)).unwrap()
    }

    fn test_push_world_1() -> (Global, Local) {
        let input = "\
--a--
10b--
####-";
        Global::read(&mut std::io::Cursor::new(input)).unwrap()
    }

    #[test]
    fn test_move() {
        let (global, mut local) = test_world_1();

        assert!(local.try_move(&global, 0, Direction::Right));
        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 2, y: 1 }, Point { x: 1, y: 1 }]
        );
    }

    #[test]
    fn test_move_up() {
        let (global, mut local) = test_world_1();

        assert!(local.try_move(&global, 0, Direction::Up));
        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 1, y: 0 }, Point { x: 1, y: 1 }]
        );
    }

    #[test]
    fn test_right_up() {
        let (global, mut local) = test_world_1();

        assert!(local.try_move(&global, 0, Direction::Right));
        assert!(local.try_move(&global, 0, Direction::Up));
        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 2, y: 0 }, Point { x: 2, y: 1 }]
        );
    }

    #[test]
    fn test_not_into_ground() {
        let (global, mut local) = test_world_1();

        assert!(!local.try_move(&global, 0, Direction::Down));
    }

    #[test]
    fn test_not_backtrack() {
        let (global, mut local) = test_world_1();

        assert!(!local.try_move(&global, 0, Direction::Left));
    }

    #[test]
    fn test_not_into_wall() {
        let (global, mut local) = test_world_2();

        assert!(!local.try_move(&global, 0, Direction::Right));
    }

    #[test]
    fn test_fall() {
        let (global, mut local) = test_world_1();

        // Move on platform
        assert!(local.try_move(&global, 0, Direction::Right));

        // Hang off platform
        assert!(local.try_move(&global, 0, Direction::Right));

        // Fall off
        assert!(!local.try_move(&global, 0, Direction::Right));
    }

    #[test]
    fn test_move_when_upright_fall() {
        let (global, mut local) = test_world_1();

        // Stand up
        assert!(local.try_move(&global, 0, Direction::Up));

        // Move over
        assert!(local.try_move(&global, 0, Direction::Right));

        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 2, y: 1 }, Point { x: 1, y: 1 }]
        );
    }

    #[test]
    fn test_turn_around() {
        let (global, mut local) = test_world_1();

        // Stand up
        assert!(local.try_move(&global, 0, Direction::Up));

        // Turn and fall back left
        assert!(local.try_move(&global, 0, Direction::Left));

        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 0, y: 1 }, Point { x: 1, y: 1 }]
        );
    }

    #[test]
    fn test_big_curl_up() {
        let (global, mut local) = test_big_world_1();

        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Left));

        assert_eq!(
            local.snakes[0].points,
            vec![
                Point { x: 2, y: 2 },
                Point { x: 3, y: 2 },
                Point { x: 3, y: 3 },
                Point { x: 2, y: 3 },
            ]
        );
    }

    #[test]
    fn test_big_curl_up_too_much() {
        let (global, mut local) = test_big_world_1();

        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Left));
        assert!(!local.try_move(&global, 0, Direction::Down));
    }

    #[test]
    fn test_big_cross_chasm() {
        let (global, mut local) = test_big_world_1();

        for _ in 0..8 {
            assert!(local.try_move(&global, 0, Direction::Right));
        }

        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Left));

        assert_eq!(
            local.snakes[0].points,
            vec![
                Point { x: 10, y: 2 },
                Point { x: 11, y: 2 },
                Point { x: 11, y: 3 },
                Point { x: 10, y: 3 },
            ]
        );
    }

    #[test]
    fn test_big_stand() {
        let (global, mut local) = test_big_world_1();

        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Up));

        assert_eq!(
            local.snakes[0].points,
            vec![
                Point { x: 3, y: 0 },
                Point { x: 3, y: 1 },
                Point { x: 3, y: 2 },
                Point { x: 3, y: 3 }
            ]
        );
    }

    #[test]
    fn test_big_stand_and_fall() {
        let (global, mut local) = test_big_world_1();

        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Up));
        assert!(local.try_move(&global, 0, Direction::Up));

        assert!(local.try_move(&global, 0, Direction::Right));
        assert_eq!(local.snakes[0].points[0], Point { x: 4, y: 1 });

        assert!(local.try_move(&global, 0, Direction::Right));
        assert_eq!(local.snakes[0].points[0], Point { x: 5, y: 2 });

        assert!(local.try_move(&global, 0, Direction::Right));
        assert_eq!(local.snakes[0].points[0], Point { x: 6, y: 3 });

        assert_eq!(
            local.snakes[0].points,
            vec![
                Point { x: 6, y: 3 },
                Point { x: 5, y: 3 },
                Point { x: 4, y: 3 },
                Point { x: 3, y: 3 }
            ]
        );
    }

    #[test]
    fn test_push() {
        let (global, mut local) = test_push_world_1();

        assert!(local.try_move(&global, 0, Direction::Right));
        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 2, y: 1 }, Point { x: 1, y: 1 }]
        );
        assert_eq!(
            local.snakes[1].points,
            vec![Point { x: 3, y: 0 }, Point { x: 3, y: 1 }]
        );
    }

    #[test]
    fn test_push_off() {
        let (global, mut local) = test_push_world_1();

        assert!(local.try_move(&global, 0, Direction::Right));
        assert!(!local.try_move(&global, 0, Direction::Right));
    }

    #[test]
    fn test_push_block() {
        let input = "\
--X-
10X-
####";
        let (global, mut local) = Global::read(&mut std::io::Cursor::new(input)).unwrap();

        assert!(local.try_move(&global, 0, Direction::Right));
        assert_eq!(
            local.snakes[0].points,
            vec![Point { x: 2, y: 1 }, Point { x: 1, y: 1 }]
        );

        assert!(local.snakes[1].points.len() == 2);
        assert!(local.snakes[1].points.contains(&Point { x: 3, y: 0 }));
        assert!(local.snakes[1].points.contains(&Point { x: 3, y: 1 }));
    }
}

impl State<Global, Step> for Local {
    fn is_valid(&self, global: &Global) -> bool {
        self.snakes.iter().all(|snake| {
            snake
                .points
                .iter()
                .all(|point| point.y <= global.height as isize)
        })
    }

    fn is_solved(&self, _map: &Global) -> bool {
        // All fruits eaten and all non-statue snakes are gone
        self.fruit.is_empty() && self.snakes.iter().all(|snake| snake.is_statue)
    }

    fn next_states(&self, global: &Global) -> Option<Vec<(i64, Step, Local)>> {
        let mut next_states = Vec::new();

        // Try to move each snake in each direction
        for (snake_index, _) in self.snakes.iter().enumerate() {
            for direction in [
                Direction::Up,
                Direction::Down,
                Direction::Left,
                Direction::Right,
            ]
            .iter()
            {
                // Generate a potential new local state
                let step = Step {
                    snake_head: self.snakes[snake_index].head,
                    direction: *direction,
                };
                let mut new_local = self.clone();

                // Try to move the snake
                if !new_local.try_move(global, snake_index, *direction) {
                    continue;
                }

                // If the move was successful, add it to the list of next states
                next_states.push((1, step, new_local));
            }
        }

        if next_states.is_empty() {
            return None;
        }
        Some(next_states)
    }

    fn heuristic(&self, global: &Global) -> i64 {
        // Guess how far we have to go yet
        // If there are no fruit, take the straight line from each snake to the exit
        // Otherwise, move each snake to a fruit, then each fruit to each other fruit, then exit
        // TODO: The other case vastly overestimates at this point because of the fruit/fruit distance

        let exit = global
            .tiles
            .iter()
            .find_map(|(point, tile)| {
                if *tile == Tile::Exit {
                    Some(*point)
                } else {
                    None
                }
            })
            .expect("There must always be a single exit");

        if self.fruit.is_empty() {
            // From each snake to the exit
            self.snakes
                .iter()
                .map(|snake| snake.points[0].manhattan_distance(exit))
                .sum::<isize>() as i64
        } else {
            // From each snake to the nearest fruit
            (self.snakes
                .iter()
                .map(|snake| {
                    self.fruit.iter().map(|fruit| snake.points[0].manhattan_distance(*fruit)).min().unwrap_or(0)
                })
                .sum::<isize>()
            // From each fruit to the nearest other fruit
            + self
                .fruit
                .iter()
                .map(|fruit| {
                    self.fruit
                        .iter()
                        .filter(|other_fruit| *other_fruit != fruit)
                        .map(|other_fruit| fruit.manhattan_distance(*other_fruit))
                        .min()
                        .unwrap_or(0)
                })
                .sum::<isize>()
            // From the exit to the nearest fruit
            + self
                .fruit
                .iter()
                .map(|fruit| exit.manhattan_distance(*fruit))
                .min()
                .unwrap_or(0)) as i64
        }
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
                    Tile::Spike => '✴',
                    Tile::Portal => '@',
                };

                for (i, snake) in self.snakes.iter().enumerate() {
                    for (j, pt) in snake.points.iter().enumerate() {
                        if pt.x == x && pt.y == y {
                            // assert!(c == ' ', "Snake in non-empty tile {c}");

                            c = if self.snakes[i].is_statue {
                                self.snakes[i].head
                            } else {
                                (self.snakes[i].head as u8 + j as u8) as char
                            };
                        }
                    }
                }

                if self.fruit.contains(&Point { x, y }) {
                    // assert!(c == ' ', "Fruit in non-empty tile");

                    c = '';
                }

                output.push(c);
            }
            output.push('\n');
        }

        output
    }
}

fn solve(global: Global, local: Local) -> Option<(Solver<Global, Local, Step>, Local)> {
    let mut solver = Solver::new(global.clone(), local);
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
        return None;
    }
    let solution = solution.unwrap();

    log::info!(
        "Solved after {} states in {} seconds:\n{}",
        solver.states_checked(),
        solver.time_spent(),
        solution.stringify(&global),
    );

    Some((solver, solution))
}

fn stringify_solution(
    solver: &Solver<Global, Local, Step>,
    initial_state: &Local,
    solved_state: &Local,
) -> String {
    let mut last_moved_snake = '\0';

    let mut path = String::new();
    for step in solver.path(initial_state, solved_state).unwrap().iter() {
        if step.snake_head != last_moved_snake {
            path.push(step.snake_head);
            last_moved_snake = step.snake_head;
        }

        path.push(match step.direction {
            Direction::Up => '↑',
            Direction::Down => '↓',
            Direction::Left => '←',
            Direction::Right => '→',
        });
    }
    path
}

fn main() -> Result<()> {
    env_logger::init();

    let stdin = std::io::stdin();
    let mut reader = stdin.lock();

    let (global, local) = Global::read(&mut reader).unwrap();

    log::info!("Initial state:\n{}", local.stringify(&global));

    // If there is an arg, assume it's a test case and run it
    if std::env::args().len() > 1 {
        std::env::args().skip(1).for_each(|instructions| {
            let mut local = local.clone();
            let mut current_snake_head = '0';

            for (step, c) in instructions.chars().enumerate() {
                println!("\n=== Step {}: {} ===", step, c);

                if c == '\n' {
                    continue;
                }

                if "0aA{".contains(c) {
                    current_snake_head = c;
                    println!("Switched to snake {c}");
                } else {
                    let direction = match c {
                        '↑' | 'U' | 'u' => Direction::Up,
                        '↓' | 'D' | 'd' => Direction::Down,
                        '←' | 'L' | 'l' => Direction::Left,
                        '→' | 'R' | 'r' => Direction::Right,
                        _ => panic!("invalid instruction: {}", c),
                    };

                    let current_snake_index = local
                        .snakes
                        .iter()
                        .position(|snake| snake.head == current_snake_head)
                        .expect("Snake not found");

                    println!("Moving snake {current_snake_head} ({current_snake_index}) {c} {direction:?}");

                    assert!(local.try_move(&global, current_snake_index, direction));
                    println!("After move:\n{}", local.stringify(&global));
                    println!("Is valid? {}", local.is_valid(&global));
                    println!("Is solved? {}", local.is_solved(&global));
                }
            }
        });
        return Ok(());
    }

    // Otherwise, try to find a new solution
    if let Some((solver, solution)) = solve(global.clone(), local.clone()) {
        let path = stringify_solution(&solver, &local, &solution);
        println!("{path}");

        // Check against known solutions
        if !global.solutions.iter().any(|s| s == &path) {
            log::warn!("Solution does not match known solution")
        }

        return Ok(());
    }

    Err(anyhow!("No solution found"))
}

#[cfg(test)]
mod test_solutions {
    use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
    use std::{fs::File, io::Read, sync::mpsc, thread};

    use super::*;

    #[test]
    fn test_all_solutions() {
        // Timeout after 1 second or SNAKEBIRD_TEST_TIMEOUT if set
        let timeout = std::time::Duration::from_secs(
            std::env::var("SNAKEBIRD_TEST_TIMEOUT")
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(1),
        );

        // Collect all tests to run in order
        let mut test_files = Vec::new();

        let folders = ["data/snakebird", "data/snakebird/primer"];

        for folder in folders.iter() {
            for entry in std::fs::read_dir(folder).unwrap() {
                let entry = entry.unwrap();
                let path = entry.path();
                if !path.is_dir() {
                    continue;
                }

                for entry in std::fs::read_dir(&path).unwrap() {
                    let entry = entry.unwrap();
                    let path = entry.path();

                    if path.extension().is_none() || path.extension().unwrap() != "txt" {
                        continue;
                    }

                    test_files.push(path);
                }
            }
        }
        test_files.sort();
        println!("Running {} tests", test_files.len());

        // Run each test with timeout
        #[derive(Debug, Clone, Eq, PartialEq)]
        enum TestResult {
            Success,
            NoSolution,
            InvalidSolution(String),
            TimedOut,
        }

        let results = test_files
            .par_iter()
            .map(move |path| {
                let mut file = File::open(&path).unwrap();
                let mut input = String::new();
                file.read_to_string(&mut input).unwrap();

                let (global, local) = Global::read(&mut input.as_bytes()).unwrap();

                let (tx, rx) = mpsc::channel();

                let solver_global = global.clone();
                let solver_local = local.clone();

                thread::spawn(move || {
                    let solution = solve(solver_global, solver_local);
                    match tx.send(solution) {
                        Ok(_) => {}
                        Err(_) => {}
                    } // I don't actually care if this succeeds, but need to consume it
                });

                match rx.recv_timeout(timeout) {
                    Ok(solution) => {
                        if solution.is_none() {
                            log::debug!("No solution: {:?}", path);
                            return TestResult::NoSolution;
                        }

                        let (solver, solution) = solution.unwrap();
                        let path = stringify_solution(&solver, &local, &solution);

                        if !global.solutions.contains(&path) {
                            log::debug!("Invalid solution: {:?}", path);
                            return TestResult::InvalidSolution(path);
                        }

                        log::debug!("Solved: {:?}", path);
                        return TestResult::Success;
                    }
                    Err(_) => {
                        log::debug!("Timed out: {:?}", path);
                        return TestResult::TimedOut;
                    }
                }
            })
            .collect::<Vec<_>>();

        // Print out the results
        if results.iter().any(|r| *r == TestResult::TimedOut) {
            println!("\nTimed out tests:");
            for (path, result) in test_files.iter().zip(results.iter()) {
                if *result == TestResult::TimedOut {
                    println!("  {:?}", path);
                }
            }
        }

        if results.iter().any(|r| *r == TestResult::NoSolution) {
            println!("\nUnsolved tests:");
            for (path, result) in test_files.iter().zip(results.iter()) {
                if *result == TestResult::NoSolution {
                    println!("  {:?}", path);
                }
            }
        }

        if results.iter().any(|r| {
            if let TestResult::InvalidSolution(_) = r {
                true
            } else {
                false
            }
        }) {
            println!("\nFailed tests:");
            for (path, result) in test_files.iter().zip(results.iter()) {
                if let TestResult::InvalidSolution(solution) = result {
                    println!("  {:?} -> {:?}", path, solution);
                }
            }
        }

        let perfect = results
            .iter()
            .all(|r| *r == TestResult::Success || *r == TestResult::TimedOut);
        if !perfect {
            println!();
        }
        assert!(perfect);
    }
}
