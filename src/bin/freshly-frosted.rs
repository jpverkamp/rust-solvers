use std::io::Read;

use bitmask_enum::bitmask;
use solver::{Direction, Point, Solver, State};

#[bitmask]
enum Toppings {
    Frosting,
    Sprinkles,
    WhippedCream,
    Cherries,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum EntityKind {
    Block,
    Source,
    Target(Toppings),
    Topper(Toppings),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Entity {
    kind: EntityKind,
    facing: Direction,
}

impl TryFrom<&str> for Entity {
    type Error = String;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        if value.is_empty() {
            return Err("Entity cannot be empty".to_owned());
        }

        fn top(c: char) -> Result<Toppings, String> {
            let v = c.to_digit(16).ok_or(format!("Invalid topping: {c}"))? as usize;
            Ok(v.into())
        }

        fn dir(c: char) -> Result<Direction, String> {
            Ok(Direction::try_from(c).map_err(|_| format!("Invalid direction: {c}"))?)
        }

        match value.chars().collect::<Vec<_>>().as_slice() {
            // A wall (doesn't actually have facing :smile:)
            &['#'] => {
                Ok(Entity {
                    kind: EntityKind::Block,
                    facing: Direction::Up,
                })
            }

            // A source of new donuts
            &['+', facing] => {
                Ok(Entity {
                    kind: EntityKind::Source,
                    facing: dir(facing)?,
                })
            }
            // A target for donuts, first without any toppings
            &['-', facing] => {
                Ok(Entity {
                    kind: EntityKind::Target(Toppings::none()),
                    facing: dir(facing)?,
                })
            }
            &['-', facing, topping] => {
                Ok(Entity {
                    kind: EntityKind::Target(top(topping)?),
                    facing: dir(facing)?,
                })
            }

            // A topping machine
            &[topping, facing] if topping.is_ascii_digit() => {
                Ok(Entity {
                    kind: EntityKind::Topper(top(topping)?),                    
                    facing: dir(facing)?,
                })
            }

            // Something we don't know how to parse
            _ => Err(format!("Invalid entity: {value}")),
        }
        
    }
}

#[derive(Debug, Clone)]
struct Global {
    width: isize,
    height: isize,
    entities: Vec<Option<Entity>>,
}

impl Global {
    fn in_bounds(&self, p: Point) -> bool {
        p.x >= 0 && p.x < self.width && p.y >= 0 && p.y < self.height
    }
}

impl From<&str> for Global {
    fn from(input: &str) -> Self {
        let mut entities = vec![];

        let mut width = 0;
        let mut height = 0;

        for line in input.lines() {
            let mut line_width = 0;

            for part in line.split_whitespace() {
                line_width += 1;

                if part == "." {
                    entities.push(None);
                    continue;
                }

                let entity = match Entity::try_from(part) {
                    Ok(entity) => entity,
                    Err(e) => panic!("Invalid entity: {e}"),
                };

                entities.push(Some(entity));
            }

            width = width.max(line_width);
            height += 1;
        }
        
        Global {
            width,
            height,
            entities,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Local {
    belts: Vec<Option<Direction>>,
}

impl std::fmt::Display for Local {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for belt in self.belts.iter() {
            if let Some(belt) = belt {
                write!(f, "{}", match belt {
                    Direction::Up => "↑",
                    Direction::Down => "↓",
                    Direction::Left => "←",
                    Direction::Right => "→",
                })?;
            } else {
                write!(f, ".")?;
            }
        }
        Ok(())
    }
}

impl Global {
    fn make_local(&self) -> Local {
        Local {
            belts: vec![None; self.width as usize * self.height as usize],
        }
    }
}

impl Local {
    fn simulate(&self, global: &Global) -> Result<Vec<(Point, Toppings)>, String> {
        let mut donuts = vec![];

        // Now we actually have to simulate each source
        for x in 0..global.width {
            for y in 0..global.height {
                let index = (y * global.width + x) as usize;
                let mut p = Point { x, y };
                let mut toppings = Toppings::none();

                // Start a source here
                if let Some(Entity { kind: EntityKind::Source, facing }) = global.entities[index] {
                    let span = tracing::debug_span!("Source", p = ?p);
                    let _enter = span.enter();

                    p = p + facing.into();

                    while !matches!(global.entities[p.index(global.width)], Some(Entity { kind: EntityKind::Target(_), .. })) {
                        tracing::debug!("Moving at {p:?} with toppings {toppings:?}");
                        if !global.in_bounds(p) {
                            return Err(format!("Attempted to move out of bounds at {p:?}"));
                        }

                        // If there is a topper adjacent to us, apply it's topping
                        for top_d in Direction::all() {
                            let p2 = p - top_d.into();
                            if !global.in_bounds(p2) {
                                continue;
                            }

                            if let Some(Entity { kind: EntityKind::Topper(new_toppings), facing }) = global.entities[p2.index(global.width)] {
                                if top_d == facing {
                                    // The toppings cannot overlap what we already have
                                    if toppings & new_toppings != Toppings::none() {
                                        return Err(format!("Attempted to add overlapping toppings {toppings:?} and {new_toppings:?} at {p:?}"));
                                    }

                                    // The toppings have to be added in order
                                    if toppings > new_toppings {
                                        return Err(format!("Attempted to add toppings out of order {toppings:?} and {new_toppings:?} at {p:?}"));
                                    }

                                    // Add the new toppings!
                                    toppings |= new_toppings;
                                }
                            }
                        }

                        // Move move along the belt
                        if let Some(belt) = self.belts[p.index(global.width)] {
                            p = p + belt.into();
                        } else {
                            break; // Ran off the end of a belt
                        }
                    }

                    donuts.push((p, toppings));
                }
            }
        }

        Ok(donuts)
    }
}

impl State<Global, ()> for Local {
    fn is_valid(&self, _global: &Global) -> bool {
        true
    }

    #[tracing::instrument(skip(self, global), fields(belts = %self), ret)]
    fn is_solved(&self, global: &Global) -> bool {
        // Each target must have a belt pointing at it
        for x in 0..global.width {
            for y in 0..global.height {
                let p = Point { x, y };

                if let Some(Entity { kind: EntityKind::Target(_), facing }) = global.entities[p.index(global.width)] {
                    let p2 = p - facing.into();
                    if let Some(facing2) = self.belts[p2.index(global.width)] {
                        if facing != facing2 {
                            tracing::debug!("- At a target but facing the wrong way");
                            return false;
                        }
                    } else {
                        tracing::debug!("- Not at a target");
                        return false;
                    }
                }
            }
        }

        // Simulate the current state
        tracing::debug!("Running simulation");
        let donuts = match self.simulate(global) {
            Ok(donuts) => donuts,
            Err(e) => {
                tracing::debug!("Simulation failed: {e}");
                return false;
            }
        };
        tracing::debug!("Simulation result: {donuts:?}");

        // Each donut must end at a matching target
        for (donuts_p, toppings) in donuts {
            tracing::debug!("Checking donut at {donuts_p:?} with toppings {toppings:?}");
            match global.entities[donuts_p.index(global.width)] {
                Some(Entity { kind: EntityKind::Target(target_toppings), .. }) if toppings == target_toppings => {},
                _ => return false,
            }
        }
        
        true
    }

    #[tracing::instrument(skip(self, global), fields(belts = %self))]
    fn next_states(&self, global: &Global) -> Option<Vec<(i64, (), Local)>> {
        let mut next_states = vec![];

        // For each 'head', add a belt in each valid direction
        // A head is a source or end of a belt pointing at nothing
        // A valid direction points at empty space or a goal
        for x in 0..global.width {
            for y in 0..global.height {
                let p = Point { x, y };

                let span = tracing::debug_span!("Point", p = ?p);
                let _enter = span.enter();

                // If this is not a head, skip it
                let is_belt = self.belts[p.index(global.width)].is_some();
                
                if !(is_belt || matches!(global.entities[p.index(global.width)], Some(Entity { kind: EntityKind::Source, .. }))) {
                    tracing::debug!("Skipping, this is not a head");
                    continue;
                }

                // p2 is the new belt, so this point must be currently empty
                let facing = if is_belt { 
                    self.belts[p.index(global.width)].unwrap()
                } else {
                    global.entities[p.index(global.width)].unwrap().facing 
                };
                let p2 = p + facing.into();

                if self.belts[p2.index(global.width)].is_some() || global.entities[p2.index(global.width)].is_some() {
                    tracing::debug!("Skipping, contains a belt or entity");
                    continue;
                }

                // Now, for each direction from *that* point, we can potentially add a belt
                for d2 in Direction::all() { 
                    let p3 = p2 + d2.into();

                    let span = tracing::debug_span!("Checking direction", d = ?d2);
                    let _enter = span.enter();

                    // The new point must be in bounds
                    if !global.in_bounds(p3) {
                        tracing::debug!("Skipping, out of bounds");
                        continue;
                    }

                    // This belt can only point at non-belt, non-target
                    if self.belts[p3.index(global.width)].is_some() {
                        tracing::debug!("Skipping, contains a belt");
                        continue;
                    }
                    match global.entities[p3.index(global.width)] {
                        None => {},
                        Some(Entity { kind: EntityKind::Target(_), facing }) if facing == d2 => {},
                        _ => {
                            tracing::debug!("Skipping, contains a non-target");
                            continue;
                        }
                    }

                    // If we made it this far, this is a valid new state
                    let mut new_state = self.clone();
                    new_state.belts[p2.index(global.width)] = Some(d2);
                    tracing::debug!("Valid new state: {}", &new_state);

                    next_states.push((1, (), new_state));                    
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
        let mut chars = vec![vec![' '; global.width as usize]; global.height as usize];

        for y in 0..global.height {
            for x in 0..global.width {
                let index = (y * global.width + x) as usize;
                chars[y as usize][x as usize] = '.';

                // An entity from the global
                let entity = &global.entities[index];
                if let Some(entity) = entity {
                    chars[y as usize][x as usize] = match entity.kind {
                        EntityKind::Block => '#',
                        EntityKind::Source => '+',
                        EntityKind::Target(_) => '-',
                        EntityKind::Topper(_) => {
                            match entity.facing {
                                Direction::Up => '╩',
                                Direction::Down => '╦',
                                Direction::Left => '╣',
                                Direction::Right => '╠',
                            }
                        }
                    };
                }

                // A belt from the local
                let belt = &self.belts[index];
                if let Some(belt) = belt {
                    chars[y as usize][x as usize] = match belt {
                        Direction::Up => '↑',
                        Direction::Down => '↓',
                        Direction::Left => '←',
                        Direction::Right => '→',
                    };
                }

                // It's an error if not at least one of these is None
                assert!(entity.is_none() || belt.is_none());
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

fn main() {
    let mut input = String::new();
    std::io::stdin().read_to_string(&mut input).unwrap();
    let global = Global::from(input.as_str());
    let local = global.make_local();

    // tracing_subscriber::fmt().without_time().with_max_level(tracing::Level::DEBUG).init();
    env_logger::init();

    log::info!("Initial state:\n{}", local.stringify(&global));

    let mut solver = Solver::new(global.clone(), local.clone());

    while let Some(state) = solver.next() {
        if solver.states_checked() % 100_000 == 0 {
            log::debug!("\n{}", state.stringify(&global));
            log::debug!("{solver}");
        }
    }
    let solution = solver.get_solution();

    if let Some(solution) = solution {
        println!("{}", solver.stringify(&solution));

        let _path = solver.path(&local, &solution).unwrap();

    } else {
        println!("No solution found");
        std::process::exit(1);
    }
}
