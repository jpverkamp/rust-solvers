use std::{
    collections::HashMap,
    io::Read,
    sync::{LazyLock, Mutex},
};

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
    Target(Option<Toppings>),
    Topper(Toppings),
    Splitter,
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
            &['#'] => Ok(Entity {
                kind: EntityKind::Block,
                facing: Direction::Up,
            }),

            // A source of new donuts
            &['+', facing] => Ok(Entity {
                kind: EntityKind::Source,
                facing: dir(facing)?,
            }),

            // A target for donuts, first without any toppings
            &['-', facing] => Ok(Entity {
                kind: EntityKind::Target(Some(Toppings::none())),
                facing: dir(facing)?,
            }),
            // Doesn't matter what toppings
            &['-', facing, '?'] => Ok(Entity {
                kind: EntityKind::Target(None),
                facing: dir(facing)?,
            }),
            // Now with specific requested toppings
            &['-', facing, topping] => Ok(Entity {
                kind: EntityKind::Target(Some(top(topping)?)),
                facing: dir(facing)?,
            }),

            // A topping machine
            &[topping, facing] if topping.is_ascii_digit() => Ok(Entity {
                kind: EntityKind::Topper(top(topping)?),
                facing: dir(facing)?,
            }),

            // A splitter
            &['X', facing] | &['x', facing] => Ok(Entity {
                kind: EntityKind::Splitter,
                facing: dir(facing)?,
            }),

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
    targets: Option<Vec<Option<Toppings>>>,
    initial_belts: Vec<Option<Direction>>,
}

impl Global {
    fn in_bounds(&self, p: Point) -> bool {
        p.x >= 0 && p.x < self.width && p.y >= 0 && p.y < self.height
    }
}

impl From<&str> for Global {
    fn from(input: &str) -> Self {
        let mut width = 0;
        let mut height = 0;

        let mut entities = vec![];
        let mut initial_belts = vec![];
        let mut targets = None;

        let mut lines = input.lines().peekable();
        while lines.peek().is_some_and(|line| line.starts_with(':')) {
            let flag = lines.next().unwrap();

            if flag.starts_with(":target") {
                targets = Some(
                    flag.split_whitespace()
                        .skip(1)
                        .map(|t| {
                            Some(
                                t.parse::<usize>()
                                    .expect("Invalid target, must be numeric")
                                    .into(),
                            )
                        })
                        .collect::<Vec<_>>(),
                );
            } else {
                panic!("Invalid/unknown flag: {flag}");
            }
        }

        for line in lines {
            let mut line_width = 0;

            for part in line.split_whitespace() {
                let mut belt = None;
                let mut entity = None;
                line_width += 1;

                if part == "." {
                    // Do nothing
                } else if let Ok(found_belt) = Direction::try_from(part) {
                    belt = Some(found_belt);
                } else {
                    entity = match Entity::try_from(part) {
                        Ok(entity) => Some(entity),
                        Err(e) => panic!("Invalid entity: {e}"),
                    };
                }

                assert!(!(entity.is_some() && belt.is_some()));

                initial_belts.push(belt);
                entities.push(entity);
            }

            width = width.max(line_width);
            height += 1;
        }

        Global {
            width,
            height,
            entities,
            initial_belts,
            targets,
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
                write!(
                    f,
                    "{}",
                    match belt {
                        Direction::Up => "↑",
                        Direction::Down => "↓",
                        Direction::Left => "←",
                        Direction::Right => "→",
                    }
                )?;
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
            belts: self.initial_belts.clone(),
        }
    }
}

static SIMULATE_CACHE: LazyLock<Mutex<HashMap<Local, Vec<(Point, Toppings)>>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

impl Local {
    fn simulate(&self, global: &Global) -> Result<Vec<(Point, Toppings)>, String> {
        // Check the cache first
        if let Some(cached) = SIMULATE_CACHE.lock().unwrap().get(self) {
            return Ok(cached.clone());
        }

        let mut donuts = vec![];
        let mut complete_donuts = vec![];

        // Now we actually have to simulate each source

        // First, find the sources:
        for x in 0..global.width {
            for y in 0..global.height {
                let index = (y * global.width + x) as usize;
                let p = Point { x, y };
                let toppings = Toppings::none();

                // Start a source here
                if let Some(Entity {
                    kind: EntityKind::Source,
                    facing,
                }) = global.entities[index]
                {
                    donuts.push((p, toppings, facing));
                }
            }
        }

        // Now, until we've simulated all of the donuts, simulate each
        while let Some((mut p, mut toppings, facing)) = donuts.pop() {
            let span = tracing::debug_span!("Simulating donut", p = ?p, toppings = ?toppings, facing = ?facing);
            let _enter = span.enter();
            tracing::debug!("Starting");

            let mut visited = vec![false; global.width as usize * global.height as usize];

            p = p + facing.into();
            visited[p.index(global.width)] = true;

            while !matches!(
                global.entities[p.index(global.width)],
                Some(Entity {
                    kind: EntityKind::Target(_),
                    ..
                })
            ) {
                if !global.in_bounds(p) {
                    return Err(format!("Attempted to move out of bounds at {p:?}"));
                }

                // If there is a topper adjacent to us, apply it's topping
                for top_d in Direction::all() {
                    let p2 = p - top_d.into();
                    if !global.in_bounds(p2) {
                        continue;
                    }

                    if let Some(Entity {
                        kind: EntityKind::Topper(new_toppings),
                        facing,
                    }) = global.entities[p2.index(global.width)]
                    {
                        if top_d == facing {
                            // Only add the toppings if we have all previous stoppings, otherwise ignore it
                            // I feel like this was frowned upon before 4/6, but so it goess bits set to add a new one
                            assert!(new_toppings.bits().count_ones() == 1);
                            let must_have = Toppings::from(new_toppings.bits() - 1);
                            if toppings & must_have != must_have {
                                continue;
                            }

                            // Add the new topping!
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

                // If we are now in a splitter, continue one way and queue the other
                if let Some(Entity {
                    kind: EntityKind::Splitter,
                    facing: splitter_facing,
                }) = global.entities[p.index(global.width)]
                {
                    // If we're coming into a splitter the wrong way
                    if facing != splitter_facing {
                        return Err(format!("Attempted to enter the splitter at {p:?} the wrong way"));
                    }
                    
                    // This is the donut we'll queue up to do later
                    tracing::debug!("Splitting donut at {p:?} with toppings {toppings:?}, queuing {:?}", facing.turn_left());
                    donuts.push((p, toppings, facing.turn_left()));

                    // This is the donut that we're continuing with now
                    // TODO: We should probably validity check this?
                    tracing::debug!("Splitting donut at {p:?} with toppings {toppings:?}, continuing {:?}", facing.turn_right());
                    let facing = facing.turn_right();
                    p = p + facing.into();
               }

                // Error on loops
                if visited[p.index(global.width)] {
                    return Err(format!("Loop detected at {p:?}"));
                } else {
                    visited[p.index(global.width)] = true;
                }
            }

            complete_donuts.push((p, toppings));
        }

        // Cache the result
        SIMULATE_CACHE
            .lock()
            .unwrap()
            .insert(self.clone(), complete_donuts.clone());

        Ok(complete_donuts)
    }

    // Is it at all possible to get from src to dst with the current belt configuration?
    // Use empty points, allowed to step on targets, and can follow belts
    #[tracing::instrument(skip(self, global), ret)]
    fn is_reachable(&self, global: &Global, src: Point, dst: Point) -> bool {
        if src == dst {
            return true;
        }

        pathfinding::prelude::bfs(
            &src,
            |p| {
                let mut neighbors = vec![];
                for d in Direction::all() {
                    let p2 = *p + d.into();

                    if !global.in_bounds(p2) {
                        continue;
                    }

                    let is_belt = self.belts[p2.index(global.width)].is_some();
                    let is_belt_in_proper_direction = self.belts[p2.index(global.width)]
                        .map(|belt_d| belt_d == d)
                        .unwrap_or(false);

                    let is_target = matches!(
                        global.entities[p2.index(global.width)],
                        Some(Entity {
                            kind: EntityKind::Target(_),
                            ..
                        })
                    );

                    if !is_belt || is_belt_in_proper_direction || (p2 == dst && is_target) {
                        neighbors.push(p2);
                    }
                }
                neighbors
            },
            |p| *p == dst,
        )
        .is_some()
    }
}

impl State<Global, ()> for Local {
    fn is_valid(&self, global: &Global) -> bool {
        // This is expensive to run on each state, but it also means that we can prune a *lot* of invalid states
        let donuts = match self.simulate(global) {
            Ok(donuts) => donuts,
            Err(e) => {
                tracing::debug!("Simulation failed: {e}");
                return false;
            }
        };
        tracing::debug!("Simulation result: {donuts:?}");

        let target_points = global
            .entities
            .iter()
            .enumerate()
            .filter_map(|(index, entity)| {
                if let Some(Entity {
                    kind: EntityKind::Target(_),
                    ..
                }) = entity
                {
                    Some(Point {
                        x: index as isize % global.width,
                        y: index as isize / global.width,
                    })
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        // Each current donut must be either on a target or be able to reach one
        tracing::debug!("Checking reachability");
        for (donut_p, toppings) in donuts.iter() {
            tracing::debug!("Checking donut at {donut_p:?} with toppings {toppings:?}");
            if !target_points
                .iter()
                .any(|target_p| self.is_reachable(global, *donut_p, *target_p))
            {
                return false;
            }
        }

        // We cannot have an over filled simulation
        // For example, if we need 1 plain donut, we cannot have frosting on all the donuts
        let mut donut_types = donuts
            .into_iter()
            .map(|(_, toppings)| Some(toppings))
            .collect::<Vec<_>>();
        donut_types.sort();

        let mut target_donuts = global
            .entities
            .iter()
            .filter_map(|e| match e {
                Some(Entity {
                    kind: EntityKind::Target(toppings),
                    ..
                }) => Some(*toppings),
                _ => None,
            })
            .collect::<Vec<_>>();
        target_donuts.sort();

        // We can't do this comparison if we have any 'any' targets
        if target_donuts.iter().all(|t| t.is_some()) {
            // This comparison also doesn't work if we have fewer targets than dummies (world 4: merge)
            if donut_types.len() == target_donuts.len() {
                // Because they're sorted, this means have at least one donut with too many toppings
                if donut_types > target_donuts {
                    tracing::debug!(
                        "Not possible. Got {donut_types:?} but needed {target_donuts:?}"
                    );
                    return false;
                }
            }
        }

        // But if we have a target in the globals, we must (also?) match that
        if let Some(target_donuts) = &global.targets {
            if &donut_types > target_donuts {
                tracing::debug!("Not possible against global types. Got {donut_types:?} but needed {target_donuts:?}");
                return false;
            }
        }

        true
    }

    #[tracing::instrument(skip(self, global), fields(belts = %self), ret)]
    fn is_solved(&self, global: &Global) -> bool {
        // Each target must have a belt pointing at it
        // TODO: Can a splitter directly point at an exit? 
        for x in 0..global.width {
            for y in 0..global.height {
                let p = Point { x, y };

                if let Some(Entity {
                    kind: EntityKind::Target(_),
                    facing,
                }) = global.entities[p.index(global.width)]
                {
                    let p2 = p - facing.into();
                    if let Some(facing2) = self.belts[p2.index(global.width)] {
                        if facing != facing2 {
                            tracing::debug!("At a target but facing the wrong way");
                            return false;
                        }
                    } else {
                        tracing::debug!("Not at a target");
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
        // A non-target accepts any donut
        for (donuts_p, toppings) in donuts.iter() {
            tracing::debug!("Checking donut at {donuts_p:?} with toppings {toppings:?}");
            match global.entities[donuts_p.index(global.width)] {
                Some(Entity {
                    kind: EntityKind::Target(target_toppings),
                    ..
                }) if target_toppings.is_none() || *toppings == target_toppings.unwrap() => {}
                _ => {
                    tracing::debug!("Invalid solution, donut at {donuts_p:?} with toppings {toppings:?} isn't at a matching target");
                    return false;
                }
            }
        }

        // If we have a toppings array, match that as well/instead
        // The check above already checks that they're all at targets, this just checks types
        if let Some(target_donuts) = &global.targets {
            let mut donut_types = donuts.iter().map(|(_, t)| Some(*t)).collect::<Vec<_>>();
            donut_types.sort();

            if &donut_types != target_donuts {
                tracing::debug!("Invalid solution, types don't match. Got {donut_types:?} but needed {target_donuts:?}");
                return false;
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
                // A head is a belt, source, or splitter
                let is_belt = self.belts[p.index(global.width)].is_some();
                let is_source = matches!(
                    global.entities[p.index(global.width)],
                    Some(Entity {
                        kind: EntityKind::Source,
                        ..
                    })
                );
                let is_splitter = matches!(
                    global.entities[p.index(global.width)],
                    Some(Entity {
                        kind: EntityKind::Splitter,
                        ..
                    })
                );

                if !(is_belt || is_source || is_splitter) {
                    continue;
                }

                // p2 is the new belt, so this point must be currently empty
                let facing = if is_belt {
                    self.belts[p.index(global.width)].unwrap()
                } else if is_source {
                    global.entities[p.index(global.width)].unwrap().facing
                } else if is_splitter {
                    let splitter_facing = global.entities[p.index(global.width)].unwrap().facing;
                    
                    // Try one way here and the other below, this is hacky to add on, but so it goes
                    let p2 = p + splitter_facing.turn_left().into();

                    // If we would continue with turn left, return turn right as facing and check it below
                    // If we wouldn't, return turn left and we'll pass that check below too
                    if self.belts[p2.index(global.width)].is_some() || global.entities[p2.index(global.width)].is_some() {
                        splitter_facing.turn_right()
                    } else {
                        splitter_facing.turn_left()
                    }

                } else {
                    unreachable!()
                };
                let p2 = p + facing.into();

                if self.belts[p2.index(global.width)].is_some()
                    || global.entities[p2.index(global.width)].is_some()
                {
                    continue;
                }

                // Now, for each direction from *that* point, we can potentially add a belt
                for d2 in Direction::all() {
                    // for d2 in [Direction::Right] {
                    // We cannot go back the way we came
                    // This wasn't a problem until world 4 allowed merging
                    if d2 == facing.flip() {
                        continue;
                    }

                    let p3 = p2 + d2.into();

                    let span = tracing::debug_span!("Checking direction", d = ?d2);
                    let _enter = span.enter();

                    // The new point must be in bounds
                    if !global.in_bounds(p3) {
                        tracing::debug!("Skipping, out of bounds");
                        continue;
                    }

                    // This belt can point to:
                    // - Empty space
                    // - Other belts (as of world 4: merging)
                    // - Targets
                    // - Splitters (as of world 6)
                    match global.entities[p3.index(global.width)] {
                        None => {}
                        Some(Entity {
                            kind: EntityKind::Target(_) | EntityKind::Splitter,
                            facing,
                        }) if facing == d2 => {}
                        _ => {
                            tracing::debug!("Skipping, directed to a non-valid space/entity");
                            continue;
                        }
                    }

                    // If we made it this far, this is a valid new state
                    let mut new_state = self.clone();
                    new_state.belts[p2.index(global.width)] = Some(d2);
                    tracing::debug!("Valid new state: {}", &new_state);

                    next_states.push((1, (), new_state));
                }

                // DEBUG: If we are here, we've expanded one head in x/y order
                if !next_states.is_empty() {
                    return Some(next_states);
                }
            }
        }

        if next_states.is_empty() {
            return None;
        }
        Some(next_states)
    }

    fn heuristic(&self, global: &Global) -> i64 {
        // For each donut, the distance to the nearest reachable target
        let donuts = match self.simulate(global) {
            Ok(donuts) => donuts,
            Err(_) => return 0,
        };

        let target_points = global
            .entities
            .iter()
            .enumerate()
            .filter_map(|(index, entity)| {
                if let Some(Entity {
                    kind: EntityKind::Target(_),
                    ..
                }) = entity
                {
                    Some(Point {
                        x: index as isize % global.width,
                        y: index as isize / global.width,
                    })
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        donuts
            .iter()
            .map(|(donut_p, _)| {
                target_points
                    .iter()
                    .filter(|target_p| self.is_reachable(global, *donut_p, **target_p))
                    .map(|target_p| donut_p.manhattan_distance(*target_p))
                    .min()
                    .unwrap_or(0) as i64
            })
            .sum()
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
                        EntityKind::Topper(_) => match entity.facing {
                            Direction::Up => '╩',
                            Direction::Down => '╦',
                            Direction::Left => '╣',
                            Direction::Right => '╠',
                        },
                        EntityKind::Splitter => 'X',
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

    let tracing_enabled = std::env::var("FRESHLY_FROSTED_TRACE").is_ok();

    if tracing_enabled {
        tracing_subscriber::fmt()
            .without_time()
            .with_max_level(tracing::Level::DEBUG)
            .init();
    } else {
        env_logger::init();
    }

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
