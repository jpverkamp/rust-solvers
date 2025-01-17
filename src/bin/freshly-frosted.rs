use std::{
    collections::{HashMap, HashSet},
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
    Bumper(Toppings),
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

            // A bumper
            &['b', facing, topping] => Ok(Entity {
                kind: EntityKind::Bumper(top(topping)?),
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
    use_tickwise: bool,
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
        let mut use_tickwise = false;

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
            } else if flag.starts_with(":comment") {
                // Not stored, just do nothing   
            } else if flag.starts_with(":tickwise") {
                use_tickwise = true;
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
            use_tickwise,
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

static SIMULATE_TICKWISE_CACHE: LazyLock<Mutex<HashMap<Local, SimulateTickwiseResult>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

#[allow(dead_code)]
#[derive(Debug, Clone)]
struct SimulateTickwiseResult {
    deliveries: HashMap<Point, HashSet<(Point, Toppings)>>,
}

impl Local {
    fn is_empty(&self, global: &Global, p: Point) -> bool {
        if !global.in_bounds(p) {
            return false;
        }

        self.belts[p.index(global.width)].is_none()
            && global.entities[p.index(global.width)].is_none()
    }

    fn simulate_tickwise(&self, global: &Global) -> Result<SimulateTickwiseResult, String> {
        // Check the cache first
        if let Some(cached) = SIMULATE_TICKWISE_CACHE.lock().unwrap().get(self) {
            return Ok(cached.clone());
        }

        let vec_size = global.width as usize * global.height as usize;

        #[derive(Debug, Clone, Default, PartialEq, Eq, Hash)]
        struct TileState {
            toppings: Option<Toppings>,
            source: Option<Point>,
            updated: bool,
            waiting_time: usize,
            split_next_right: bool,
        }

        let mut state = vec![TileState::default(); vec_size];
        let mut deliveries = HashMap::new();
        let mut states_seen = HashSet::new();

        macro_rules! state_at {
            ($p:expr) => {
                state[$p.index(global.width)]
            };
        }

        // We have a very expensive tracing option; so don't calculate it if we're not going to print it
        let tracing_enabled = std::env::var("FRESHLY_FROSTED_TRACE").is_ok();

        // Advance the simulation one tick
        'tick: loop {
            state.iter_mut().for_each(|s| s.updated = false);

            // Debugging ticking
            if tracing_enabled {
                let mut map = self.stringify(global).chars().collect::<Vec<_>>();

                for y in 0..global.height {
                    for x in 0..global.width {
                        let p = Point { x, y };
                        if let Some(toppings) = state_at!(p).toppings {
                            map[p.index(global.width + 1)] =
                                toppings.bits.to_string().chars().next().unwrap();
                        }
                    }
                }

                let map = map.iter().collect::<String>();
                let cache_size = states_seen.len();
                let max_waiting_time = state.iter().map(|s| s.waiting_time).max().unwrap_or(0);

                tracing::debug!(
                    "\
=== Starting tick ===
Deliveries: {deliveries:?}
States seen: {cache_size}
Max waiting time: {max_waiting_time}

{map}
",
                );
            }

            // Cache which exact states we've seen; break once we see the same more than once
            // TODO: This is expensive..

            // TODO: I don't think we should be able to get away with just caching the toppings, but it's working so far
            if !states_seen.insert(state.clone()) {
                tracing::debug!("Loop detected, breaking");
                tracing::debug!("Deliveries are: {deliveries:#?}");
                break 'tick;
            }

            // Look for a space that can be updated
            // This is always a destination space, find a donut that can be put there
            'find_update: loop {
                for x in 0..global.width {
                    for y in 0..global.height {
                        let p = Point { x, y };

                        // Only update each point once
                        if state_at!(p).updated {
                            continue;
                        }

                        // If the spot is occupied, skip it
                        if state_at!(p).toppings.is_some() {
                            continue;
                        }

                        // Find the potential donuts that can move into this space
                        let mut potential_donuts = vec![];
                        for d in Direction::all() {
                            let p_from = p - d.into();
                            if !global.in_bounds(p_from) {
                                continue;
                            }

                            // TODO: We can only move into a target or splitter the correct way

                            // A belt pointing at us, containing a donut, and didn't just get that donut
                            if let Some(belt) = self.belts[p_from.index(global.width)] {
                                if belt == d && !state_at!(p_from).updated {
                                    if let Some(toppings) = state_at!(p_from).toppings {
                                        potential_donuts.push((p_from, toppings, d));
                                    }
                                }
                                continue;
                            }

                            // A source pointing at us
                            if let Some(Entity {
                                kind: EntityKind::Source,
                                facing,
                            }) = global.entities[p_from.index(global.width)]
                            {
                                if facing == d {
                                    potential_donuts.push((p_from, Toppings::none(), d));
                                }
                                continue;
                            }

                            // A splitter pointing at us + on the correct cycle
                            if let Some(Entity {
                                kind: EntityKind::Splitter,
                                facing,
                            }) = global.entities[p_from.index(global.width)]
                            {
                                let next_left = !state_at!(p_from).split_next_right;

                                if next_left && facing.turn_left() == d
                                    || !next_left && facing.turn_right() == d
                                {
                                    if let Some(toppings) = state_at!(p_from).toppings {
                                        potential_donuts.push((p_from, toppings, d));
                                    }
                                }

                                continue;
                            }

                            // A bumper two spaces away matching a donut on a belt one tile away
                            let p_bumper = p_from - d.into();
                            if global.in_bounds(p_bumper) {
                                if let Some(Entity {
                                    kind: EntityKind::Bumper(bumper_toppings),
                                    facing,
                                }) = global.entities[p_bumper.index(global.width)]
                                {
                                    if facing == d {
                                        if let Some(toppings) = state_at!(p_from).toppings {
                                            if toppings == bumper_toppings {
                                                potential_donuts.push((p_from, toppings, d));
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        // If we have no potential donuts, try the next position
                        // TODO: Is this correct to treat it as updated?
                        if potential_donuts.is_empty() {
                            state_at!(p).updated = true;
                            continue;
                        }

                        tracing::debug!(
                            "Found potential donuts at {p:?}: {potential_donuts:?}",
                            p = p,
                            potential_donuts = potential_donuts
                        );

                        // If we have at least one donut, make sure we're moving it onto something valid
                        if let Some(entity) = global.entities[p.index(global.width)] {
                            match entity.kind {
                                // Can never be moved onto
                                EntityKind::Block
                                | EntityKind::Source
                                | EntityKind::Topper(_)
                                | EntityKind::Bumper(_) => {
                                    return Err(format!(
                                        "Attempted to move donut onto a {:?} at {:?}",
                                        entity.kind, p
                                    ))
                                }
                                // Can only be moved onto matching the facing
                                EntityKind::Target(_) | EntityKind::Splitter => {
                                    if entity.facing != potential_donuts[0].2 {
                                        return Err(format!("Attempted to move donut onto a {:?} at {:?} facing the wrong way (expected {:?}, got {:?})", entity.kind, p, entity.facing, potential_donuts[0].2));
                                    }
                                }
                            };
                        }

                        // If we have multiple potential donuts, choose the one waiting longest
                        // TODO: What if there is a tie?
                        potential_donuts.sort_by(|(p1, _, _), (p2, _, _)| {
                            state_at!(p2).waiting_time.cmp(&state_at!(p1).waiting_time)
                        });

                        // Potential donut 0 moves here and resets wait time
                        state_at!(p).toppings = Some(potential_donuts[0].1);
                        state_at!(potential_donuts[0].0).toppings = None;
                        state_at!(potential_donuts[0].0).waiting_time = 0;

                        // Update the source as well, if the previous had a source, that's what we get, otherwise it's the source point
                        if let Some(source) = state_at!(potential_donuts[0].0).source {
                            state_at!(p).source = Some(source);
                        } else {
                            state_at!(p).source = Some(potential_donuts[0].0);
                        }

                        // If potential 0 was a splitter, toggle it
                        if let Some(Entity {
                            kind: EntityKind::Splitter,
                            ..
                        }) = global.entities[potential_donuts[0].0.index(global.width)]
                        {
                            state_at!(potential_donuts[0].0).split_next_right =
                                !state_at!(potential_donuts[0].0).split_next_right;
                        }

                        // Any others increment their wait time
                        potential_donuts.iter().skip(1).for_each(|(p_from, _, _)| {
                            state_at!(p_from).waiting_time += 1;
                        });

                        // If there is a valid topper pointing at us, apply it
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
                                    assert!(new_toppings.bits().count_ones() == 1);
                                    let must_have = Toppings::from(new_toppings.bits() - 1);
                                    if state_at!(p).toppings.unwrap() & must_have != must_have {
                                        continue;
                                    }

                                    // Add the new topping!
                                    state_at!(p).toppings =
                                        Some(state_at!(p).toppings.unwrap() | new_toppings);
                                }
                            }
                        }

                        // If we are now sitting on a target, either deliver the donut or break
                        if let Some(Entity {
                            kind: EntityKind::Target(target_toppings),
                            ..
                        }) = global.entities[p.index(global.width)]
                        {
                            // Always remove if there's no target
                            let toppings = state_at!(p).toppings.unwrap();
                            let valid = match target_toppings {
                                None => true,
                                Some(target_toppings) => toppings == target_toppings,
                            };

                            if valid {
                                state_at!(p).toppings = None;
                                deliveries
                                    .entry(p)
                                    .or_insert_with(HashSet::new)
                                    .insert((state_at!(p).source.unwrap(), toppings));
                            } else {
                                return Err(format!("Attempted to deliver {toppings:?} to {p:?} but needed {target_toppings:?}"));
                            }
                        }

                        // If we made it this far, there was an update to this space, log it and continue
                        state_at!(p).updated = true;
                        continue 'find_update;
                    }
                }

                // If we make it here without continuing, there are no more updates
                break 'find_update;
            }
        }

        let result = SimulateTickwiseResult { deliveries };

        // Cache the result
        SIMULATE_TICKWISE_CACHE
            .lock()
            .unwrap()
            .insert(self.clone(), result.clone());

        Ok(result)
    }

    fn simulate(&self, global: &Global) -> Result<Vec<(Point, Toppings)>, String> {
        // Check the cache first
        if let Some(cached) = SIMULATE_CACHE.lock().unwrap().get(self) {
            return Ok(cached.clone());
        }

        let mut donuts = vec![];
        let mut complete_donuts = vec![];

        // Now we actually have to simulate each source
        let mut visited = vec![];
        for _ in 0..(Toppings::all_flags().bits + 1) {
            visited.push(vec![false; global.width as usize * global.height as usize]);
        }

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
                    donuts.push((p, toppings, facing, visited.clone()));
                }
            }
        }

        // Now, until we've simulated all of the donuts, simulate each
        'each_donut: while let Some((mut p, mut toppings, initial_facing, mut visited)) =
            donuts.pop()
        {
            let span = tracing::debug_span!("Simulating donut", p = ?p, t = ?toppings, f = ?initial_facing);
            let _enter = span.enter();

            p = p + initial_facing.into();
            if !global.in_bounds(p) {
                return Err(format!(
                    "Donut started at {p:?} but immediately went out of bounds"
                ));
            }
            
            visited[toppings.bits][p.index(global.width)] = true;

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

                // tracing::debug!("Step {p:?}, visited: {}",
                //     visited
                //         .iter()
                //         .enumerate()
                //         .flat_map(|(i, v)|
                //             if *v {
                //                 Some(format!("({x}, {y})",
                //                     x = i as isize % global.width,
                //                     y = i as isize / global.width,
                //                 ))
                //             } else {
                //                 None
                //             }
                //         )
                //         .collect::<Vec<_>>()
                //         .join(" ")
                //     );

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

                // If we are at a splitter, split the donut
                if let Some(Entity {
                    kind: EntityKind::Splitter,
                    facing,
                }) = global.entities[p.index(global.width)]
                {
                    // TODO: Check if we came into the splitter the wrong way?

                    // Queue the two new donuts
                    donuts.push((p, toppings, facing.turn_left(), visited.clone()));
                    donuts.push((p, toppings, facing.turn_right(), visited.clone()));

                    // Do not simulate this path any more
                    continue 'each_donut;
                }

                // If we're adjacent to a matching bumper, bump
                // TODO: Assume that the space we move onto is valid
                for bump_d in Direction::all() {
                    let p2 = p - bump_d.into();
                    if !global.in_bounds(p2) {
                        continue;
                    }

                    if let Some(Entity {
                        kind: EntityKind::Bumper(bump_toppings),
                        facing,
                    }) = global.entities[p2.index(global.width)]
                    {
                        tracing::debug!("Bumper at {p2:?} with facing {bump_d:?} and toppings {bump_toppings:?}");
                        tracing::debug!("Donut at {p:?} with toppings {toppings:?}");

                        if bump_d == facing && toppings == bump_toppings {
                            p = p + bump_d.into();
                            break;
                        }
                    }
                }

                // Move move along the belt
                if let Some(belt) = self.belts[p.index(global.width)] {
                    p = p + belt.into();
                } else {
                    break; // Ran off the end of a belt
                }

                // Error on loops
                if visited[toppings.bits][p.index(global.width)] {
                    return Err(format!("Loop detected at {p:?}"));
                } else {
                    visited[toppings.bits][p.index(global.width)] = true;
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
    // Use empty points, allowed to step on targets, and can follow belts + splitters
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
                    // If we're expanding along a belt, we have to follow the belt
                    if let Some(belt_direction) = self.belts[p.index(global.width)] {
                        if belt_direction != d {
                            continue;
                        }
                    }

                    // Don't expand the bfs out of bounds
                    let p2 = *p + d.into();
                    if !global.in_bounds(p2) {
                        continue;
                    }

                    // We can always step onto a belt
                    let is_belt = self.belts[p2.index(global.width)].is_some();

                    // TODO: This ignore direction for splitters for the time being
                    // This is technically correct, but could be optimized
                    let is_splitter = matches!(
                        global.entities[p2.index(global.width)],
                        Some(Entity {
                            kind: EntityKind::Splitter,
                            ..
                        })
                    );

                    // We can only step onto a splitter in the proper direction
                    let is_target = matches!(
                        global.entities[p2.index(global.width)],
                        Some(Entity {
                            kind: EntityKind::Target(_),
                            ..
                        })
                    );
                    let is_target_at_dst = p2 == dst && is_target;

                    if self.is_empty(global, p2) || is_belt || is_splitter || is_target_at_dst {
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
    #[tracing::instrument(skip(self, global), fields(belts = %self), ret)]
    fn is_valid(&self, global: &Global) -> bool {
        if global.use_tickwise {
            self.simulate_tickwise(global).is_ok()
        } else {
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
    }

    #[tracing::instrument(skip(self, global), fields(belts = %self), ret)]
    fn is_solved(&self, global: &Global) -> bool {
        if global.use_tickwise {
            let simulation_result = match self.simulate_tickwise(global) {
                Ok(simulation_result) => simulation_result,
                Err(e) => {
                    tracing::debug!("Simulation failed: {e}");
                    return false;
                }
            };
            tracing::debug!("Simulation result: {simulation_result:#?}");

            // All sources must have been delivered from
            let all_delivered_from = simulation_result
                .deliveries
                .values()
                .flatten()
                .map(|(p, _)| *p)
                .collect::<HashSet<_>>();

            for (index, entity) in global.entities.iter().enumerate() {
                if let Some(Entity {
                    kind: EntityKind::Source,
                    ..
                }) = entity
                {
                    let p = Point {
                        x: index as isize % global.width,
                        y: index as isize / global.width,
                    };

                    if !all_delivered_from.contains(&p) {
                        tracing::debug!("Source at {index} was not delivered from");
                        return false;
                    }
                }
            }

            // All targets must have been delivered to
            for (index, entity) in global.entities.iter().enumerate() {
                if let Some(Entity {
                    kind: EntityKind::Target(_),
                    ..
                }) = entity
                {
                    let p = Point {
                        x: index as isize % global.width,
                        y: index as isize / global.width,
                    };

                    if !simulation_result.deliveries.contains_key(&p) {
                        tracing::debug!("Target at {index} was not delivered to");
                        return false;
                    }
                }
            }

            // If a global target list is set, check that the donuts match
            if let Some(global_targets) = &global.targets {
                // For this check, a single source delivering to multiple targets counts as a single donut
                // That's why we HashSet first
                let mut all_delivered_donuts = simulation_result
                    .deliveries
                    .values()
                    .flatten()
                    .map(|pt| pt)
                    .collect::<HashSet<_>>()
                    .into_iter()
                    .map(|(_, t)| Some(*t))
                    .collect::<Vec<_>>();

                all_delivered_donuts.sort();

                if &all_delivered_donuts != global_targets {
                    tracing::debug!("Invalid solution, types don't match. Got {all_delivered_donuts:?} but needed {global_targets:?}");
                    return false;
                }
            }

            // We passed all conditions, we're SOLVED!
            true
        } else {
            // Simulate the current state
            tracing::debug!("Running simulation");
            tracing::debug!("\n{}", self.stringify(global));
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

                    // If we would continue with turn left, return turn right as facing and check it below
                    // If we wouldn't, return turn left and we'll pass that check below too
                    if self.is_empty(global, p + splitter_facing.turn_left().into()) {
                        splitter_facing.turn_left()
                    } else if self.is_empty(global, p + splitter_facing.turn_right().into()) {
                        splitter_facing.turn_right()
                    } else {
                        // Both branches of the splitter are already filled in
                        continue;
                    }
                } else {
                    unreachable!()
                };
                let p2 = p + facing.into();

                if !self.is_empty(global, p2) {
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
        // TODO: This uses the non-tickwise solver in both modes, is this okay?

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
                        EntityKind::Bumper(_) => match entity.facing {
                            Direction::Up => '┴',
                            Direction::Down => '┬',
                            Direction::Left => '┤',
                            Direction::Right => '├',
                        },
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
        if solver.states_checked() % 10_000 == 0 {
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
