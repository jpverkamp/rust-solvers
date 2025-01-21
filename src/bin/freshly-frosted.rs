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
    TeleporterIn,
    TeleporterOut,
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

            // Teleporters (in doesn't currently have a facing)
            &['T', '-'] => Ok(Entity {
                kind: EntityKind::TeleporterIn,
                facing: Direction::default(),
            }),
            &['T', '+', facing] => Ok(Entity {
                kind: EntityKind::TeleporterOut,
                facing: dir(facing)?,
            }),

            // Something we don't know how to parse
            _ => Err(format!("Invalid entity: {value}")),
        }
    }
}

#[derive(Debug, Clone, Default)]
struct Global {
    width: isize,
    height: isize,
    entities: Vec<Option<Entity>>,
    targets: Option<Vec<Option<Toppings>>>,
    initial_belts: Vec<Option<Direction>>,
    use_tickwise: bool,
    allow_invalid_deliveries: bool,
    loop_threshold: Option<usize>,
}

impl Global {
    fn in_bounds(&self, p: Point) -> bool {
        p.x >= 0 && p.x < self.width && p.y >= 0 && p.y < self.height
    }
}

impl From<&str> for Global {
    fn from(input: &str) -> Self {
        let mut global = Global::default();

        for line in input.lines() {
            let mut line_width = 0;
            let line = line.trim();

            if line.is_empty() {
                continue;
            }

            if line.starts_with(':') {
                if line.starts_with(":target") {
                    global.targets = Some(
                        line.split_whitespace()
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
                } else if line.starts_with(":comment") {
                    // Not stored, just do nothing
                } else if line.starts_with(":tickwise") {
                    global.use_tickwise = true;
                } else if line.starts_with(":allow-invalid-deliveries") {
                    global.allow_invalid_deliveries = true;
                } else if line.starts_with(":loop-threshold") {
                    global.loop_threshold = Some(
                        line.split_whitespace()
                            .skip(1)
                            .next()
                            .expect("Missing loop threshold")
                            .parse()
                            .expect("Invalid loop threshold"),
                    );
                } else {
                    panic!("Invalid/unknown flag: {line}");
                }
                continue;
            }

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

                global.initial_belts.push(belt);
                global.entities.push(entity);
            }

            global.width = global.width.max(line_width);
            global.height += 1;
        }

        global
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
    extras: Vec<(Point, Toppings)>,
}

impl Local {
    fn is_empty(&self, global: &Global, p: Point) -> bool {
        if !global.in_bounds(p) {
            return false;
        }

        self.belts[p.index(global.width)].is_none()
            && global.entities[p.index(global.width)].is_none()
    }

    #[tracing::instrument(skip(self, global), ret)]
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
            waiting_time: usize,
            split_next_right: bool,
        }

        #[derive(Debug, Clone, PartialEq, Eq, Hash)]
        struct Update {
            move_from: Point,
            move_to: Point,
            toppings: Toppings,
            source: Point,
        }

        let mut state = vec![TileState::default(); vec_size];

        let mut deliveries = HashMap::new();
        let mut states_seen = HashMap::new();

        macro_rules! state_at {
            ($p:expr) => {
                state[$p.index(global.width)]
            };
        }

        // We have a very expensive tracing option; so don't calculate it if we're not going to print it
        let tracing_enabled = std::env::var("FRESHLY_FROSTED_TRACE").is_ok();
        let step_tracing_enabled = std::env::var("FRESHLY_FROSTED_STEP_TRACE").is_ok();

        // Advance the simulation one tick
        'tick: loop {
            let mut updates = vec![];
            let max_waiting_time = state.iter().map(|s| s.waiting_time).max().unwrap_or(0);

            // This shouldn't generally happen; but if it does we have an infinite loop so don't hang
            if max_waiting_time > global.width as usize * global.height as usize {
                panic!("Waiting time exceeded maximum");
            }

            // Debugging ticking
            if tracing_enabled {
                if step_tracing_enabled {
                    let initial_map = self.stringify(global).chars().collect::<Vec<_>>();
                    let mut maps = [
                        initial_map.clone(), // Map 0: Default
                        initial_map.clone(), // Map 1: Current toppings
                        initial_map.clone(), // Map 2: Waiting times
                        initial_map.clone(), // Map 3: Splitters
                    ];

                    for y in 0..global.height {
                        for x in 0..global.width {
                            // Map 0 does nothing

                            // Map 1 shows the toppings
                            // Map 2 shows the wait times (only for donuts)
                            let p = Point { x, y };
                            if let Some(toppings) = state_at!(p).toppings {
                                maps[1][p.index(global.width + 1)] =
                                    toppings.bits.to_string().chars().next().unwrap();

                                maps[2][p.index(global.width + 1)] = state_at!(p)
                                    .waiting_time
                                    .to_string()
                                    .chars()
                                    .next()
                                    .unwrap();
                            } else {
                                maps[1][p.index(global.width + 1)] = '.';
                                maps[2][p.index(global.width + 1)] = '.';
                            }

                            // Map 3 shows current splitter state
                            if let Some(Entity {
                                kind: EntityKind::Splitter,
                                ..
                            }) = global.entities[p.index(global.width)]
                            {
                                if state_at!(p).split_next_right {
                                    maps[3][p.index(global.width + 1)] = 'R';
                                } else {
                                    maps[3][p.index(global.width + 1)] = 'L';
                                }
                            } else {
                                maps[3][p.index(global.width + 1)] = '.';
                            }
                        }
                    }

                    // We want to render them side by side

                    // Convert into a string
                    let maps = maps
                        .iter()
                        .map(|m| m.iter().collect::<String>())
                        .collect::<Vec<_>>();

                    // Convert into a list of lines
                    let maps = maps
                        .iter()
                        .map(|m| m.split("\n").collect::<Vec<_>>())
                        .collect::<Vec<_>>();

                    // Combine them line by line
                    let mut final_map = String::new();
                    for y in 0..global.height {
                        for map in maps.iter() {
                            final_map.push_str(map[y as usize]);
                            final_map.push_str("   ");
                        }
                        final_map.push('\n');
                    }

                    let cache_size = states_seen.len();

                    let max_cache_value = states_seen
                        .iter()
                        .map(|(_, v)| *v)
                        .max()
                        .unwrap_or(0);

                    tracing::debug!(
                        "\
=== Starting tick ===
Deliveries: {deliveries:?}
States seen: {cache_size} (max: {max_cache_value})
Max waiting time: {max_waiting_time}

Maps (belts, toppings, waiting times, splitters):
{final_map}
",
                    );
                }
            }

            // Cache which exact states we've seen; break once we see the same more than once (or a set threshold of times)
            if let Some(count) = states_seen.get(&state) {
                let threshold = global.loop_threshold.unwrap_or(1);
                if *count >= threshold {
                    tracing::debug!("Loop detected (threshold={threshold}), breaking");
                    tracing::debug!("Deliveries are: {deliveries:?}");
                    break 'tick;
                }
            }
            states_seen.entry(state.clone()).and_modify(|e| *e += 1).or_insert(1);

            // Calculate all requested updates
            for x in 0..global.width {
                'next_point: for y in 0..global.height {
                    let p = Point { x, y };

                    // Create a potential update for sources
                    if let Some(Entity {
                        kind: EntityKind::Source,
                        facing,
                    }) = global.entities[p.index(global.width)]
                    {
                        updates.push(Update {
                            move_from: p,
                            move_to: p + facing.into(),
                            toppings: Toppings::none(),
                            source: p,
                        });
                        continue 'next_point;
                    }

                    // No updates for spaces that do not have a donut
                    if state_at!(p).toppings.is_none() {
                        continue;
                    }

                    // Bumpers take priority over belts
                    // TODO: Can you push into a merge situation? Then we'd need both.
                    // TODO: Assume there's only one bumper per space
                    for d in Direction::all() {
                        let p2 = p - d.into();
                        if !global.in_bounds(p2) {
                            continue;
                        }

                        if let Some(Entity {
                            kind: EntityKind::Bumper(bumper_toppings),
                            facing,
                        }) = global.entities[p2.index(global.width)]
                        {
                            if facing == d && bumper_toppings == state_at!(p).toppings.unwrap() {
                                tracing::warn!("Bumper at {p2:?} with facing {d:?} and toppings {bumper_toppings:?}");
                                updates.push(Update {
                                    move_from: p,
                                    move_to: p + d.into(),
                                    toppings: state_at!(p).toppings.unwrap(),
                                    source: state_at!(p).source.unwrap(),
                                });
                                continue 'next_point;
                            }
                        }
                    }

                    // Belts are easy!
                    if let Some(belt) = self.belts[p.index(global.width)] {
                        // We have to be able to move onto that space
                        let p2 = p + belt.into();
                        if !global.in_bounds(p2) {
                            return Err(format!("Attempted to move out of bounds at {p:?}"));
                        }

                        // Some entities cannot be moved onto (at all or in a specific direction)
                        if let Some(entity) = global.entities[p2.index(global.width)] {
                            match entity.kind {
                                EntityKind::Block
                                | EntityKind::Source
                                | EntityKind::Topper(_)
                                | EntityKind::Bumper(_)
                                | EntityKind::TeleporterOut => {
                                    return Err(format!(
                                        "Donut at {p:?} tried to move onto a {:?}",
                                        entity.kind
                                    ));
                                }
                                EntityKind::Target(_) | EntityKind::Splitter => {
                                    if entity.facing != belt {
                                        return Err(format!("Donut at {p:?} tried to move onto a {:?} facing the wrong way", entity.kind));
                                    }
                                }
                                EntityKind::TeleporterIn => {
                                    // No facing, we can always move onto this
                                }
                            }
                        }

                        // If we make it this far, this is a valid potential move
                        updates.push(Update {
                            move_from: p,
                            move_to: p + belt.into(),
                            toppings: state_at!(p).toppings.unwrap(),
                            source: state_at!(p).source.unwrap(),
                        });
                        continue 'next_point;
                    }

                    // Any entities we could be standing on
                    // We should never have moved onto invalid ones (see above)
                    if let Some(entity) = global.entities[p.index(global.width)] {
                        match entity.kind {
                            EntityKind::Block | EntityKind::Topper(_) | EntityKind::Bumper(_) | EntityKind::TeleporterOut => {
                                unreachable!("Donut at {p:?} is on a {:?}", entity.kind);
                            }
                            // Try to create a (potential) new donut
                            EntityKind::Source => unreachable!("Sources are handled earlier"),
                            // If we're on a target, matching done/not error
                            EntityKind::Target(toppings) => {
                                if toppings.is_none()
                                    || state_at!(p).toppings.unwrap() == toppings.unwrap()
                                    || global.allow_invalid_deliveries
                                {
                                    deliveries.entry(p).or_insert_with(HashSet::new).insert((
                                        state_at!(p).source.unwrap(),
                                        state_at!(p).toppings.unwrap(),
                                    ));
                                } else {
                                    return Err(format!(
                                        "Donut at {p:?} is on a target with the wrong toppings"
                                    ));
                                }
                            }
                            // Splitters try to move in the next direction
                            EntityKind::Splitter => {
                                if state_at!(p).split_next_right {
                                    updates.push(Update {
                                        move_from: p,
                                        move_to: p + entity.facing.turn_right().into(),
                                        toppings: state_at!(p).toppings.unwrap(),
                                        source: state_at!(p).source.unwrap(),
                                    });
                                } else {
                                    updates.push(Update {
                                        move_from: p,
                                        move_to: p + entity.facing.turn_left().into(),
                                        toppings: state_at!(p).toppings.unwrap(),
                                        source: state_at!(p).source.unwrap(),
                                    });
                                }
                            }
                            // On a teleporter, try to teleport
                            EntityKind::TeleporterIn => {
                                todo!();
                            }
                        }
                    }
                }
            }

            // Okay, now for any update that has multiple choices, we have to choose one
            // Choose the one that has the largest waiting_time
            // Then we have to increment the waiting time for the rest and wind back any updates depending on those
            let mut will_update = vec![true; updates.len()];
            for x in 0..global.width {
                for y in 0..global.height {
                    let p = Point { x, y };

                    let mut updates = updates
                        .iter()
                        .enumerate()
                        .filter(|(_, u)| u.move_to == p)
                        .collect::<Vec<_>>();

                    if updates.len() <= 1 {
                        continue;
                    }

                    // This is a kludge for 6/8; it requires that Bottom goes before Left
                    // So we'll initially sort bottom up and then because they're tied, it will win
                    updates.reverse();

                    // Sort, the longest waiting will end up first
                    updates.sort_by(|(_, a), (_, b)| {
                        state_at!(b.move_from)
                            .waiting_time
                            .cmp(&state_at!(a.move_from).waiting_time)
                    });

                    // The waiting time for the source of the one that moves is 0'ed
                    // The rest incremented
                    state[updates[0].1.move_from.index(global.width)].waiting_time = 0;
                    for (_, u) in updates.iter().skip(1) {
                        state[u.move_from.index(global.width)].waiting_time += 1;
                    }

                    // All the rest of them are flagged as not updating
                    for (index, _) in updates.iter().skip(1) {
                        will_update[*index] = false;
                    }
                }
            }

            // Now propagate that backwards
            'still_did_not_updating: loop {
                for (i, ui) in updates.iter().enumerate() {
                    for (j, uj) in updates.iter().enumerate() {
                        if i != j && will_update[i] && !will_update[j] && ui.move_to == uj.move_from
                        {
                            will_update[i] = false;
                            continue 'still_did_not_updating;
                        }
                    }
                }

                // If we make it through the loops without updating anything, we're done
                break;
            }

            // Now apply each update that is still in the list

            // First, remove from the source to make space for the destinations
            for (i, u) in updates.iter().enumerate() {
                if will_update[i] {
                    state_at!(u.move_from).toppings = None;
                    state_at!(u.move_from).source = None;

                    // Moving from a splitter means it toggles
                    // This updates the flag for everything, even non-splitters, but we never read it otherwise
                    state_at!(u.move_from).split_next_right =
                        !state_at!(u.move_from).split_next_right;
                }
            }

            // And then set each destination
            for (i, u) in updates.iter().enumerate() {
                if will_update[i] {
                    state_at!(u.move_to).toppings = Some(u.toppings);
                    state_at!(u.move_to).source = Some(u.source);

                    // Moving to a splitter does *not* toggle it
                }
            }

            // Finally, apply toppers
            for x in 0..global.width {
                for y in 0..global.height {
                    let p = Point { x, y };

                    if let Some(Entity {
                        kind: EntityKind::Topper(topping),
                        facing,
                    }) = global.entities[p.index(global.width)]
                    {
                        let p2 = p + facing.into();
                        if !global.in_bounds(p2) {
                            continue;
                        }

                        if let Some(state) = state.get_mut(p2.index(global.width)) {
                            if let Some(toppings) = state.toppings {
                                // We have to have exactly all previous toppings; if so apply the new one
                                let must_have = Toppings::from(topping.bits() - 1);
                                if toppings & must_have == must_have {
                                    state.toppings = Some(toppings | topping);
                                }
                            }
                        }
                    }
                }
            }
        }

        // At the end, extra donuts are any on an empty space
        let extras = state
            .iter()
            .enumerate()
            .filter_map(|(index, s)| {
                let p = Point {
                    x: index as isize % global.width,
                    y: index as isize / global.width,
                };

                if self.is_empty(global, p) && s.toppings.is_some() {
                    Some((p, s.toppings.unwrap()))
                } else {
                    None
                }
            })
            .collect();

        let result = SimulateTickwiseResult { deliveries, extras };

        // Cache the result
        SIMULATE_TICKWISE_CACHE
            .lock()
            .unwrap()
            .insert(self.clone(), result.clone());

        Ok(result)
    }

    #[tracing::instrument(skip(self, global), ret)]
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

                            // Already applied, not an error but we don't want to trace it
                            if toppings | new_toppings == toppings {
                                continue;
                            }

                            // Add the new topping!
                            tracing::debug!(
                                "Adding topping {new_toppings:?} from {p2:?} / {facing:?}"
                            );
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
                        if bump_d == facing && toppings == bump_toppings {
                            p = p + bump_d.into();
                            break;
                        }
                    }
                }

                // Move move along the belt
                if let Some(belt) = self.belts[p.index(global.width)] {
                    p = p + belt.into();

                    // Some entities cannot be run into at all
                    // Some require that you are moving the right direction
                    if let Some(entity) = global.entities[p.index(global.width)] {
                        match entity.kind {
                            EntityKind::Block
                            | EntityKind::Source
                            | EntityKind::Topper(_)
                            | EntityKind::Bumper(_)
                            | EntityKind::TeleporterOut => {
                                return Err(format!("Ran into a {:?} at {p:?}", entity.kind))
                            }
                            EntityKind::Target(_) | EntityKind::Splitter => {
                                if entity.facing != belt {
                                    return Err(format!("Ran into a {:?} at {p:?} facing the wrong way (expected {:?}, got {belt:?})", entity.kind, entity.facing));
                                }
                            }
                            EntityKind::TeleporterIn => {
                                // No facing, we can always move onto this
                            }
                        }
                    }
                } else {
                    break; // Ran off the end of a belt
                }

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

            // If we're allowing invalid deliveries, it's still not solved if there are any that don't match
            if global.allow_invalid_deliveries {
                for (dst, donuts) in simulation_result.deliveries.iter() {
                    for (_, toppings) in donuts.iter() {
                        match global.entities[dst.index(global.width)] {
                            // Each destination must be none or matching
                            Some(Entity {
                                kind: EntityKind::Target(target_toppings),
                                ..
                            }) if target_toppings.is_none()
                                || *toppings == target_toppings.unwrap() => {}
                            // If not, this is not a valid solution
                            _ => {
                                tracing::debug!(
                                    "Invalid delivery at {dst:?} with toppings {toppings:?}"
                                );
                                return false;
                            }
                        }
                    }
                }
            }

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
                        tracing::debug!("Source at {p:?} was not delivered from");
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
                        tracing::debug!("Target at {p:?} was not delivered to");
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

            // Each target must have at least one donut
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

                    if !donuts.iter().any(|(donut_p, _)| *donut_p == p) {
                        tracing::debug!("No donut delivered to {p:?}");
                        return false;
                    }
                }
            }

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

    #[allow(unreachable_code)]
    #[tracing::instrument(skip(self, global), fields(belts = %self))]
    fn next_states(&self, global: &Global) -> Option<Vec<(i64, (), Local)>> {
        let mut donuts = if global.use_tickwise {
            match self.simulate_tickwise(global) {
                Ok(tickwise_result) => tickwise_result.extras,
                Err(e) => {
                    tracing::debug!("Simulation failed: {e}");
                    return None;
                }
            }
        } else {
            match self.simulate(global) {
                Ok(donuts) => donuts,
                Err(e) => {
                    tracing::debug!("Simulation failed: {e}");
                    return None;
                }
            }
        };
        donuts.sort();

        // Find the first empty point
        for (p, _) in donuts {
            tracing::debug!("Checking for expansion at {p:?}");

            if !self.is_empty(global, p) {
                continue;
            }

            tracing::debug!("Expanding new_state at {p:?}");

            // This is necessary in the (rare, only 6/8 so far) case that one expansion is invalid now but will become valid after another state is added
            let maybe_next_states = Direction::all()
                .iter()
                .flat_map(|&d| {
                    if global.in_bounds(p + d.into()) {
                        let mut new_state = self.clone();
                        new_state.belts[p.index(global.width)] = Some(d);
                        Some((1, (), new_state))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();

            if !maybe_next_states.is_empty() {
                return Some(maybe_next_states);
            }
        }

        None
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
                        EntityKind::TeleporterIn => '○',
                        EntityKind::TeleporterOut => match entity.facing {
                            Direction::Up => '◒',
                            Direction::Down => '◓',
                            Direction::Left => '◑',
                            Direction::Right => '◐',
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

#[cfg(test)]
mod freshly_frosted_tests {
    use super::*;
    use Direction::{Down, Left, Right, Up};

    macro_rules! test_next_states_inner {
        ($name:ident, $source:expr, $expected_states:expr) => {
            // Initial stringify
            let global = Global::from($source);
            let local = global.make_local();
            let mut failures = vec![];

            match local.next_states(&global) {
                None => failures.push("No next states found".to_owned()),
                Some(next_states) => {
                    let expected_len = $expected_states
                        .iter()
                        .map(|(_, ds)| ds.len())
                        .sum::<usize>();

                    if next_states.len() != expected_len {
                        let next_states_stringy = next_states
                            .iter()
                            .map(|(_, _, state)| state.stringify(&global))
                            .collect::<Vec<_>>()
                            .join("\n\n");

                        failures.push(format!(
                            "Expected {} state(s), got {}:\n{}",
                            expected_len,
                            next_states.len(),
                            next_states_stringy
                        ));
                    }

                    for ((x, y), ds) in $expected_states.iter() {
                        let i = (y * global.width + x) as usize;
                        for d in ds.iter() {
                            if !next_states
                                .iter()
                                .any(|(_, _, state)| state.belts[i].is_some_and(|bd| bd == *d))
                            {
                                failures.push(format!(
                                    "({x}, {y}) + {d:?} not found",
                                    x = x,
                                    y = y,
                                    d = d
                                ));
                            }
                        }
                    }
                }
            }

            if !failures.is_empty() {
                println!("Failures:");
                for failure in failures {
                    println!("  {}", failure);
                }

                println!("Simulate results: {:#?}", local.simulate(&global));
                println!(
                    "Tickwise simulate results: {:#?}",
                    local.simulate_tickwise(&global)
                );

                panic!();
            }
        };
    }

    macro_rules! test_next_states {
        ($name:ident, $source:expr, $expected_states:expr) => {
            paste::paste! {
                #[test]
                fn [<test_next_states_ $name>] () {
                    test_next_states_inner!($name, $source, $expected_states);
                }

                #[test]
                fn [<test_next_states_tickwise_ $name>] () {
                    let source = format!(":tickwise\n{}", $source);
                    test_next_states_inner!($name, source.as_str(), $expected_states);
                }
            }
        };
    }

    test_next_states! {
        source,
        "
        .   .   .
        +>  .   .
        .   .   .
        ",
        [
            ((1, 1), [Up, Down, Left, Right]),
        ]
    }

    test_next_states! {
        single_belt,
        "
        .   .   .
        .   .   .
        +>  ^   .
        .   .   .
        ",
        [
            ((1, 1), [Up, Down, Left, Right]),
        ]
    }

    test_next_states! {
        splitter_left, // Note: Splitters should default left
        "
        .   .   .
        .   .   .
        +>  x>  .
        .   .   .
        ",
        [
            ((1, 1), [Up, Down, Left, Right]),
        ]
    }

    test_next_states! {
        splitter_right,
        "
        .   .   .
        .   -^? .
        +>  x>  .
        .   .   .
        .   .   .
        ",
        [
            ((1, 3), [Up, Down, Left, Right]),
        ]
    }

    test_next_states! {
        bumper,
        "
        .   .   .
        .   bv0 .
        +>  >   .
        .   .   .
        ",
        [
            ((1, 3), [Left, Right, Up]),
        ]
    }

    test_next_states! {
        non_matching_bumper,
        "
        .   .   .
        .   bv1 .
        +>  >   .
        .   .   .
        ",
        [
            ((2, 2), [Up, Down, Left]),
        ]
    }

    test_next_states! {
        bumper_loop,
        "
        1>  v   <   <
        .   v   bv1 ^
        +>  >   >   ^
        .   .   .   .
        ",
        [
            ((2, 3), [Up, Left, Right]),
        ]
    }

    test_next_states! {
        bumper_double_loop,
        "
        2>  v   <   <
        1>  v   bv3 ^
        +>  >   >   ^
        .   .   .   .
        ",
        [
            ((2, 3), [Up, Left, Right]),
        ]
    }
}
