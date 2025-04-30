use std::{env, fmt::Display};

use serde::{Deserialize, Serialize};
use solver::{Solver, State};

const ANTENNA_RADIUS: f32 = 2.0;

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
struct Map {
    nodes: Vec<Node>,

    #[serde(default)]
    mode_target: bool,

    #[serde(default)]
    mode_no_cross: bool,
}

#[derive(Copy, Clone, Debug, Deserialize, Serialize, Eq, PartialEq, Hash)]
struct Node {
    x: u8,
    y: u8,
    kind: Kind,

    #[serde(default)]
    size: u8,

    #[serde(default)]
    color: Color,

    #[serde(default)]
    target: bool,

    #[serde(default)]
    initial_signal: u8,
}

impl Node {
    fn distance(&self, other: &Node) -> f32 {
        ((other.x as f32 - self.x as f32).powi(2) + (other.y as f32 - self.y as f32).powi(2)).sqrt()
    }
}

impl Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut buf = String::new();

        buf.push_str(format!("{:?}", self.kind).as_str());

        match self.kind {
            Kind::Transmitter | Kind::Receiver | Kind::Transceiver => {
                buf.push_str(self.size.to_string().as_str())
            }
            Kind::Antenna => {}
        };

        buf.push_str(format!("({}, {})", self.x, self.y).as_str());

        if self.target {
            buf.push_str("+Target");
        }

        if self.color != Color::White {
            buf.push_str(format!("+{:?}", self.color).as_str());
        }

        write!(f, "{}", buf)
    }
}

#[derive(Copy, Clone, Debug, Deserialize, Serialize, Eq, PartialEq, Hash)]
enum Color {
    White,
    Yellow,
    Green,
}

impl Default for Color {
    fn default() -> Self {
        Color::White
    }
}

#[derive(Copy, Clone, Debug, Deserialize, Serialize, Eq, PartialEq, Hash)]
enum Kind {
    Transmitter,
    Receiver,
    Transceiver,
    Antenna,
}

#[derive(Clone, Debug, Hash, Eq, PartialEq)]
struct Energy {
    current: Vec<u8>,
    filled: Vec<u8>,
    steps: Vec<(usize, usize)>,
}

impl State<Map, ()> for Energy {
    fn next_states(&self, global: &Map) -> Option<Vec<(i64, (), Self)>>
    where
        Self: Sized,
    {
        let mut next = Vec::new();

        for src_index in 0..global.nodes.len() {
            for dst_index in 0..global.nodes.len() {
                // Cannot send to self
                if src_index == dst_index {
                    continue;
                }
                let src = global.nodes[src_index];
                let dst = global.nodes[dst_index];

                // Cannot send from some types (receivers, antennas)
                match src.kind {
                    Kind::Transmitter => {}
                    Kind::Receiver => continue,
                    Kind::Transceiver => {}
                    Kind::Antenna => continue,
                }

                // Cannot receive for some types (receivers)
                match dst.kind {
                    Kind::Transmitter => continue,
                    Kind::Receiver => {}
                    Kind::Transceiver => {}
                    Kind::Antenna => {}
                }

                // Cannot send to a different color
                // TODO: Transformers
                if src.color != dst.color {
                    continue;
                }

                // Cannot send if we don't have any energy
                if self.current[src_index] == 0 {
                    continue;
                }

                // Cannot send if the dst is full
                // Antennas as always sendable
                if dst.kind != Kind::Antenna && self.filled[dst_index] == dst.size {
                    continue;
                }

                // Cannot reuse or reverse a previous connection
                if self.steps.contains(&(src_index, dst_index))
                    || self.steps.contains(&(dst_index, src_index))
                {
                    continue;
                }

                // Cannot go through another node
                if global.nodes.iter().enumerate().any(|(other_index, other)| {
                    if src_index == other_index || dst_index == other_index {
                        return false;
                    }

                    src.distance(other) + other.distance(&dst) == src.distance(&dst)
                }) {
                    continue;
                }

                // Update energy with the new connection
                let mut new_state = self.clone();
                new_state.steps.push((src_index, dst_index));

                // Antennas are treated as maximum size
                let dst_size = if dst.kind == Kind::Antenna {
                    u8::MAX
                } else {
                    dst.size
                };

                let to_move = u8::min(
                    new_state.current[src_index],
                    dst_size - new_state.filled[dst_index],
                );

                new_state.current[src_index] -= to_move;
                new_state.current[dst_index] += to_move;
                new_state.filled[dst_index] += to_move;

                // Re-use loops/previous paths
                // + Broadcast from antennas
                'update_loops: loop {
                    // Update antennas with energy
                    for (i, src) in global.nodes.iter().enumerate() {
                        if new_state.current[i] > 0 && global.nodes[i].kind == Kind::Antenna {
                            // Broadcast to anyone not connected to us with capacity
                            for (j, dst) in global.nodes.iter().enumerate() {
                                // Don't self broadcast
                                if i == j {
                                    continue;
                                }

                                // TODO: Don't broadcast to transmitters?

                                // TODO: Antenna loops
                                if dst.kind == Kind::Antenna {
                                    continue;
                                }

                                // Don't broadcast to anyone outside of range
                                if src.distance(&dst) > ANTENNA_RADIUS {
                                    continue;
                                }

                                // Don't broadcast to anyone connected to me
                                if new_state.steps.contains(&(j, i)) {
                                    continue;
                                }

                                // Otherwise, fill up to maximum
                                let dst_size = if dst.kind == Kind::Antenna {
                                    u8::MAX
                                } else {
                                    dst.size
                                };

                                let to_broadcast =
                                    u8::min(new_state.current[i], dst_size - new_state.filled[j]);

                                // println!("broadcasting from {i}={src} to {j}={dst}, sending: {to_broadcast}, state is {new_state:?}");

                                if to_broadcast > 0 {
                                    new_state.current[j] += to_broadcast;
                                    new_state.filled[j] += to_broadcast;
                                }

                                // println!("after state is {new_state:?}");
                            }

                            // Now that we're done, remove the energy from the antenna
                            new_state.current[i] = 0;
                            continue 'update_loops;
                        }
                    }

                    // Otherwise, try to update using previous links
                    for (i, j) in new_state.steps.iter() {
                        // Antennas are treated as maximum size
                        let j_size = if global.nodes[*j].kind == Kind::Antenna {
                            u8::MAX
                        } else {
                            global.nodes[*j].size
                        };

                        let to_move = u8::min(new_state.current[*i], j_size - new_state.filled[*j]);

                        if to_move > 0 {
                            new_state.current[*i] -= to_move;
                            new_state.current[*j] += to_move;
                            new_state.filled[*j] += to_move;

                            continue 'update_loops;
                        }
                    }

                    break;
                }

                next.push((1, (), new_state));
            }
        }

        if next.len() == 0 {
            None
        } else {
            Some(next)
        }
    }

    fn is_valid(&self, _global: &Map) -> bool {
        return true;
    }

    fn is_solved(&self, global: &Map) -> bool {
        for index in 0..global.nodes.len() {
            if self.filled[index] < global.nodes[index].size {
                return false;
            }

            // If mode_target, we aren't solved unless there's >= 1 currently at each target
            if global.mode_target && global.nodes[index].target && self.current[index] == 0 {
                return false;
            }
        }

        // If mode_no_cross, we aren't solved if any two lines cross
        if global.mode_no_cross {
            for (i, line1) in self.steps.iter().enumerate() {
                for (j, line2) in self.steps.iter().enumerate() {
                    if i == j {
                        continue;
                    }

                    let pa = global.nodes[line1.0];
                    let pb = global.nodes[line1.1];
                    let pc = global.nodes[line2.0];
                    let pd = global.nodes[line2.1];

                    if pa == pc || pa == pd || pb == pc || pb == pd {
                        continue;
                    }

                    let a1 = pb.y as f32 - pa.y as f32;
                    let b1 = pa.x as f32 - pb.x as f32;
                    let c1 = a1 * (pa.x as f32) + b1 * (pa.y as f32);

                    let a2 = pd.y as f32 - pc.y as f32;
                    let b2 = pc.x as f32 - pd.x as f32;
                    let c2 = a2 * (pc.x as f32) + b2 * (pc.y as f32);

                    let determinant = a1 * b2 - a2 * b1;

                    if determinant == 0f32 {
                        continue;
                    }

                    let x = (b2 * c1 - b1 * c2) / determinant;
                    let y = (a1 * c2 - a2 * c1) / determinant;

                    let x_in_range = (pa.x as f32) < x && x < (pb.x as f32)
                        || (pb.x as f32) < x && x < (pa.x as f32);

                    let y_in_range = (pa.y as f32) < y && y < (pb.y as f32)
                        || (pb.y as f32) < y && y < (pa.y as f32);

                    if x_in_range && y_in_range {
                        // println!("{l1_1} -> {l1_2} crosses {l2_1} -> {l2_2}");

                        return false;
                    }
                }
            }
        }

        true
    }

    fn stringify(&self, global: &Map) -> String {
        format!("global: {global:?}, energy: {self:?}")
    }

    fn heuristic(&self, _global: &Map) -> i64 {
        10000 // TODO
    }
}

impl From<Map> for Energy {
    fn from(map: Map) -> Self {
        let initial = map
            .nodes
            .iter()
            .map(|node| node.initial_signal)
            .collect::<Vec<_>>();
        Energy {
            current: initial.clone(),
            filled: initial.clone(),
            steps: vec![],
        }
    }
}

fn main() {
    let stdin = std::io::stdin().lock();
    let mut map: Map = serde_json::from_reader(stdin).unwrap();

    let args = env::args().collect::<Vec<_>>();
    if args.contains(&"--target".to_string()) {
        map.mode_target = true;
    }
    tracing::info!("Target mode: {}", map.mode_target);
    if args.contains(&"--no-cross".to_string()) {
        map.mode_no_cross = true;
    }
    tracing::info!("No cross mode: {}", map.mode_no_cross);

    let energy = Energy::from(map.clone());
    let mut solver = Solver::new(map.clone(), energy.clone());

    while let Some(_state) = solver.next() {}

    let solution = solver.get_solution();

    // Output answer as text (if we have one)
    if let Some(path) = solution.clone() {
        for (i, (src, dst)) in path.clone().steps.into_iter().enumerate() {
            let i = i + 1;

            let src = map.nodes[src];
            let dst = map.nodes[dst];

            println!("[{i}] {src} -> {dst}");
        }
    }

    // Output answer as graph
    println!();
    println!("digraph G {{");
    for (i, node) in map.nodes.iter().enumerate() {
        let mut properties = vec![];

        // Set to the proper location
        properties.push(format!("pos=\"{},{}!\"", node.x, -(node.y as i32)));

        // Shape is determined by node type
        properties.push(format!(
            "shape=\"{}\"",
            match node.kind {
                Kind::Transmitter => "triangle",
                Kind::Receiver => "box",
                Kind::Transceiver =>
                    if node.target {
                        "doublecircle"
                    } else {
                        "circle"
                    },
                Kind::Antenna => "star",
            }
        ));

        // Label is determined by size and optionally initial signal
        properties.push(if node.initial_signal > 0 || node.size > 0 {
            if node.initial_signal > 0 {
                format!("label=\"{}/{}\"", node.initial_signal, node.size)
            } else {
                format!("label=\"{}\"", node.size)
            }
        } else {
            format!("label=\"\"")
        });

        properties.push(format!("style=\"filled\""));
        properties.push(format!("fillcolor=\"{:?}\"", node.color));

        println!("\tnode{i}[{}]", properties.join(" "));
    }

    // Write lines between attached nodes (if we have a solution)
    if let Some(path) = solution.clone() {
        println!();

        // Connections with the ordering
        for (i, (src, dst)) in path.clone().steps.into_iter().enumerate() {
            let i = i + 1;
            println!("\tnode{src} -> node{dst} [label=\"{i}\"]");
        }

        // Potential broadcast partners
        for (i, src) in map.nodes.iter().enumerate() {
            for (j, dst) in map.nodes.iter().enumerate() {
                if i == j {
                    continue;
                }

                if src.kind != Kind::Antenna {
                    continue;
                }

                // TODO: Antenna loops
                if dst.kind == Kind::Antenna {
                    continue;
                }

                // if path.steps.contains(&(j, i)) {
                //     continue;
                // }

                if src.distance(&dst) > ANTENNA_RADIUS {
                    continue;
                }

                println!("\tnode{i} -> node{j} [style=\"dashed\"]");
            }
        }
    }

    println!("}}");

    tracing::info!(
        "{} states, {} seconds",
        solver.states_checked(),
        solver.time_spent()
    );
}
