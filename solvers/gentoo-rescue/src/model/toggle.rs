use direction::Direction;
use point::Point;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum ToggleKind {
    Wall(Point, Direction),
    Floor(Point),
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ToggleRule {
    pub(crate) key: char,
    pub(crate) kind: ToggleKind,
}

impl From<&str> for ToggleRule {
    fn from(value: &str) -> Self {
        let parts: Vec<&str> = value.split_whitespace().collect();
        let key = parts[1].chars().next().unwrap();
        assert!(
            key.is_ascii_uppercase(),
            "toggle key must be an uppercase letter"
        );
        let row = parts[3].parse::<isize>().expect("Invalid row in toggle") - 1;
        let col = parts[4].parse::<isize>().expect("Invalid column in toggle") - 1;

        let kind = match parts[2] {
            "floor" => ToggleKind::Floor(Point { x: col, y: row }),
            "wall" => ToggleKind::Wall(
                Point { x: col, y: row },
                Direction::try_from(parts[5]).expect("Invalid direction in toggle"),
            ),
            _ => unimplemented!("unknown toggle kind {}", parts[2]),
        };
        ToggleRule { key, kind }
    }
}
