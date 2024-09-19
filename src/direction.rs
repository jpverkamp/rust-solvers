use anyhow::{anyhow, Result};

use crate::point::Point;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Direction {
    Up,
    Down,
    Left,
    Right,
}

impl From<Direction> for Point {
    fn from(direction: Direction) -> Point {
        match direction {
            Direction::Up => Point { x: 0, y: -1 },
            Direction::Down => Point { x: 0, y: 1 },
            Direction::Left => Point { x: -1, y: 0 },
            Direction::Right => Point { x: 1, y: 0 },
        }
    }
}

impl From<Point> for Direction {
    fn from(point: Point) -> Direction {
        match point {
            Point { x: 0, y: -1 } => Direction::Up,
            Point { x: 0, y: 1 } => Direction::Down,
            Point { x: -1, y: 0 } => Direction::Left,
            Point { x: 1, y: 0 } => Direction::Right,
            _ => panic!("Invalid point: {:?}", point),
        }
    }
}

impl TryFrom<&str> for Direction {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        match value {
            "up" => Ok(Direction::Up),
            "down" => Ok(Direction::Down),
            "left" => Ok(Direction::Left),
            "right" => Ok(Direction::Right),
            _ => Err(anyhow!("Invalid direction: {value}")),
        }
    }
}

impl TryFrom<char> for Direction {
    type Error = anyhow::Error;

    fn try_from(value: char) -> Result<Self> {
        match value {
            'U' | 'u' | 'N' | 'n' => Ok(Direction::Up),
            'D' | 'd' | 'S' | 's' => Ok(Direction::Down),
            'L' | 'l' | 'W' | 'w' => Ok(Direction::Left),
            'R' | 'r' | 'E' | 'e' => Ok(Direction::Right),
            _ => Err(anyhow!("Invalid direction: {value}")),
        }
    }
}

impl Direction {
    pub fn all() -> Vec<Direction> {
        vec![
            Direction::Up,
            Direction::Down,
            Direction::Left,
            Direction::Right,
        ]
    }

    pub fn flip(&self) -> Direction {
        match self {
            Direction::Up => Direction::Down,
            Direction::Down => Direction::Up,
            Direction::Left => Direction::Right,
            Direction::Right => Direction::Left,
        }
    }
}
