use crate::Global;
use point::Point;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Local {
    pub(crate) worm: Vec<Point>,
    pub(crate) blocks: Vec<Vec<Point>>,
}

impl Global {
    pub(crate) fn make_local(&self) -> Local {
        // Start at the lower left 3 in length
        let worm = vec![
            Point {
                x: -1,
                y: self.height - 1,
            },
            Point {
                x: -2,
                y: self.height - 1,
            },
            Point {
                x: -3,
                y: self.height - 1,
            },
        ];

        // Create a single initial block
        let mut block = vec![];
        for y in 0..self.height {
            for x in 0..self.width {
                block.push(Point { x, y });
            }
        }
        let blocks = vec![block];

        Local { worm, blocks }
    }
}

impl std::fmt::Display for Local {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut result = String::new();
        result.push_str("Worm{");
        for (i, p) in self.worm.iter().enumerate() {
            if i > 0 {
                result.push(' ');
            }
            result.push_str(&format!("{},{}", p.x, p.y));
        }
        result.push_str("} ");

        result.push_str("Blocks{");
        for b in self.blocks.iter() {
            result.push_str(&format!(
                "[{}]",
                b.iter()
                    .map(|p| format!("{},{}", p.x, p.y))
                    .collect::<Vec<_>>()
                    .join(" ")
            ));
        }

        write!(f, "{}", result)
    }
}
