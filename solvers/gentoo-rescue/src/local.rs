use crate::Global;
use crate::global::{Critter, Tile, WallKind};
use direction::Direction;
use point::Point;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Local {
    pub(crate) critters: Vec<Critter>,
}

impl Local {
    // Try to move the given critter in the given direction
    // Returns the point the critter moves to (if it moves)
    #[tracing::instrument(skip(self, global), ret)]
    pub(crate) fn try_move(
        &self,
        global: &Global,
        index: usize,
        direction: Direction,
    ) -> Option<Point> {
        let mut pt = self.critters[index].location;
        let mut moved = false;

        loop {
            // If we're on water, stop moving
            if global.tile_at(pt) == Tile::Water {
                tracing::debug!("{pt:?} stopped at water");
                break;
            }

            // Bumped into a wall
            if global.wall_at(pt, direction) != WallKind::Empty {
                tracing::debug!("{pt:?} stopped at wall");
                break;
            }

            // Bumped into any other critter
            if self
                .critters
                .iter()
                .any(|c| c.location == pt + direction.into())
            {
                tracing::debug!("{pt:?} stopped at critter");
                break;
            }

            pt = pt + direction.into();
            tracing::debug!("moved to {pt:?}");
            moved = true;
        }

        // If we didn't move, this is invalid location
        // TODO: Handle bouncing etc
        if !moved {
            return None;
        }

        Some(pt)
    }
}

impl std::fmt::Display for Local {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for critter in &self.critters {
            let Critter {
                kind,
                color,
                location,
            } = critter;
            f.write_fmt(format_args!("{color:?} {kind:?} as {location:?}\n"))?
        }
        Ok(())
    }
}
