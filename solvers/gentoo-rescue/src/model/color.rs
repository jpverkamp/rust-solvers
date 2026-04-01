#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum Color {
    #[default]
    Gray,
    Red,
    Yellow,
    Green,
    Blue,
}

impl From<&str> for Color {
    fn from(value: &str) -> Self {
        match value {
            "gray" | "grey" => Color::Gray,
            "red" => Color::Red,
            "yellow" => Color::Yellow,
            "green" => Color::Green,
            "blue" => Color::Blue,
            _ => unimplemented!("Unknown color {value}"),
        }
    }
}

impl TryFrom<char> for Color {
    type Error = ();

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'r' => Ok(Color::Red),
            'y' => Ok(Color::Yellow),
            'g' => Ok(Color::Green),
            'b' => Ok(Color::Blue),
            'w' => Ok(Color::Gray),
            _ => Err(()),
        }
    }
}
