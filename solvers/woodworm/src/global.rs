#[derive(Debug, Clone, Default)]
pub(crate) struct Global {
    // Map settings
    pub(crate) width: isize,
    pub(crate) height: isize,
    pub(crate) cells: Vec<bool>,
}

impl From<&str> for Global {
    fn from(input: &str) -> Self {
        let mut global = Global::default();

        for line in input.lines() {
            if line.trim().is_empty() || line.starts_with("//") {
                continue;
            }

            let chars = line.chars();
            let mut width = 0;
            global.height += 1;

            for c in chars {
                width += 1;
                if c == '#' {
                    global.cells.push(true);
                } else {
                    global.cells.push(false);
                }
            }

            if global.width == 0 {
                global.width = width;
            } else if global.width != width {
                panic!("Map width mismatch: {} != {}", global.width, width);
            }
        }

        global
    }
}
