use std::fmt::Display;
use std::fs::OpenOptions;
use std::io::{Error, Write};

use crate::canvas::Canvas;
use crate::color::Color;

pub struct PpmHeader {
    pub magic_number: u8,
    pub width: usize,
    pub height: usize,
    maximum_color: u8,
}

impl PpmHeader {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            magic_number: 3,
            width,
            height,
            maximum_color: 255,
        }
    }
}

impl Display for PpmHeader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let first_line = format!("P{}", self.magic_number);
        let second_line = format!("{} {}", self.width, self.height);
        write!(f, "{}\n{}\n{}", first_line, second_line, self.maximum_color)
    }
}

pub struct PpmFile {
    header: PpmHeader,
    body: Vec<Vec<Color>>,
}

impl PpmFile {
    pub fn max_color(&self) -> u8 {
        self.header.maximum_color
    }

    pub fn from_canvas(canvas: &Canvas) -> Self {
        let header = PpmHeader::new(canvas.width, canvas.height);
        Self {
            header,
            body: canvas.pixels.clone()
        }
    }

    pub fn save(&self, path: &str) -> Result<(), Error> {
        let mut file = OpenOptions::new()
            .write(true)
            .create(true)
            .open(path)
            .expect("Cannot write to file");

        let header = self.header.to_string();
        let body = self.prepare_body();
        writeln!(file, "{}", header).expect("Cannot write to file");
        writeln!(file, "{}", body)
    }

    // render canvas as Vector of strings
    pub fn prepare_body(&self) -> String {
        let strings = self
            .body
            .clone()
            .into_iter()
            .flat_map(|row| self.render_as_strings(row))
            .collect::<Vec<String>>();
        strings.join("\n")
    }

    // render canvas row as set of strings, each no more than 70 chars
    fn render_as_strings(&self, row: Vec<Color>) -> Vec<String> {
        let mut i = 0;
        let mut row_as_multiline = vec![];
        'outer: loop {
            let mut row_as_string: String = "".into();
            while row_as_string.len() < 65 {
                row_as_string.push_str(&self.render_pixel_as_string(row[i]));

                i += 1;
                if i == row.len() {
                    row_as_multiline.push(row_as_string);
                    break 'outer;
                }
            }
            row_as_multiline.push(row_as_string);
        }
        row_as_multiline
    }

    fn render_pixel_as_string(&self, c: Color) -> String {
        format!(
            "{} {} {} ",
            self.scale(c.red()),
            self.scale(c.green()),
            self.scale(c.blue())
        )
    }

    fn scale(&self, value: f64) -> u8 {
        match value {
            x if x < 0.0 => 0,
            x if x > 1.0 => self.max_color(),
            x => (x * (self.max_color() as f64)).floor() as u8,
        }
    }
}
