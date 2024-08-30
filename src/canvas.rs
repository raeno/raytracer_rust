use crate::color::{color, Color, Colors};

struct Canvas {
    width: usize,
    height: usize,
    pixels: Vec<Vec<Color>>,
}

impl Canvas {
    pub fn new(width: usize, height: usize) -> Self {
        let row = vec![color(0.0, 0.0, 0.0); width];
        let pixels = vec![row; height];

        Self {
            width,
            height,
            pixels,
        }
    }

    pub fn is_solid(&self, c: Color) -> bool {
        self.pixels
            .iter()
            .all(|row| row.iter().all(|pixel| *pixel == c))
    }

    pub fn write_pixel(&mut self, x: usize, y: usize, c: Color) {
        self.pixels[y][x] = c;
    }

    pub fn pixel_at(&self, x: usize, y: usize) -> Color {
        self.pixels[y][x]
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::color::color_i8;

    #[test]
    fn test_empty_canvas_has_black_pixels() {
        let canvas = Canvas::new(10, 20);

        assert_eq!(10, canvas.width);
        assert_eq!(20, canvas.height);

        let is_all_black = canvas.is_solid(Colors::black());

        assert_eq!(true, is_all_black)
    }

    #[test]
    fn test_writing_and_reading_pixel_to_canvas() {
        let mut canvas = Canvas::new(10, 20);
        let red = color_i8(1, 0, 0);

        assert_eq!(Colors::black(), canvas.pixel_at(5, 5));
        canvas.write_pixel(5, 5, red);
        assert_eq!(Colors::red(), canvas.pixel_at(5, 5));
    }
}
