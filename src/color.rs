use crate::tuple::Tuple;
use std::ops::Mul;

type Color = Tuple<f64>;

impl Color {
    pub fn red(&self) -> f64 {
        self.x
    }
    pub fn green(&self) -> f64 {
        self.y
    }
    pub fn blue(&self) -> f64 {
        self.z
    }
}

impl Mul for Color {
    type Output = Color;
    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.red() * rhs.red(),
            self.green() * rhs.green(),
            self.blue() * rhs.blue(),
            0.0
        )
    }
}

pub fn color(red: f64, green: f64, blue: f64) -> Color {
    Color::new(red, green, blue, 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_colors() {
        let c = color(-0.5, 0.4, 1.7);

        assert_eq!(-0.5, c.red());
        assert_eq!(0.4, c.green());
        assert_eq!(1.7, c.blue());
    }

    #[test]
    fn test_adding_colors() {
        let c1 = color(0.9, 0.6, 0.75);
        let c2 = color(0.7, 0.1, 0.25);
        assert_relative_eq!(color(1.6, 0.7, 1.0), c1 + c2);
    }

    #[test]
    fn test_substracting_colors() {
        let c1 = color(0.9, 0.6, 0.75);
        let c2 = color(0.7, 0.1, 0.25);
        assert_relative_eq!(color(0.2, 0.5, 0.5), c1 - c2);
    }

    #[test]
    fn test_multiplying_color_by_a_scalar() {
        let c = color(0.2, 0.3, 0.4);
        assert_relative_eq!(color(0.4, 0.6, 0.8), c * 2);
    }

    #[test]
    fn test_multiplying_colors() {
        let c1 = color(1.0, 0.2, 0.4);
        let c2 = color(0.9, 1.0, 0.1);
        assert_relative_eq!(color(0.9, 0.2, 0.04), c1 * c2);
    }
}
