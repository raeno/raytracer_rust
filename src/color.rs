use approx::{AbsDiffEq, RelativeEq};
use num::zero;
use std::ops::{Add, Div, Mul, Sub};

use crate::tuple::Tuple;

#[derive(PartialEq, Clone, Debug, Copy)]
pub struct Color(Tuple<f64>);

impl Color {
    pub fn new(red: f64, green: f64, blue: f64) -> Self {
        Self(Tuple::new(red, green, blue, 0.0))
    }

    pub fn red(&self) -> f64 {
        self.0.x
    }

    pub fn green(&self) -> f64 {
        self.0.y
    }
    pub fn blue(&self) -> f64 {
        self.0.z
    }

    pub fn is_black(&self) -> bool {
        self.red() == zero() && self.green() == zero() && self.blue() == zero()
    }
}

impl Add for Color {
    type Output = Color;
    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}

impl Sub for Color {
    type Output = Color;
    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl Mul<Color> for Color {
    type Output = Color;
    fn mul(self, rhs: Color) -> Self::Output {
        Self::new(
            self.red() * rhs.red(),
            self.green() * rhs.green(),
            self.blue() * rhs.blue(),
        )
    }
}

impl Div<f64> for Color {
    type Output = Color;
    fn div(self, rhs: f64) -> Self::Output {
        Self(self.0 / rhs)
    }
}

impl Div<i32> for Color {
    type Output = Color;
    fn div(self, rhs: i32) -> Self::Output {
        Self(self.0 / rhs)
    }
}

impl Mul<f64> for Color {
    type Output = Color;
    fn mul(self, rhs: f64) -> Self::Output {
        Self(self.0 * rhs)
    }
}

impl Mul<i32> for Color {
    type Output = Color;
    fn mul(self, rhs: i32) -> Self::Output {
        Self(self.0 * rhs)
    }
}

impl AbsDiffEq for Color {
    type Epsilon = <Tuple<f64> as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        Tuple::abs_diff_eq(&self.0, &other.0, epsilon)
    }
}

impl RelativeEq for Color {
    fn default_max_relative() -> Self::Epsilon {
        f64::default_max_relative()
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        Tuple::relative_eq(&self.0, &other.0, epsilon, max_relative)
    }
}

pub fn color_i8(red: i8, green: i8, blue: i8) -> Color {
    Color::new(red.into(), blue.into(), green.into())
}

pub fn color(red: f64, green: f64, blue: f64) -> Color {
    Color::new(red, green, blue)
}

pub struct Colors;

impl Colors {
    pub fn black() -> Color {
        Color::new(0.0, 0.0, 0.0)
    }

    pub fn red() -> Color {
        Color::new(1.0, 0.0, 0.0)
    }

    pub fn green() -> Color {
        Color::new(0.0, 1.0, 0.0)
    }

    pub fn blue() -> Color {
        Color::new(0.0, 0.0, 1.0)
    }

    pub fn white() -> Color {
        Color::new(1.0, 1.0, 1.0)
    }
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
