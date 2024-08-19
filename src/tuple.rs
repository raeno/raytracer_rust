use num::{one, zero, Num};
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(PartialEq, Eq, Clone, Debug, Copy)]
pub struct Tuple<T> {
    x: T,
    y: T,
    z: T,
    w: T,
}

pub trait Magnitude {
    fn magnitude(&self) -> f64;
}

pub trait Normalize {
    type Output;
    fn normalize(&self) -> Self::Output;
}

pub trait Dot<Rhs = Self> {
    type Output;
    fn dot(&self, rhs: &Rhs) -> Self::Output;
}

pub trait CrossProduct<Rhs = Self> {
    type Output;
    fn cross(&self, rhs: &Rhs) -> Self::Output;
}

impl From<Tuple<i32>> for Tuple<f64> {
    fn from(tuple: Tuple<i32>) -> Self {
        Self::new(
            tuple.x.into(),
            tuple.y.into(),
            tuple.z.into(),
            tuple.w.into(),
        )
    }
}

impl<T: Num + Copy + Into<f64>> Tuple<T> {
    pub fn new(x: T, y: T, z: T, w: T) -> Self {
        Self { x, y, z, w }
    }

    pub fn is_point(&self) -> bool {
        let fw: f64 = self.w.into();
        fw == one()
    }

    pub fn is_vector(&self) -> bool {
        let fw: f64 = self.w.into();
        fw == zero()
    }
}

impl<T: Num + Copy + Into<f64>> Add<Tuple<T>> for Tuple<T> {
    type Output = Tuple<T>;

    fn add(self, rhs: Tuple<T>) -> Self::Output {
        Self::new(
            self.x + rhs.x,
            self.y + rhs.y,
            self.z + rhs.z,
            self.w + rhs.w,
        )
    }
}

impl<T: Num + Copy + Into<f64>> Sub<Tuple<T>> for Tuple<T> {
    type Output = Self;

    fn sub(self, rhs: Tuple<T>) -> Self::Output {
        Self::new(
            self.x - rhs.x,
            self.y - rhs.y,
            self.z - rhs.z,
            self.w - rhs.w,
        )
    }
}

impl<T: Num + Copy + Into<f64>> Mul<T> for Tuple<T> {
    type Output = Tuple<T>;
    fn mul(self, rhs: T) -> Self::Output {
        Tuple::new(self.x * rhs, self.y * rhs, self.z * rhs, self.w * rhs)
    }
}

impl Mul<f64> for Tuple<i32> {
    type Output = Tuple<f64>;
    fn mul(self, rhs: f64) -> Self::Output {
        let tuple: Tuple<f64> = self.into();
        tuple * rhs
    }
}

impl Mul<i32> for Tuple<f64> {
    type Output = Tuple<f64>;

    fn mul(self, rhs: i32) -> Self::Output {
        self * rhs as f64
    }
}

impl<T> Div<T> for Tuple<f64>
where
    T: Into<f64>,
{
    type Output = Tuple<f64>;

    fn div(self, rhs: T) -> Self::Output {
        let right: f64 = rhs.into();
        Tuple::new(
            self.x / right,
            self.y / right,
            self.z / right,
            self.w / right,
        )
    }
}

impl<T> Div<T> for Tuple<i32>
where
    T: Into<f64>,
{
    type Output = Tuple<f64>;
    fn div(self, rhs: T) -> Self::Output {
        let tuple: Tuple<f64> = self.into();
        tuple / rhs
    }
}

impl<T> Neg for Tuple<T>
where
    T: Num + Copy + Into<f64> + Neg<Output = T>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y, -self.z, -self.w)
    }
}

impl<T: Num + Copy + Into<f64>> Magnitude for Tuple<T> {
    fn magnitude(&self) -> f64 {
        let square_sum =
            (self.x * self.x + self.y * self.y + self.z * self.z + self.w * self.w).into();
        square_sum.sqrt()
    }
}

impl<T: Num + Copy + Into<f64>> Normalize for Tuple<T> {
    type Output = Tuple<f64>;
    fn normalize(&self) -> Self::Output {
        let magnitude = self.magnitude();
        let f64_tuple: Tuple<f64> =
            Tuple::new(self.x.into(), self.y.into(), self.z.into(), self.w.into());
        f64_tuple / magnitude
    }
}

impl<T: Num + Copy + Into<f64>> Dot for Tuple<T> {
    type Output = T;
    fn dot(&self, vec: &Tuple<T>) -> Self::Output {
        self.x * vec.x + self.y * vec.y + self.z * vec.z + self.w * vec.w
    }
}

impl<T: Num + Copy + Into<f64>> CrossProduct for Tuple<T> {
    type Output = Tuple<T>;
    fn cross(&self, vec: &Self) -> Self::Output {
        Self::Output::new(
            self.y * vec.z - self.z * vec.y,
            self.z * vec.x - self.x * vec.z,
            self.x * vec.y - self.y * vec.x,
            zero()
        )
    }
}

pub fn tuple<T: Num + Copy + Into<f64>>(x: T, y: T, z: T, w: T) -> Tuple<T> {
    Tuple::new(x, y, z, w)
}

pub fn point<T: Num + Copy + Into<f64>>(x: T, y: T, z: T) -> Tuple<T> {
    Tuple::new(x, y, z, one())
}

pub fn vector<T: Num + Copy + Into<f64>>(x: T, y: T, z: T) -> Tuple<T> {
    Tuple::new(x, y, z, zero())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tuple_with_w_eq_1_is_a_point() {
        let a = tuple(4.3, -4.2, 3.1, 1.0);

        assert_eq!(4.3, a.x);
        assert_eq!(-4.2, a.y);
        assert_eq!(3.1, a.z);
        assert_eq!(1.0, a.w);
        assert!(a.is_point());
        assert!(!a.is_vector());
    }

    #[test]
    fn test_tuple_with_w_eq_0_is_a_vector() {
        let a = tuple(4.3, -4.2, 3.1, 0.0);

        assert_eq!(4.3, a.x);
        assert_eq!(-4.2, a.y);
        assert_eq!(3.1, a.z);
        assert_eq!(0.0, a.w);
        assert!(!a.is_point());
        assert!(a.is_vector());
    }

    #[test]
    fn test_point_create_tuple_with_w_eq_1() {
        let a = point(4, -4, 3);
        assert_eq!(a, Tuple::new(4, -4, 3, 1));
        assert!(a.is_point());
    }
    #[test]
    fn test_vector_create_tuple_with_w_eq_0() {
        let a = vector(4, -4, 3);
        assert_eq!(a, Tuple::new(4, -4, 3, 0));
        assert!(a.is_vector());
    }

    #[test]
    fn test_adding_point_and_vector() {
        let p = point(3, -2, 5);
        let v = vector(-2, 3, 1);

        let sum = p + v;
        assert_eq!(sum, point(1, 1, 6));
    }

    #[test]
    fn test_adding_two_vectors() {
        let v1 = vector(3, -2, 5);
        let v2 = vector(-2, 3, 1);

        let sum = v1 + v2;
        assert_eq!(sum, vector(1, 1, 6));
    }

    #[test]
    fn test_adding_two_points_makes_no_sense() {
        let p1 = point(3, -2, 5);
        let p2 = point(-2, 3, 1);

        let meaningless_tuple = p1 + p2;
        assert!(!meaningless_tuple.is_point(), "not a point");
        assert!(!meaningless_tuple.is_vector(), "not a vector");
    }

    #[test]
    fn test_substracting_two_points() {
        let p1 = point(3, 2, 1);
        let p2 = point(5, 6, 7);

        assert_eq!(p1 - p2, vector(-2, -4, -6));
    }

    #[test]
    fn test_substracting_vector_from_a_point() {
        let p = point(3, 2, 1);
        let v = vector(5, 6, 7);

        assert_eq!(p - v, point(-2, -4, -6));
    }

    #[test]
    fn test_substracting_two_vectors() {
        let v1 = vector(3, 2, 1);
        let v2 = vector(5, 6, 7);

        assert_eq!(v1 - v2, vector(-2, -4, -6));
    }

    #[test]
    fn test_substracting_vector_from_zero_vector() {
        let zero = vector(0, 0, 0);
        let v = vector(1, 2, 3);

        assert_eq!(zero - v, vector(-1, -2, -3));
    }

    #[test]
    fn negating_a_tuple() {
        let a = tuple(1, 2, 3, -4);
        assert_eq!(-a, tuple(-1, -2, -3, 4))
    }

    #[test]
    fn multiplying_integer_tuple_by_integer_scalar() {
        let a = tuple(1, -2, 3, -4);
        assert_eq!(a * 2, tuple(2, -4, 6, -8));
    }

    #[test]
    fn multiplying_integer_tuple_by_float() {
        let a = tuple(1, -2, 3, -4);

        assert_eq!(a * 0.5, tuple(0.5, -1.0, 1.5, -2.0));
    }

    #[test]
    fn multiplying_float_tuple_by_integer() {
        let a = tuple(1.0, 2.0, 3.0, 4.0);

        assert_eq!(a * 2, tuple(2.0, 4.0, 6.0, 8.0))
    }

    #[test]
    fn multiplying_float_tuple_by_float() {
        let a = tuple(1.0, 2.0, 3.0, 4.0);

        assert_eq!(a * 2.0, tuple(2.0, 4.0, 6.0, 8.0))
    }

    #[test]
    fn dividing_integer_tuple_by_integer_scalar() {
        let a = tuple(2, -2, 3, 5);
        assert_eq!(a / 2, tuple(1.0, -1.0, 1.5, 2.5));
    }

    #[test]
    fn dividing_integer_tuple_by_float_scalar() {
        let a = tuple(2, -2, 3, 5);
        assert_eq!(a / 2.0, tuple(1.0, -1.0, 1.5, 2.5));
    }

    #[test]
    fn dividing_float_tuple_by_integer_scalar() {
        let a = tuple(2.0, -2.0, 3.0, 5.0);
        assert_eq!(a / 2, tuple(1.0, -1.0, 1.5, 2.5));
    }

    #[test]
    fn dividing_float_tuple_by_float_scalar() {
        let a = tuple(2.0, -2.0, 3.0, 5.0);
        assert_eq!(a / 2.0, tuple(1.0, -1.0, 1.5, 2.5));
    }

    #[test]
    fn magnitude_of_1_0_0_vector_is_1() {
        let vector = vector(1, 0, 0);
        assert_eq!(1.0, vector.magnitude());
    }

    #[test]
    fn magnitude_of_0_1_0_vector_is_1() {
        let vector = vector(0, 1, 0);
        assert_eq!(1.0, vector.magnitude());
    }

    #[test]
    fn magnitude_of_0_0_1_vector_is_1() {
        let vector = vector(0, 0, 1);
        assert_eq!(1.0, vector.magnitude());
    }

    #[test]
    fn computes_properly_magnitude_of_int_vectors() {
        let vector = vector(-1, -2, -3);
        assert_eq!(14f64.sqrt(), vector.magnitude());
    }

    #[test]
    fn computes_properly_magnitude_of_floart_vectors() {
        let vector = vector(10.0, 20.0, -30.0);
        assert_eq!(1400f64.sqrt(), vector.magnitude());
    }

    #[test]
    fn normalizing_vector_4_0_0_gives_1_0_0() {
        let v1 = vector(4, 0, 0);
        let unit: Tuple<f64> = vector(1.0, 0.0, 0.0);
        assert_eq!(unit, v1.normalize());
    }

    #[test]
    fn normalizing_vector_1_2_3() {
        let v1 = vector(1, 2, 3);

        let length = 1.0 / 14f64.sqrt();
        let expected = vector(length, length * 2.0, length * 3.0);
        assert_eq!(expected, v1.normalize());
    }

    #[test]
    fn magnitude_of_normalized_vector() {
        let v = vector(1, 2, 3);
        let norm = v.normalize();

        assert_eq!(1.0, norm.magnitude());
    }

    #[test]
    fn dot_product_of_two_vectors() {
        let a = vector(1, 2, 3);
        let b = vector(2, 3, 4);

        assert_eq!(20, a.dot(&b));
    }

    #[test]
    fn dot_product_of_two_f64_vectors() {
        let a = vector(10.0, 20.0, 30.0);
        let b = vector(1.0, 1.0, 1.0);

        assert_eq!(60.0, a.dot(&b));
    }

    #[test]
    fn cross_product_of_integer_vectors() {
        let a = vector(1,2,3);
        let b = vector(2,3,4);
        assert_eq!(vector(-1,2,-1), a.cross(&b));
        assert_eq!(vector(1,-2,1), b.cross(&a));
    }

    #[test]
    fn cross_product_of_float_vectors() {
        let a = vector(1.0, 2.0, 3.0);
        let b = vector(2.0, 3.0, 4.0);
        assert_eq!(vector(-1.0, 2.0, -1.0), a.cross(&b));
        assert_eq!(vector(1.0,-2.0,1.0), b.cross(&a));
    }
}
