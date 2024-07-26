use std::ops::{Add, Mul, Neg, Sub};

#[derive(PartialEq, Debug)]
pub struct Tuple {
    x: f64,
    y: f64,
    z: f64,
    w: f64,
}

impl Tuple {
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Self { x, y, z, w }
    }

    pub fn is_point(&self) -> bool {
        self.w == 1.0
    }

    pub fn is_vector(&self) -> bool {
        self.w == 0.0
    }
}

impl Add<Tuple> for Tuple {
    type Output = Tuple;

    fn add(self, rhs: Tuple) -> Self::Output {
        Tuple::new(
            self.x + rhs.x,
            self.y + rhs.y,
            self.z + rhs.z,
            self.w + rhs.w,
        )
    }
}

impl Sub<Tuple> for Tuple {
    type Output = Self;

    fn sub(self, rhs: Tuple) -> Self::Output {
        Tuple::new(
            self.x - rhs.x,
            self.y - rhs.y,
            self.z - rhs.z,
            self.w - rhs.w,
        )
    }
}

impl Neg for Tuple {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Tuple::new(-self.x, -self.y, -self.z, -self.w)
    }
}



// TODO: doesnt work, need to understand how to do multiplication both for int and float together
impl<T> Mul<T> for Tuple where T: Mul<Output=T> + Copy
{
    type Output = Tuple;
    fn mul(self, rhs: T) -> Self::Output {
        Tuple::new(
            self.x * rhs,
            self.y * rhs,
            self.z * rhs,
            self.w * rhs
        )
    }
}

pub fn tuple<T: Into<f64>>(x: T, y: T, z: T, w: T) -> Tuple {
    Tuple::new(x.into(), y.into(), z.into(), w.into())
}

pub fn point<T: Into<f64>>(x: T, y: T, z: T) -> Tuple {
    Tuple::new(x.into(), y.into(), z.into(), 1.0)
}

pub fn vector<T: Into<f64>>(x: T, y: T, z: T) -> Tuple {
    Tuple::new(x.into(), y.into(), z.into(), 0.0)
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
        assert_eq!(a, Tuple::new(4.0, -4.0, 3.0, 1.0));
        assert!(a.is_point());
    }
    #[test]
    fn test_vector_create_tuple_with_w_eq_0() {
        let a = vector(4, -4, 3);
        assert_eq!(a, Tuple::new(4.0, -4.0, 3.0, 0.0));
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
    fn multiplying_a_tuple_by_scalar() {
        let a = tuple(1, -2, 3, -4);
        let result = a * 3.5;

        assert_eq!(result, tuple(3.5, -7.0, 10.5, -14.0));
    }
}
