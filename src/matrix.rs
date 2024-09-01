use std::ops::Index;

pub trait Identity {
    fn identity() -> Self;
}

#[derive(PartialEq, Clone, Debug)]
pub struct Matrix2 {
    rows: [[f64; 2]; 2],
}

#[derive(PartialEq, Clone, Debug)]
pub struct Matrix3 {
    rows: [[f64; 3]; 3],
}

#[derive(PartialEq, Clone, Debug)]
pub struct Matrix4 {
    rows: [[f64; 4]; 4],
}

impl Matrix2 {
    pub fn new(m11: f64, m12: f64, m21: f64, m22: f64) -> Self {
        Self {
            rows: [[m11, m12], [m21, m22]],
        }
    }
}

impl Index<usize> for Matrix2 {
    type Output = [f64; 2];
    fn index(&self, index: usize) -> &Self::Output {
        &self.rows[index]
    }
}

impl Identity for Matrix2 {
    fn identity() -> Self {
        Self::new(1.0, 0.0, 0.0, 1.0)
    }
}

impl Index<usize> for Matrix3 {
    type Output = [f64; 3];
    fn index(&self, index: usize) -> &Self::Output {
        &self.rows[index]
    }
}

impl Matrix3 {
    pub fn new(
        m11: f64,
        m12: f64,
        m13: f64,
        m21: f64,
        m22: f64,
        m23: f64,
        m31: f64,
        m32: f64,
        m33: f64,
    ) -> Self {
        Self {
            rows: [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]],
        }
    }
}

impl Identity for Matrix3 {
    fn identity() -> Self {
        Self::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    }
}

impl Matrix4 {
    pub fn new(
        m11: f64,
        m12: f64,
        m13: f64,
        m14: f64,
        m21: f64,
        m22: f64,
        m23: f64,
        m24: f64,
        m31: f64,
        m32: f64,
        m33: f64,
        m34: f64,
        m41: f64,
        m42: f64,
        m43: f64,
        m44: f64,
    ) -> Self {
        Self {
            rows: [
                [m11, m12, m13, m14],
                [m21, m22, m23, m24],
                [m31, m32, m33, m34],
                [m41, m42, m43, m44],
            ],
        }
    }
}

impl Identity for Matrix4 {
    fn identity() -> Self {
        Self::new(
            1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        )
    }
}

impl Index<usize> for Matrix4 {
    type Output = [f64; 4];
    fn index(&self, index: usize) -> &Self::Output {
        &self.rows[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constructing_2x2_matrix() {
        let m = Matrix2::new(1.0, 2.0, 3.0, 4.0);

        assert_eq!(1.0, m[0][0]);
        assert_eq!(2.0, m[0][1]);
        assert_eq!(3.0, m[1][0]);
        assert_eq!(4.0, m[1][1]);
    }

    #[test]
    fn test_constructing_3x3_matrix() {
        let m = Matrix3::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

        assert_eq!(1.0, m[0][0]);
        assert_eq!(2.0, m[0][1]);
        assert_eq!(3.0, m[0][2]);
        assert_eq!(4.0, m[1][0]);
        assert_eq!(5.0, m[1][1]);
        assert_eq!(6.0, m[1][2]);
        assert_eq!(7.0, m[2][0]);
        assert_eq!(8.0, m[2][1]);
        assert_eq!(9.0, m[2][2]);
    }

    #[test]
    fn test_constructing_4x4_matrix() {
        let m = Matrix4::new(
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
        );

        assert_eq!(1.0, m[0][0]);
        assert_eq!(2.0, m[0][1]);
        assert_eq!(3.0, m[0][2]);
        assert_eq!(4.0, m[0][3]);

        assert_eq!(5.0, m[1][0]);
        assert_eq!(6.0, m[1][1]);
        assert_eq!(7.0, m[1][2]);
        assert_eq!(8.0, m[1][3]);

        assert_eq!(9.0, m[2][0]);
        assert_eq!(10.0, m[2][1]);
        assert_eq!(11.0, m[2][2]);
        assert_eq!(12.0, m[2][3]);

        assert_eq!(13.0, m[3][0]);
        assert_eq!(14.0, m[3][1]);
        assert_eq!(15.0, m[3][2]);
        assert_eq!(16.0, m[3][3]);
    }

    #[test]
    fn test_comparing_2x2_matrices() {
        let a = Matrix2::new(1.0, 2.0, 3.0, 4.0);
        let b = Matrix2::new(1.0, 2.0, 3.0, 4.0);
        let c = Matrix2::identity();

        assert_eq!(a, b);
        assert_ne!(a, c);
        assert_ne!(b, c);
    }

    #[test]
    fn test_comparing_3x3_matrices() {
        let a = Matrix3::new(1.0, 2.0, 3.0, -1.0, -2.0, -3.0, 0.0, 0.0, 0.0);
        let b = Matrix3::new(1.0, 2.0, 3.0, -1.0, -2.0, -3.0, 0.0, 0.0, 0.0);
        let c = Matrix3::identity();

        assert_eq!(a, b);
        assert_ne!(a, c);
        assert_ne!(b, c);
    }

    #[test]
    fn test_comparing_4x4_matrices() {
        let a = Matrix4::new(
            1.0, 2.0, 3.0, 4.0, -1.0, -2.0, -3.0, -4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
        );

        let b = Matrix4::new(
            1.0, 2.0, 3.0, 4.0, -1.0, -2.0, -3.0, -4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
        );
        let c = Matrix4::identity();

        assert_eq!(a, b);
        assert_ne!(a, c);
        assert_ne!(b, c);
    }
}
