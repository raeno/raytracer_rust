use crate::tuple::Tuple;
use std::ops::{Index, Mul};

pub trait Identity {
    fn identity() -> Self;
}

#[derive(PartialEq, Clone, Debug)]
pub struct Matrix2 {
    cells: [[f64; 2]; 2],
}

#[derive(PartialEq, Clone, Debug)]
pub struct Matrix3 {
    cells: [[f64; 3]; 3],
}

#[derive(PartialEq, Clone, Debug)]
pub struct Matrix4 {
    cells: [[f64; 4]; 4],
}

impl Matrix2 {
    pub fn new(m11: f64, m12: f64, m21: f64, m22: f64) -> Self {
        Self {
            cells: [[m11, m12], [m21, m22]],
        }
    }
}

impl Index<usize> for Matrix2 {
    type Output = [f64; 2];
    fn index(&self, index: usize) -> &Self::Output {
        &self.cells[index]
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
        &self.cells[index]
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
            cells: [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]],
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
            cells: [
                [m11, m12, m13, m14],
                [m21, m22, m23, m24],
                [m31, m32, m33, m34],
                [m41, m42, m43, m44],
            ],
        }
    }

    pub fn row(&self, index: usize) -> [f64; 4] {
        self.cells[index]
    }

    pub fn col(&self, index: usize) -> [f64; 4] {
        self.cells.map(|row| row[index])
    }

    pub fn cell(&self, row: usize, col: usize) -> f64 {
        self.cells[row][col]
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
        &self.cells[index]
    }
}

impl Mul<Matrix4> for Matrix4 {
    type Output = Matrix4;

    fn mul(self, rhs: Matrix4) -> Self::Output {
        let mut cells = [[0.0; 4]; 4];
        for row_i in 0..4 {
            for col_i in 0..4 {
                cells[row_i][col_i] = self.cell(row_i, 0) * rhs.cell(0, col_i)
                    + self.cell(row_i, 1) * rhs.cell(1, col_i)
                    + self.cell(row_i, 2) * rhs.cell(2, col_i)
                    + self.cell(row_i, 3) * rhs.cell(3, col_i)
            }
        }
        Self::Output { cells }
    }
}

impl Mul<Tuple<f64>> for Matrix4 {
    type Output = Tuple<f64>;
    fn mul(self, rhs: Tuple<f64>) -> Self::Output {
        let sums = (0..4)
            .map(|i| {
                self.cell(i, 0) * rhs.x
                    + self.cell(i, 1) * rhs.y
                    + self.cell(i, 2) * rhs.z
                    + self.cell(i, 3) * rhs.w
            })
            .collect::<Vec<f64>>();
        Self::Output::new(sums[0], sums[1], sums[2], sums[3])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tuple::tuple;

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

    #[test]
    fn test_multiplying_4x4_matrices() {
        let a = Matrix4::new(
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0,
        );

        let b = Matrix4::new(
            -2.0, 1.0, 2.0, 3.0, 3.0, 2.0, 1.0, -1.0, 4.0, 3.0, 6.0, 5.0, 1.0, 2.0, 7.0, 8.0,
        );

        let result = Matrix4::new(
            20.0, 22.0, 50.0, 48.0, 44.0, 54.0, 114.0, 108.0, 40.0, 58.0, 110.0, 102.0, 16.0, 26.0,
            46.0, 42.0,
        );
        assert_eq!(result, a * b);
    }

    #[test]
    fn test_multiplying_4x4_matrice_by_identity() {
        let a = Matrix4::new(
            0.0, 1.0, 2.0, 4.0, 1.0, 2.0, 4.0, 8.0, 2.0, 4.0, 8.0, 16.0, 4.0, 8.0, 16.0, 32.0,
        );

        assert_eq!(a.clone(), a * Matrix4::identity());
    }

    #[test]
    fn test_multiplying_identity_matrix_by_tuple() {
        let a = tuple(1.0, 2.0, 3.0, 4.0);
        assert_eq!(a.clone(), Matrix4::identity() * a);
    }
}
