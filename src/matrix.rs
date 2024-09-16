use approx::{AbsDiffEq, RelativeEq};
use num::{zero, ToPrimitive};

use crate::tuple::Tuple;
use std::ops::{Div, Index, Mul};

pub trait Identity {
    fn identity() -> Self;
}
#[derive(PartialEq, Clone, Debug)]
pub struct Matrix<const N: usize> {
    cells: [[f64; N]; N],
}

pub type Matrix2 = Matrix<2>;
pub type Matrix3 = Matrix<3>;
pub type Matrix4 = Matrix<4>;

pub trait Transpose {
    fn transpose(&self) -> Self;
}

pub trait Determinant {
    fn determinant(&self) -> f64;
}

pub trait Submatrix {
    type Output;
    fn submatrix(&self, row: usize, col: usize) -> Self::Output;
}

pub trait Inverse {
    fn is_invertible(&self) -> bool;
    fn inverse(&self) -> Self;
}

pub trait Minor {
    fn minor(&self, row: usize, col: usize) -> f64;
    fn cofactor(&self, row: usize, col: usize) -> f64;
}

impl<const N: usize> Matrix<N> {
    pub fn row(&self, index: usize) -> [f64; N] {
        self.cells[index]
    }

    pub fn col(&self, index: usize) -> [f64; N] {
        self.cells.map(|row| row[index])
    }

    pub fn cell(&self, row: usize, col: usize) -> f64 {
        self.cells[row][col]
    }

    pub fn round(&self, digits: u32) -> Self {
        let ten_power = 10_u32.pow(digits).to_f64().unwrap();
        let cells = self
            .cells
            .map(|row| row.map(|cell| (cell * ten_power).round() / ten_power));
        Matrix { cells }
    }
}

impl<const N: usize> Index<usize> for Matrix<N> {
    type Output = [f64; N];
    fn index(&self, index: usize) -> &Self::Output {
        &self.cells[index]
    }
}

impl Matrix<2> {
    pub fn new(m00: f64, m01: f64, m10: f64, m11: f64) -> Self {
        Self {
            cells: [[m00, m01], [m10, m11]],
        }
    }
}

impl<const N: usize> Identity for Matrix<N> {
    fn identity() -> Self {
        let rows = (0..N)
            .map(|i| {
                (0..N)
                    .map(|j| if i == j { 1.0 } else { 0.0 })
                    .collect::<Vec<f64>>()
                    .try_into()
                    .unwrap()
            })
            .collect::<Vec<[f64; N]>>()
            .try_into()
            .unwrap();
        Self { cells: rows }
    }
}

impl<const N: usize> AbsDiffEq for Matrix<N> {
    type Epsilon = f64;
    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.cells
            .into_iter()
            .zip(other.cells)
            .all(|(row, other_row)| row.abs_diff_eq(&other_row, epsilon))
    }
}

impl<const N: usize> RelativeEq for Matrix<N> {
    fn default_max_relative() -> Self::Epsilon {
        f64::default_max_relative()
    }
    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        self.cells
            .into_iter()
            .zip(other.cells)
            .all(|(row, other_row)| row.relative_eq(&other_row, epsilon, max_relative))
    }
}

impl Matrix<3> {
    pub fn new(
        m00: f64,
        m01: f64,
        m02: f64,
        m10: f64,
        m11: f64,
        m12: f64,
        m20: f64,
        m21: f64,
        m22: f64,
    ) -> Self {
        Self {
            cells: [[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]],
        }
    }
}

impl Matrix<4> {
    pub fn new(
        m00: f64,
        m01: f64,
        m02: f64,
        m03: f64,
        m10: f64,
        m11: f64,
        m12: f64,
        m13: f64,
        m20: f64,
        m21: f64,
        m22: f64,
        m23: f64,
        m30: f64,
        m31: f64,
        m32: f64,
        m33: f64,
    ) -> Self {
        Self {
            cells: [
                [m00, m01, m02, m03],
                [m10, m11, m12, m13],
                [m20, m21, m22, m23],
                [m30, m31, m32, m33],
            ],
        }
    }
}

impl<const N: usize> Mul<Matrix<N>> for Matrix<N> {
    type Output = Matrix<N>;

    fn mul(self, rhs: Matrix<N>) -> Self::Output {
        let mut cells = [[0.0; N]; N];
        for row_i in 0..N {
            for col_i in 0..N {
                let sum = (0..N)
                    .map(|i| self.cell(row_i, i) * rhs.cell(i, col_i))
                    .sum();
                cells[row_i][col_i] = sum
            }
        }
        Self::Output { cells }
    }
}

impl<const N: usize> Mul<f64> for Matrix<N> {
    type Output = Matrix<N>;

    fn mul(self, rhs: f64) -> Self::Output {
        let cells = self.cells.map(|row| row.map(|cell| cell * rhs));
        Matrix { cells }
    }
}

impl<const N: usize> Div<f64> for Matrix<N> {
    type Output = Matrix<N>;

    fn div(self, rhs: f64) -> Self::Output {
        let cells = self.cells.map(|row| row.map(|cell| cell / rhs));
        Matrix { cells }
    }
}

impl Mul<Tuple<f64>> for Matrix<4> {
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

impl<const N: usize> Transpose for Matrix<N> {
    fn transpose(&self) -> Self {
        let mut cells = [[0.0; N]; N];
        for i in 0..N {
            cells[i] = self.col(i)
        }
        Self { cells }
    }
}

impl Determinant for Matrix<2> {
    fn determinant(&self) -> f64 {
        self.cell(0, 0) * self.cell(1, 1) - self.cell(0, 1) * self.cell(1, 0)
    }
}

impl Determinant for Matrix<3> {
    fn determinant(&self) -> f64 {
        self.row(0)
            .into_iter()
            .enumerate()
            .map(|(index, cell)| cell * self.cofactor(0, index))
            .sum()
    }
}

impl Determinant for Matrix<4> {
    fn determinant(&self) -> f64 {
        self.row(0)
            .into_iter()
            .enumerate()
            .map(|(index, cell)| cell * self.cofactor(0, index))
            .sum()
    }
}

impl Submatrix for Matrix<3> {
    type Output = Matrix<2>;
    fn submatrix(&self, row: usize, col: usize) -> Self::Output {
        let mut rows: Vec<[f64; 2]> = Vec::new();
        for i in 0..3 {
            if i == row {
                continue; // skip row when it matches passed column
            }
            let mut cells: Vec<f64> = Vec::new();
            for j in 0..3 {
                if j == col {
                    continue; // skip cell when it matches both row and column
                }
                cells.push(self.cell(i, j))
            }
            let new_row = cells.try_into().unwrap();
            rows.push(new_row)
        }
        Matrix2 {
            cells: rows.try_into().unwrap(),
        }
    }
}

impl Submatrix for Matrix<4> {
    type Output = Matrix<3>;
    fn submatrix(&self, row: usize, col: usize) -> Self::Output {
        let mut rows: Vec<[f64; 3]> = Vec::new();
        for i in 0..4 {
            if i == row {
                continue; // skip row when it matches passed column
            }
            let mut cells: Vec<f64> = Vec::new();
            for j in 0..4 {
                if j == col {
                    continue; // skip cell when it matches both row and column
                }
                cells.push(self.cell(i, j));
            }
            rows.push(cells.try_into().unwrap())
        }
        Matrix3 {
            cells: rows.try_into().unwrap(),
        }
    }
}

impl Minor for Matrix<3> {
    fn minor(&self, row: usize, col: usize) -> f64 {
        self.submatrix(row, col).determinant()
    }

    fn cofactor(&self, row: usize, col: usize) -> f64 {
        let minor = self.minor(row, col);
        if (row + col) % 2 == 0 {
            minor
        } else {
            -minor
        }
    }
}

impl Minor for Matrix<4> {
    fn minor(&self, row: usize, col: usize) -> f64 {
        self.submatrix(row, col).determinant()
    }

    fn cofactor(&self, row: usize, col: usize) -> f64 {
        let minor = self.minor(row, col);
        if (row + col) % 2 == 0 {
            minor
        } else {
            -minor
        }
    }
}

impl Inverse for Matrix<4> {
    fn inverse(&self) -> Self {
        let mut cells = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                cells[i][j] = self.cofactor(i, j)
            }
        }
        let cofactors = Matrix { cells }.transpose();
        cofactors / self.determinant()
    }

    fn is_invertible(&self) -> bool {
        self.determinant() != zero()
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

    #[test]
    fn transposing_2x2_matrix() {
        let m = Matrix2::new(1.0, 2.0, 3.0, 4.0);
        let transposed = Matrix2::new(1.0, 3.0, 2.0, 4.0);

        assert_eq!(transposed, m.transpose());
    }

    #[test]
    fn transposing_4x4_matrix() {
        let m = Matrix4::new(
            0.0, 9.0, 3.0, 0.0, 9.0, 8.0, 3.0, 8.0, 1.0, 8.0, 5.0, 3.0, 0.0, 0.0, 5.0, 8.0,
        );
        let transposed = Matrix4::new(
            0.0, 9.0, 1.0, 0.0, 9.0, 8.0, 8.0, 0.0, 3.0, 3.0, 5.0, 5.0, 0.0, 8.0, 3.0, 8.0,
        );

        assert_eq!(transposed, m.transpose());
    }

    #[test]
    fn transposing_identity_matrix() {
        assert_eq!(Matrix2::identity(), Matrix2::identity().transpose());
        assert_eq!(Matrix3::identity(), Matrix3::identity().transpose());
        assert_eq!(Matrix4::identity(), Matrix4::identity().transpose());
    }

    #[test]
    fn determinant_of_2x2_matrix() {
        let m = Matrix2::new(1.0, 5.0, -3.0, 2.0);
        assert_eq!(17.0, m.determinant());
    }

    #[test]
    fn submatrix_of_3x3_matrix_is_2x2_matrix() {
        let m = Matrix3::new(1.0, 5.0, 0.0, -3.0, 2.0, 7.0, 0.0, 6.0, -3.0);
        let submatrix = Matrix2::new(-3.0, 2.0, 0.0, 6.0);

        assert_eq!(submatrix, m.submatrix(0, 2));
    }

    #[test]
    fn submatrix_of_4x4_matrix_is_3x3_matrix() {
        let m = Matrix4::new(
            -6.0, 1.0, 1.0, 6.0, -8.0, 5.0, 8.0, 6.0, -1.0, 0.0, 8.0, 2.0, -7.0, 1.0, -1.0, 1.0,
        );
        let submatrix = Matrix3::new(-6.0, 1.0, 6.0, -8.0, 8.0, 6.0, -7.0, -1.0, 1.0);
        assert_eq!(submatrix, m.submatrix(2, 1))
    }

    #[test]
    fn minor_of_3x3_submatrix() {
        let m = Matrix3::new(3.0, 5.0, 0.0, 2.0, -1.0, -7.0, 6.0, -1.0, 5.0);
        let b = m.submatrix(1, 0);
        assert_eq!(25.0, b.determinant());
        assert_eq!(25.0, m.minor(1, 0));
    }

    #[test]
    fn cofactor_of_3x3_matrix() {
        let m = Matrix3::new(3.0, 5.0, 0.0, 2.0, -1.0, -7.0, 6.0, -1.0, 5.0);
        assert_eq!(-12.0, m.cofactor(0, 0));
        assert_eq!(-25.0, m.cofactor(1, 0));
    }

    #[test]
    fn calculating_determinant_of_3x3_matrix() {
        let m = Matrix3::new(1.0, 2.0, 6.0, -5.0, 8.0, -4.0, 2.0, 6.0, 4.0);

        assert_eq!(56.0, m.cofactor(0, 0));
        assert_eq!(12.0, m.cofactor(0, 1));
        assert_eq!(-46.0, m.cofactor(0, 2));
        assert_eq!(-196.0, m.determinant());
    }

    #[test]
    fn calculating_determinant_of_4x4_matrix() {
        let m = Matrix4::new(
            -2.0, -8.0, 3.0, 5.0, -3.0, 1.0, 7.0, 3.0, 1.0, 2.0, -9.0, 6.0, -6.0, 7.0, 7.0, -9.0,
        );
        assert_eq!(690.0, m.cofactor(0, 0));
        assert_eq!(447.0, m.cofactor(0, 1));
        assert_eq!(210.0, m.cofactor(0, 2));
        assert_eq!(51.0, m.cofactor(0, 3));
        assert_eq!(-4071.0, m.determinant());
    }

    #[test]
    fn testing_invertible_matrix_for_invertibility() {
        let m = Matrix4::new(
            6.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 6.0, 4.0, -9.0, 3.0, -7.0, 9.0, 1.0, 7.0, -6.0,
        );
        assert_eq!(-2120.0, m.determinant());
        assert!(m.is_invertible());
    }

    #[test]
    fn testing_not_invertible_matrix_for_invertibility() {
        let m = Matrix4::new(
            -4.0, 2.0, -2.0, -3.0, 9.0, 6.0, 2.0, 6.0, 0.0, -5.0, 1.0, -5.0, 0.0, 0.0, 0.0, 0.0,
        );
        assert_eq!(0.0, m.determinant());
        assert!(!m.is_invertible());
    }

    #[test]
    fn calculating_inverse_of_matrix_4x4() {
        let m = Matrix4::new(
            -5.0, 1.0, 7.0, 1.0, 2.0, -5.0, 7.0, -3.0, 6.0, 1.0, -6.0, 7.0, -8.0, 8.0, -7.0, 4.0,
        );

        let expected = Matrix4::new(
            0.21805, -0.80827, -0.07895, -0.52256, 0.45113, -1.45677, -0.22368, -0.81391, 0.24060,
            -0.44361, -0.05263, -0.30075, -0.04511, 0.52068, 0.19737, 0.30639,
        );

        assert_relative_eq!(expected, m.inverse().round(5))
    }

    #[test]
    fn inverse_of_identity_matrix_is_identity() {
        let m = Matrix4::identity();

        assert_eq!(Matrix4::identity(), m.inverse())
    }

    #[test]
    fn inverse_of_transposed_matrix_is_same_as_traspose_of_inversed() {
        let m = Matrix4::new(
            -5.0, 1.0, 7.0, 1.0, 2.0, -5.0, 7.0, -3.0, 6.0, 1.0, -6.0, 7.0, -8.0, 8.0, -7.0, 4.0,
        );

        let expected = m.transpose().inverse();

        assert_eq!(expected, m.inverse().transpose())
    }
}
