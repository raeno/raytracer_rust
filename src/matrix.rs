use crate::tuple::Tuple;
use std::ops::{Index, Mul};

pub trait Identity {
    fn identity() -> Self;
}
#[derive(PartialEq, Clone, Debug)]
pub struct Matrix<const N: usize> {
    cells: [[f64; N]; N]
}

pub type Matrix2 = Matrix<2>;
pub type Matrix3 = Matrix<3>;
pub type Matrix4 = Matrix<4>;

pub trait Transpose {
    fn transpose(&self) -> Self;
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
        let rows = (0..N).map( |i| 
            (0..N).map( |j| if i == j { 1.0 } else { 0.0 } ).collect::<Vec<f64>>().try_into().unwrap()
        ).collect::<Vec<[f64; N]>>().try_into().unwrap();
        Self { cells: rows }
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
                let sum = (0..N).map( |i| self.cell(row_i, i) * rhs.cell(i, col_i) ).sum();
                cells[row_i][col_i] = sum
            }
        }
        Self::Output { cells }
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
            0.0, 9.0, 3.0, 0.0,
            9.0, 8.0, 3.0, 8.0,
            1.0, 8.0, 5.0, 3.0,
            0.0, 0.0, 5.0, 8.0
        );
        let transposed = Matrix4::new(
            0.0, 9.0, 1.0, 0.0,
            9.0, 8.0, 8.0, 0.0,
            3.0, 3.0, 5.0, 5.0,
            0.0, 8.0, 3.0, 8.0
        );

        assert_eq!(transposed, m.transpose());
    }

    #[test]
    fn transposing_identity_matrix() {
        assert_eq!(Matrix2::identity(), Matrix2::identity().transpose());
        assert_eq!(Matrix3::identity(), Matrix3::identity().transpose());
        assert_eq!(Matrix4::identity(), Matrix4::identity().transpose());
    }
}
