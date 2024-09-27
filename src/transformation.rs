use crate::matrix::{Identity, Matrix4};

fn translation(x: f64, y: f64, z: f64) -> Matrix4 {
    let mut m = Matrix4::identity();
    m[0][3] = x;
    m[1][3] = y;
    m[2][3] = z;
    m
}

fn scaling(x: f64, y: f64, z: f64) -> Matrix4 {
    let mut m = Matrix4::identity();
    m[0][0] = x;
    m[1][1] = y;
    m[2][2] = z;
    m
}

fn rotation_x(r: f64) -> Matrix4 {
    let mut m = Matrix4::identity();
    m[1][1] = r.cos();
    m[1][2] = -r.sin();
    m[2][1] = r.sin();
    m[2][2] = r.cos();
    m
}

fn rotation_y(r: f64) -> Matrix4 {
    let mut m = Matrix4::identity();
    m[0][0] = r.cos();
    m[0][2] = r.sin();
    m[2][0] = -r.sin();
    m[2][2] = r.cos();
    m
}

fn rotation_z(r: f64) -> Matrix4 {
    let mut m = Matrix4::identity();
    m[0][0] = r.cos();
    m[0][1] = -r.sin();
    m[1][0] = r.sin();
    m[1][1] = r.cos();
    m
}

fn shearing(xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) -> Matrix4 {
    let mut m = Matrix4::identity();
    m[0][1] = xy;
    m[0][2] = xz;
    m[1][0] = yx;
    m[1][2] = yz;
    m[2][0] = zx;
    m[2][1] = zy;
    m
}

pub trait Transformation {
    fn translate(&self, x: f64, y: f64, z: f64) -> Self;
    fn scale(&self, x: f64, y: f64, z: f64) -> Self;
    fn rotate_x(&self, r: f64) -> Self;
    fn rotate_y(&self, r: f64) -> Self;
    fn rotate_z(&self, r: f64) -> Self;
    fn skew(&self, xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) -> Self;
}

impl Transformation for Matrix4 {
    fn translate(&self, x: f64, y: f64, z: f64) -> Self {
        self.clone() * translation(x, y, z)
    }

    fn scale(&self, x: f64, y: f64, z: f64) -> Self {
        self.clone() * scaling(x, y, z)
    }

    fn rotate_x(&self, r: f64) -> Self {
        self.clone() * rotation_x(r)
    }

    fn rotate_y(&self, r: f64) -> Self {
        self.clone() * rotation_y(r)
    }

    fn rotate_z(&self, r: f64) -> Self {
        self.clone() * rotation_z(r)
    }

    fn skew(&self, xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) -> Self {
        self.clone() * shearing(xy, xz, yx, yz, zx, zy)
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, FRAC_PI_4};

    use super::*;
    use crate::matrix::Inverse;
    use crate::tuple::{point, vector, Round, Tuple};

    #[test]
    fn multiplying_point_by_translation_matrix() {
        let transform = translation(5.0, -3.0, 2.0);
        let p: Tuple<f64> = point(-3, 4, 5).into();
        assert_ulps_eq!(point(2.0, 1.0, 7.0), transform * p)
    }

    #[test]
    fn multiplying_point_by_inverse_of_translation_matrix() {
        let transform = translation(5.0, -3.0, 2.0);
        let p: Tuple<f64> = point(-3, 4, 5).into();

        assert_ulps_eq!(point(-8.0, 7.0, 3.0), transform.inverse() * p)
    }

    #[test]
    fn translation_does_not_affect_vectors() {
        let transform = translation(5.0, -3.0, 2.0);
        let v = vector(-3.0, 4.0, 5.0);
        assert_ulps_eq!(v, transform * v)
    }

    #[test]
    fn scaling_matrix_applied_to_point() {
        let transform = scaling(2.0, 3.0, 4.0);
        let p: Tuple<f64> = point(-4, 6, 8).into();
        assert_ulps_eq!(point(-8.0, 18.0, 32.0), transform * p)
    }

    #[test]
    fn scaling_matrix_applied_to_vector() {
        let transform = scaling(2.0, 3.0, 4.0);
        let p: Tuple<f64> = vector(-4, 6, 8).into();
        assert_ulps_eq!(vector(-8.0, 18.0, 32.0), transform * p)
    }

    #[test]
    fn multiplying_vector_by_inverse_of_scaling_matrix() {
        let transform = scaling(2.0, 3.0, 4.0);
        let p: Tuple<f64> = vector(-4, 6, 8).into();
        assert_ulps_eq!(vector(-2.0, 2.0, 2.0), transform.inverse() * p)
    }

    #[test]
    fn reflection_is_scaling_by_negative_value() {
        let transform = scaling(-1.0, 1.0, 1.0);
        let p = point(2.0, 3.0, 4.0);
        assert_ulps_eq!(point(-2.0, 3.0, 4.0), transform * p)
    }

    #[test]
    fn rotating_a_point_around_x_axis() {
        let p = point(0.0, 1.0, 0.0);
        let half_quarter = rotation_x(FRAC_PI_4);
        let full_quarter = rotation_x(FRAC_PI_2);

        let rotated_half_quarter = point(0.0, FRAC_1_SQRT_2, FRAC_1_SQRT_2);
        assert_ulps_eq!(rotated_half_quarter, half_quarter * p);
        assert_ulps_eq!(point(0_f64, 0_f64, 1_f64), full_quarter * p)
    }

    #[test]
    fn inverse_of_x_rotation_rotates_in_opposite_direction() {
        let p = point(0.0, 1.0, 0.0);
        let half_quarter = rotation_x(FRAC_PI_4).inverse();

        let rotated_half_quarter = point(0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2);
        assert_ulps_eq!(rotated_half_quarter, half_quarter * p);
    }

    #[test]
    fn rotating_a_point_around_y_axis() {
        let p = point(0.0, 0.0, 1.0);
        let half_quarter = rotation_y(FRAC_PI_4);
        let full_quarter = rotation_y(FRAC_PI_2);

        let rotated_half_quarter = point(FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2);
        assert_ulps_eq!(rotated_half_quarter, half_quarter * p);
        assert_ulps_eq!(point(1_f64, 0_f64, 0_f64), full_quarter * p)
    }

    #[test]
    fn rotating_a_point_around_z_axis() {
        let p = point(0.0, 1.0, 0.0);

        let half_quarter = rotation_z(FRAC_PI_4);
        let full_quarter = rotation_z(FRAC_PI_2);

        let rotated_half_quarter = point(-FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0);
        assert_ulps_eq!(rotated_half_quarter, half_quarter * p);
        assert_ulps_eq!(point(-1_f64, 0_f64, 0_f64), full_quarter * p)
    }

    #[test]
    fn shearing_transformation_moves_x_in_proportion_to_y() {
        let transform = shearing(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let p = point(2.0, 3.0, 4.0);
        assert_eq!(point(5.0, 3.0, 4.0), transform * p)
    }

    #[test]
    fn shearing_transformation_moves_x_in_proportion_to_z() {
        let transform = shearing(0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
        let p = point(2.0, 3.0, 4.0);
        let expected = point(6.0, 3.0, 4.0);
        assert_eq!(expected, transform * p);
    }

    #[test]
    fn shearing_transformation_moves_y_in_proportion_to_x() {
        let transform = shearing(0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
        let p = point(2.0, 3.0, 4.0);
        let expected = point(2.0, 5.0, 4.0);
        assert_eq!(expected, transform * p);
    }

    #[test]
    fn shearing_transformation_moves_y_in_proportion_to_z() {
        let transform = shearing(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        let p = point(2.0, 3.0, 4.0);
        let expected = point(2.0, 7.0, 4.0);
        assert_eq!(expected, transform * p);
    }

    #[test]
    fn shearing_transformation_moves_z_in_proportion_to_x() {
        let transform = shearing(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        let p = point(2.0, 3.0, 4.0);
        let expected = point(2.0, 3.0, 6.0);
        assert_eq!(expected, transform * p);
    }

    #[test]
    fn shearing_transformation_moves_z_in_proportion_to_y() {
        let transform = shearing(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let p = point(2.0, 3.0, 4.0);
        let expected = point(2.0, 3.0, 7.0);
        assert_eq!(expected, transform * p);
    }

    #[test]
    fn individual_transformations_applied_in_sequence() {
        let p = point(1_f64, 0_f64, 1_f64);
        let a = rotation_x(FRAC_PI_2);
        let b = scaling(5.0, 5.0, 5.0);
        let c = translation(10.0, 5.0, 7.0);

        let p2 = a * p;
        assert_ulps_eq!(point(1.0, -1.0, 0.0), p2);
        let p3 = b * p2;
        assert_ulps_eq!(point(5.0, -5.0, 0.0), p3.round(10));

        let p4 = c * p3;
        assert_ulps_eq!(point(15.0, 0.0, 7.0), p4);
    }

    #[test]
    fn chained_transformations_must_be_applied_in_reverse_order() {
        let p = point(1_f64, 0_f64, 1_f64);
        let a = rotation_x(FRAC_PI_2);
        let b = scaling(5.0, 5.0, 5.0);
        let c = translation(10.0, 5.0, 7.0);
        let transform = c * b * a;

        let p4 = transform * p;
        assert_eq!(point(15.0, 0.0, 7.0), p4.round(10));
    }
}
