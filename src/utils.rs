use nalgebra::{allocator::Allocator, DefaultAllocator, Vector, Vector2, Vector3, U1};

pub type Vec3 = Vector3<f64>;
pub type Vec2 = Vector2<f64>;

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use std::io::Write;

    pub fn assert_feq(a: f64, b: f64) {
        if !float_compare(a, b, 1e-4) {
            assert_eq!(a, b);
        }
    }

    pub fn float_compare(expect: f64, actual: f64, epsilon: f64) -> bool {
        return (expect - actual).abs() < epsilon;
        // let average = (expect + actual) / 2.0;
        //
        // if average < f64::EPSILON.sqrt() {
        //     return expect == actual;
        // }
        //
        // return (expect - actual) / average < epsilon;
    }

    pub fn write_data_to_file(data: &Vec<f64>, filename: &str) {
        let mut file = std::fs::File::create(filename).unwrap();
        for point in data.iter().enumerate() {
            writeln!(file, "{} {}", point.0, point.1).unwrap();
        }
    }
}
