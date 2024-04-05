#[cfg(test)]
extern crate quickcheck;

#[cfg(test)]
extern crate quickcheck_macros;

mod utils;
mod scurve;
mod spline;
mod polynomials;

use nalgebra::{DimName, SVector, Scalar, Vector2, Vector3};


struct TrajectorySample {
    time: f64,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
}

struct DriveSpecs {
    max_speed: f64,
    max_angular_vel: f64,

    max_acceleration: f64,
    max_angular_accel: f64,
}


// fn main() {
//     println!("Hello, world!");
// }
