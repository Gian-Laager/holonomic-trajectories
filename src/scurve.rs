use nalgebra::{ComplexField, SVector, Scalar};
#[cfg(test)]
use quickcheck::Arbitrary;
#[cfg(test)]
use rand::prelude::*;
#[cfg(test)]
use rand::rngs::SmallRng;

#[derive(Clone, Debug)]
pub struct SCurve {
    start_time: f64,
    start_pos: f64,
    start_vel: f64,
    end_pos: f64,
    end_vel: f64,
    max_vel: f64,
    max_accel: f64,
}

#[cfg(test)]
impl Arbitrary for SCurve {
    fn arbitrary(g: &mut quickcheck::Gen) -> Self {
        let mut rng1 = SmallRng::seed_from_u64(u64::arbitrary(g));
        let mut rng2 = SmallRng::seed_from_u64(u64::arbitrary(g));
        let mut rng3 = SmallRng::seed_from_u64(u64::arbitrary(g));
        let mut arb_float = || rng1.gen_range(-10000.0..10000.0_f64);

        let mut arb_pos_float =
            || -> f64 { rng2.gen_range(f32::EPSILON.sqrt() as f64..10000.0_f64) };

        let max_vel = arb_pos_float();

        let start_time = arb_float();
        let start_pos = arb_float();
        let end_pos = arb_float();
        let start_vel = rng3.gen_range(-max_vel..max_vel);
        let max_accel = arb_pos_float();

        let distance = (end_pos - start_pos).abs();
        let t = (-start_vel + (start_vel.powi(2) + 2.0 * max_accel * distance).sqrt()) / max_accel;
        let max_end_vel = (start_vel + max_accel * t).abs() * 0.99;

        let end_vel = rng3.gen_range(-max_end_vel..max_end_vel);
        let maybe_curve = SCurve::try_new(
            start_time, start_pos, start_vel, end_pos, end_vel, max_vel, max_accel,
        );

        if let Ok(curve) = maybe_curve {
            return curve;
        } else {
            return SCurve::try_new(0.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5).unwrap();
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum SCurveGenerationError {
    StartVelTooHigh,
    EndVelTooHigh,
    NotEnoughTimeToReachEndVel,
    MaxVelMustBePositive,
    MaxAccelMustBePositive,
}

impl SCurve {
    pub fn try_new(
        start_time: f64,
        start_pos: f64,
        start_vel: f64,
        end_pos: f64,
        end_vel: f64,
        max_vel: f64,
        max_accel: f64,
    ) -> Result<Self, SCurveGenerationError> {
        let mut curve = SCurve {
            start_time,
            start_pos,
            start_vel,
            end_pos,
            end_vel,
            max_vel,
            max_accel,
        };

        let accel_distance = curve.accel_distance();
        let deccel_distance = curve.deccel_distance();
        let distance = curve.distance();

        if distance < accel_distance + deccel_distance {
            let max_vel_new = ((2.0 * curve.max_accel * distance
                + curve.start_vel.powi(2)
                + curve.end_vel.powi(2))
                / 2.0_f64)
                .sqrt()
                * 0.99;
            println!("Updated max_vel from {} to {}", curve.max_vel, max_vel_new);
            curve.max_vel = max_vel_new;
        }

        let time_for_d_with_full_accel =
            (-start_vel + (start_vel.powi(2) + 2.0 * max_accel * distance).sqrt()) / max_accel;
        if start_vel + max_accel * time_for_d_with_full_accel < end_vel {
            return Err(SCurveGenerationError::NotEnoughTimeToReachEndVel);
        }

        if curve.max_vel <= 0.0 {
            return Err(SCurveGenerationError::MaxVelMustBePositive);
        }

        if curve.max_accel <= 0.0 {
            return Err(SCurveGenerationError::MaxAccelMustBePositive);
        }

        if curve.start_vel.abs() > curve.max_vel {
            return Err(SCurveGenerationError::StartVelTooHigh);
        }

        if curve.end_vel.abs() > curve.max_vel {
            return Err(SCurveGenerationError::EndVelTooHigh);
        }

        println!("accel_duration: {}", curve.accel_duration());
        println!("cruise_duration: {}", curve.cruise_duration());
        println!("deccel_duration: {}", curve.deccel_duration());
        println!("distance: {}", curve.distance());
        println!("accel_distance: {}", curve.accel_distance());
        println!("deccel_distance: {}", curve.deccel_distance());
        println!("duration: {}", curve.duration());
        println!("accel_target: {}", curve.accel_target());
        println!("deccel_target: {}", curve.deccel_target());

        return Ok(curve);
    }

    fn accel_target(&self) -> f64 {
        if self.start_vel < self.vel_target() {
            return self.max_accel;
        } else {
            return -self.max_accel;
        }
    }

    fn deccel_target(&self) -> f64 {
        if self.end_vel > self.vel_target() {
            return self.max_accel;
        } else {
            return -self.max_accel;
        }
    }

    fn accel_duration(&self) -> f64 {
        (self.vel_target().abs() - self.start_vel.abs()) / self.max_accel
    }

    fn deccel_duration(&self) -> f64 {
        (self.vel_target().abs() - self.end_vel.abs()) / self.max_accel
    }

    fn deccel_distance(&self) -> f64 {
        let dt = self.deccel_duration();
        (self.vel_target() * dt + 0.5 * self.deccel_target() * dt.powi(2)).abs()
    }

    fn pos_after_cruise(&self) -> f64 {
        let deccel_direction = (self.end_pos - self.start_pos).signum();
        self.end_pos - self.deccel_distance() * deccel_direction
    }

    fn cruise_distance(&self) -> f64 {
        (self.pos_after_cruise() - self.pos_after_accel()).abs()
    }

    fn cruise_duration(&self) -> f64 {
        (self.cruise_distance() / self.max_vel).abs()
    }

    fn pos_after_accel(&self) -> f64 {
        let dt = self.accel_duration();
        self.start_pos + self.start_vel * dt + 0.5 * self.accel_target() * dt.powi(2)
    }

    fn accel_distance(&self) -> f64 {
        (self.pos_after_accel() - self.start_pos).abs() 
    }

    fn vel_target(&self) -> f64 {
        if self.end_pos > self.start_pos {
            return self.max_vel;
        } else {
            return -self.max_vel;
        }
    }

    pub fn distance(&self) -> f64 {
        (self.end_pos - self.start_pos).abs()
    }

    pub fn duration(&self) -> f64 {
        self.accel_duration() + self.cruise_duration() + self.deccel_duration()
    }

    pub fn sample_accel(&self, time: f64) -> f64 {
        let t = time - self.start_time;

        if t < 0.0 {
            return 0.0;
        } else if t < self.accel_duration() {
            return self.accel_target();
        } else if t < self.accel_duration() + self.cruise_duration() {
            return 0.0;
        } else if t < self.accel_duration() + self.cruise_duration() + self.deccel_duration() {
            return self.deccel_target();
        } else {
            return 0.0;
        }
    }

    pub fn sample_vel(&self, time: f64) -> f64 {
        let t = time - self.start_time;

        let vel = if t < 0.0 {
            self.start_vel
        } else if t < self.accel_duration() {
            self.start_vel + self.accel_target() * t
        } else if t < self.accel_duration() + self.cruise_duration() {
            self.vel_target()
        } else if t < self.accel_duration() + self.cruise_duration() + self.deccel_duration() {
            let t_prime = t - self.accel_duration() - self.cruise_duration();
            self.vel_target() + self.deccel_target() * t_prime
        } else {
            self.end_vel
        };

        return vel.min(self.max_vel).max(-self.max_vel);
    }

    pub fn sample_pos(&self, time: f64) -> f64 {
        let t = time - self.start_time;

        if t <= 0.0 {
            return self.start_pos;
        } else if t <= self.accel_duration() {
            return self.start_pos + self.start_vel * t + 0.5 * self.accel_target() * t.powi(2);
        } else if t <= self.accel_duration() + self.cruise_duration() {
            let t_prime = t - self.accel_duration();
            return self.pos_after_accel() + self.vel_target() * t_prime;
        } else if t < self.accel_duration() + self.cruise_duration() + self.deccel_duration() {
            let t_prime = t - self.accel_duration() - self.cruise_duration();
            return self.pos_after_cruise()
                + self.vel_target() * t_prime
                + 0.5 * self.deccel_target() * t_prime.powi(2);
        } else {
            return self.end_pos;
        }
    }
}

#[cfg(test)]
mod tests {
    use quickcheck_macros::quickcheck;

    use super::*;
    use crate::utils::test_utils::{assert_feq, float_compare, write_data_to_file};

    #[test]
    fn test_scurve() {
        let start_pos = 0.0;
        let start_vel = 0.0;
        let end_pos = 30.0;
        let end_vel = 0.0;
        let max_speed = 10.0;
        let max_accel = 5.0;

        const N_SAMPLES: usize = 100;
        const N_SAMPLES_F: f64 = N_SAMPLES as f64;

        let curve = SCurve::try_new(
            0.0, start_pos, start_vel, end_pos, end_vel, max_speed, max_accel,
        )
        .unwrap();

        let duration = curve.duration();
        let sample_pos = (0..=N_SAMPLES)
            .map(|i| curve.sample_pos(i as f64 * duration / N_SAMPLES_F))
            .collect::<Vec<f64>>();

        let sample_vel = (0..=N_SAMPLES)
            .map(|i| curve.sample_vel(i as f64 * duration / N_SAMPLES_F))
            .collect::<Vec<f64>>();

        let sample_accel = (0..=N_SAMPLES)
            .map(|i| curve.sample_accel(i as f64 * duration / N_SAMPLES_F))
            .collect::<Vec<f64>>();

        sample_accel.iter().for_each(|a| assert!(*a <= max_accel));
        sample_vel.iter().for_each(|v| assert!(*v <= max_speed));

        assert_feq(sample_pos[0], start_pos);
        assert_feq(*sample_pos.last().unwrap(), end_pos);
    }

    #[test]
    fn test_scurve_to_short() {
        let start_pos = 0.0;
        let start_vel = 0.0;
        let end_pos = 1.0;
        let end_vel = 0.0;
        let max_speed = 10.0;
        let max_accel = 5.0;

        const N_SAMPLES: usize = 100;
        const N_SAMPLES_F: f64 = N_SAMPLES as f64;

        let curve = SCurve::try_new(
            0.0, start_pos, start_vel, end_pos, end_vel, max_speed, max_accel,
        )
        .unwrap();

        let duration = curve.duration();
        let sample_pos = (0..=N_SAMPLES)
            .map(|i| curve.sample_pos(i as f64 * duration / N_SAMPLES_F))
            .collect::<Vec<f64>>();

        let pos_deriv = (0..N_SAMPLES)
            .map(|i| {
                let p1 = curve.sample_pos(i as f64 * duration / N_SAMPLES_F);
                let p2 = curve.sample_pos((i + 1) as f64 * duration / N_SAMPLES_F);
                (p2 - p1) / (duration / N_SAMPLES_F)
            })
            .collect::<Vec<f64>>();

        let pos_deriv_deriv = (0..N_SAMPLES)
            .map(|i| {
                let p1 = curve.sample_pos(i as f64 * duration / N_SAMPLES_F);
                let p2 = curve.sample_pos((i + 1) as f64 * duration / N_SAMPLES_F);
                let p3 = curve.sample_pos((i + 2) as f64 * duration / N_SAMPLES_F);
                (p3 - 2.0 * p2 + p1) / (duration / N_SAMPLES_F).powi(2)
            })
            .collect::<Vec<f64>>();

        let sample_vel = (0..=N_SAMPLES)
            .map(|i| curve.sample_vel(i as f64 * duration / N_SAMPLES_F))
            .collect::<Vec<f64>>();

        let vel_deriv = (0..N_SAMPLES)
            .map(|i| {
                let v1 = curve.sample_vel(i as f64 * duration / N_SAMPLES_F);
                let v2 = curve.sample_vel((i + 1) as f64 * duration / N_SAMPLES_F);
                (v2 - v1) / (duration / N_SAMPLES_F)
            })
            .collect::<Vec<f64>>();

        let sample_accel = (0..=N_SAMPLES)
            .map(|i| curve.sample_accel(i as f64 * duration / N_SAMPLES_F))
            .collect::<Vec<f64>>();

        sample_accel.iter().for_each(|a| assert!(*a <= max_accel));
        sample_vel.iter().for_each(|v| assert!(*v <= max_speed));

        pos_deriv
            .iter()
            .for_each(|v| assert!(v.abs() <= max_speed * 1.01));
        pos_deriv_deriv
            .iter()
            .for_each(|v| assert!(v.abs() <= max_accel * 1.01));
        vel_deriv
            .iter()
            .for_each(|v| assert!(v.abs() <= max_accel * 1.01));


        assert_feq(sample_pos[0], start_pos);
        assert_feq(*sample_pos.last().unwrap(), end_pos);
    }

    // #[ignore]
    #[quickcheck]
    fn scurve_quickcheck(mut curve: SCurve) -> bool {
        curve.start_time = 0.0;
        let curve = curve;
        let duration = curve.duration();
        let sample_pos = (0..=100)
            .map(|i| curve.sample_pos(i as f64 * duration / 100.0))
            .collect::<Vec<f64>>();

        let sample_vel = (0..=100)
            .map(|i| curve.sample_vel(i as f64 * duration / 100.0))
            .collect::<Vec<f64>>();

        let sample_accel = (0..=100)
            .map(|i| curve.sample_accel(i as f64 * duration / 100.0))
            .collect::<Vec<f64>>();

        // if !sample_accel.iter().all(|a| *a <= curve.max_accel * 1.1) {
        //     println!("Failed accel check");
        //     return false;
        // }
        //
        // if !sample_vel.iter().all(|v| *v <= curve.max_vel * 1.1) {
        //     println!("Failed vel check");
        //     return false;
        // }

        if !float_compare(sample_pos[0], curve.start_pos, 1e-2) {
            println!(
                "Failed start pos check, got: {}, expected: {}",
                sample_pos[0], curve.start_pos
            );
            write_data_to_file(&sample_pos, "scurve_pos.dat");
            write_data_to_file(&sample_vel, "scurve_vel.dat");
            write_data_to_file(&sample_accel, "scurve_acc.dat");
            return false;
        }

        if !float_compare(*sample_pos.last().unwrap(), curve.end_pos, 1e-2) {
            println!(
                "Failed end pos check, got: {}, expected: {}",
                *sample_pos.last().unwrap(),
                curve.end_pos
            );
            write_data_to_file(&sample_pos, "scurve_pos.dat");
            write_data_to_file(&sample_vel, "scurve_vel.dat");
            write_data_to_file(&sample_accel, "scurve_acc.dat");
            return false;
        }

        println!("################################################################################");

        return true;
    }
}
