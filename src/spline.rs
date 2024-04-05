use std::usize;

use nalgebra::SMatrix;
use nalgebra::SVector;

use crate::polynomials::*;
use crate::utils::*;

pub trait C2Spline<const D: usize> {
    fn evaluate(&self, x: f64) -> SVector<f64, D>;
    fn derivative(&self, x: f64) -> SVector<f64, D>;
    fn derivative2nd(&self, x: f64) -> SVector<f64, D>;
    fn range(&self) -> Interval;
}

#[derive(Debug, Clone)]
struct QuarticSpline<const D: usize> {
    curve: Polynomial<D, 5>,
}

impl<const D: usize> QuarticSpline<D> {
    fn new(
        start: SVector<f64, D>,
        end: SVector<f64, D>,
        start_vel: SVector<f64, D>,
        end_vel: SVector<f64, D>,
        start_acc: SVector<f64, D>,
    ) -> QuarticSpline<D> {
        let b = SMatrix::<f64, D, 5>::from_columns(
            [start, end, start_vel, end_vel, start_acc].as_ref(),
        )
        .transpose();

        let a_inv = SMatrix::<f64, 5, 5>::from_row_slice(
            [
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -4.0,
                4.0, -3.0, -1.0, -1.0, 3.0, -3.0, 2.0, 1.0, 0.5,
            ]
            .as_ref(),
        );

        let curve = a_inv * b;

        return QuarticSpline {
            curve: curve.transpose(),
        };
    }
}

impl<const D: usize> C2Spline<D> for QuarticSpline<D> {
    fn evaluate(&self, x: f64) -> SVector<f64, D> {
        evaluate_polynomial(&self.curve, x)
    }

    fn derivative(&self, x: f64) -> SVector<f64, D> {
        let deriv = derivative_polynomial(&self.curve);
        evaluate_polynomial(&deriv, x)
    }

    fn derivative2nd(&self, x: f64) -> SVector<f64, D> {
        let deriv = derivative2nd_polynomial(&self.curve);
        evaluate_polynomial(&deriv, x)
    }

    fn range(&self) -> Interval {
        Interval::new_closed(0.0, 1.0)
    }
}

#[derive(Debug, Clone)]
struct CubicSpline<const D: usize> {
    curve: Polynomial<D, 4>,
}

impl<const D: usize> CubicSpline<D> {
    fn new(
        start: SVector<f64, D>,
        end: SVector<f64, D>,
        start_vel: SVector<f64, D>,
        end_vel: SVector<f64, D>,
    ) -> CubicSpline<D> {
        let b = SMatrix::<f64, D, 4>::from_columns([start, end, start_vel, end_vel].as_ref())
            .transpose();

        let a_inv = SMatrix::<f64, 4, 4>::from_row_slice(
            [
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -3.0, 3.0, -2.0, -1.0, 2.0, -2.0, 1.0, 1.0,
            ]
            .as_ref(),
        );

        let curve = a_inv * b;

        return CubicSpline {
            curve: curve.transpose(),
        };
    }
}

impl<const D: usize> C2Spline<D> for CubicSpline<D> {
    fn evaluate(&self, x: f64) -> SVector<f64, D> {
        evaluate_polynomial(&self.curve, x)
    }

    fn derivative(&self, x: f64) -> SVector<f64, D> {
        let deriv = derivative_polynomial(&self.curve);
        evaluate_polynomial(&deriv, x)
    }

    fn derivative2nd(&self, x: f64) -> SVector<f64, D> {
        let deriv = derivative2nd_polynomial(&self.curve);
        evaluate_polynomial(&deriv, x)
    }

    fn range(&self) -> Interval {
        Interval::new_closed(0.0, 1.0)
    }
}

pub struct Waypoint<const D: usize> {
    pub position: SVec<D>,
    pub velocity: SVec<D>,
}

impl<const D: usize> Waypoint<D> {
    pub fn new(position: SVec<D>, velocity: SVec<D>) -> Waypoint<D> {
        Waypoint { position, velocity }
    }
}

impl<const D: usize> PartialEq for Waypoint<D> {
    fn eq(&self, other: &Self) -> bool {
        self.position == other.position && self.velocity == other.velocity
    }
}

struct ShiftedSpline<S: C2Spline<D>, const D: usize> {
    pub spline: S,
    pub interval: Interval,
}

impl<const D: usize, S: C2Spline<D>> ShiftedSpline<S, D> {
    pub fn new(spline: S, new_interval: Interval) -> ShiftedSpline<S, D> {
        ShiftedSpline {
            spline,
            interval: new_interval,
        }
    }

    fn fwrap_deriv(&self) -> f64 {
        let old_range = self.spline.range();
        (old_range.suprem.val() - old_range.infim.val())
            / (self.interval.suprem.val() - self.interval.infim.val())
    }
}

impl<const D: usize, S: C2Spline<D>> C2Spline<D> for ShiftedSpline<S, D> {
    fn evaluate(&self, x: f64) -> SVector<f64, D> {
        let old_range = self.spline.range();
        self.spline.evaluate(fwrap(
            x,
            self.interval.infim.val(),
            self.interval.suprem.val(),
            old_range.infim.val(),
            old_range.suprem.val(),
        ))
    }

    fn derivative(&self, x: f64) -> SVector<f64, D> {
        let old_range = self.spline.range();
        self.spline.derivative(fwrap(
            x,
            self.interval.infim.val(),
            self.interval.suprem.val(),
            old_range.infim.val(),
            old_range.suprem.val(),
        )) * self.fwrap_deriv()
    }

    fn derivative2nd(&self, x: f64) -> SVector<f64, D> {
        let old_range = self.spline.range();
        self.spline.derivative2nd(fwrap(
            x,
            self.interval.infim.val(),
            self.interval.suprem.val(),
            old_range.infim.val(),
            old_range.suprem.val(),
        )) * self.fwrap_deriv().powi(2)
    }

    fn range(&self) -> Interval {
        self.interval.clone()
    }
}

pub struct WaypointSpline<const D: usize> {
    head: ShiftedSpline<CubicSpline<D>, D>,
    rest: Vec<ShiftedSpline<QuarticSpline<D>, D>>,
}

impl<const D: usize> WaypointSpline<D> {
    fn find_spline(&self, t: f64) -> &dyn C2Spline<D> {
        if self.rest.is_empty()
            || self.head.range().is_lower_bound(t)
            || self.head.range().contains(t)
        {
            return &self.head;
        }

        if self.rest.last().unwrap().range().contains(t)
            || self.rest.last().unwrap().range().is_upper_bound(t)
        {
            return self.rest.last().unwrap();
        }

        let mut i = (self.rest.len() - 1) / 2;

        while i > 0 && i < self.rest.len() - 1 {
            if self.rest[i].range().contains(t) {
                return &self.rest[i];
            }

            if self.rest[i].range().is_lower_bound(t) {
                i /= 2;
            } else {
                i = (i + self.rest.len()) / 2;
            }
        }

        unreachable!()
    }
}

impl<const D: usize> C2Spline<D> for WaypointSpline<D> {
    fn evaluate(&self, x: f64) -> SVector<f64, D> {
        let spline = self.find_spline(x);
        spline.evaluate(x)
    }

    fn derivative(&self, x: f64) -> SVector<f64, D> {
        let spline = self.find_spline(x);
        spline.derivative(x)
    }

    fn derivative2nd(&self, x: f64) -> SVector<f64, D> {
        let spline = self.find_spline(x);
        spline.derivative2nd(x)
    }

    fn range(&self) -> Interval {
        let infim = self.head.range().infim;

        let suprem = match self.rest.last() {
            Some(spline) => spline.range().suprem,
            None => self.head.range().suprem,
        };

        Interval::new(infim, suprem)
    }
}

pub struct SplineBuilder<const D: usize> {
    pub waypoints: Vec<Waypoint<D>>,
}

impl<const D: usize> SplineBuilder<D> {
    pub fn new() -> SplineBuilder<D> {
        SplineBuilder {
            waypoints: Vec::new(),
        }
    }

    pub fn add_waypoint(&mut self, waypoint: Waypoint<D>) {
        self.waypoints.push(waypoint);
    }

    pub fn build(&self) -> Result<WaypointSpline<D>, ()> {
        if self.waypoints.len() < 2 {
            return Err(());
        }

        let head = CubicSpline::new(
            self.waypoints[0].position,
            self.waypoints[1].position,
            self.waypoints[0].velocity,
            self.waypoints[1].velocity,
        );
        let head_interval = Interval::new(
            IntervalEnd::Closed(head.range().infim.val()),
            IntervalEnd::Open(head.range().suprem.val()),
        );
        let head = ShiftedSpline::new(head, head_interval);

        let mut spline = WaypointSpline {
            head,
            rest: Vec::with_capacity(self.waypoints.len() - 2),
        };

        let mut prev_spline: &dyn C2Spline<D> = &spline.head;
        for i in 2..self.waypoints.len() {
            let waypoint = &self.waypoints[i];
            let prev_spline_end_point = prev_spline.evaluate(prev_spline.range().suprem.val());
            let prev_spline_end_vel = prev_spline.derivative(prev_spline.range().suprem.val());
            let prev_spline_end_accel = prev_spline.derivative2nd(prev_spline.range().suprem.val());

            let section = QuarticSpline::new(
                prev_spline_end_point,
                waypoint.position,
                prev_spline_end_vel,
                waypoint.velocity,
                prev_spline_end_accel,
            );

            let interval = Interval::new(
                IntervalEnd::Closed(prev_spline.range().suprem.val()),
                IntervalEnd::Open(prev_spline.range().suprem.val() + section.range().length()),
            );

            let section = ShiftedSpline::new(section, interval);
            spline.rest.push(section);
            prev_spline = &spline.rest[spline.rest.len() - 1];
        }

        Ok(spline)
    }
}

struct ReparametrizedSpline<S: C2Spline<D>, R: C2Spline<1>, const D: usize> {
    spline: S,
    reparam: R,
}

type DynReparametrizedSpline<S, const D: usize> = ReparametrizedSpline<S, Box<dyn C2Spline<1>>, D>;

impl<S: C2Spline<D>, R: C2Spline<1>, const D: usize> ReparametrizedSpline<S, R, D> {
    fn new(spline: S, reparam: R) -> ReparametrizedSpline<S, R, D> {
        ReparametrizedSpline { spline, reparam }
    }
}

impl<S: C2Spline<D>, R: C2Spline<1>, const D: usize> C2Spline<D> for ReparametrizedSpline<S, R, D> {
    fn evaluate(&self, x: f64) -> SVector<f64, D> {
        let t = self.reparam.evaluate(x)[0];
        self.spline.evaluate(t)
    }

    fn derivative(&self, x: f64) -> SVector<f64, D> {
        let t = self.reparam.evaluate(x)[0];
        self.spline.derivative(t) * self.reparam.derivative(x)
    }

    fn derivative2nd(&self, x: f64) -> SVector<f64, D> {
        let t = self.reparam.evaluate(x)[0];
        self.spline.derivative2nd(t) * self.reparam.derivative(x)[0].powi(2)
            + self.spline.derivative(t) * self.reparam.derivative2nd(x)
    }

    fn range(&self) -> Interval {
        self.reparam.range()
    }
}

#[cfg(test)]
mod tests {
    use quickcheck::{Arbitrary, Gen};
    use quickcheck_macros::quickcheck;
    use rand::distributions::Uniform;
    use rand::prelude::*;
    use rand::rngs::SmallRng;
    use std::io::Write;

    use super::*;
    use crate::assert_vec_eq;
    use crate::utils::test_utils::assert_feq;
    use crate::utils::Vec3;

    #[test]
    fn test_polynomial() {
        {
            let poly = Polynomial::<1, 3>::from_column_slice([1.0, 2.0, 1.0].as_ref());
            let result = evaluate_polynomial(&poly, -1.0);
            assert_eq!(result, SVector::<f64, 1>::new(0.0));
        }

        {
            let poly = Polynomial::<1, 3>::from_column_slice([2.0, 2.0, 1.0].as_ref());
            let result = evaluate_polynomial(&poly, -1.0);
            assert_eq!(result, SVector::<f64, 1>::new(1.0));
        }
    }

    #[test]
    fn test_derivative() {
        let poly = Polynomial::<1, 3>::from_column_slice([1.0, 2.0, 1.0].as_ref());
        let result = derivative_polynomial(&poly);
        assert_eq!(
            result,
            Polynomial::<1, 3>::from_column_slice([2.0, 2.0, 0.0].as_ref())
        );
    }

    #[test]
    fn test_derivative2nd() {
        let poly = Polynomial::<1, 3>::from_column_slice([1.0, 2.0, 1.0].as_ref());
        let result = derivative2nd_polynomial(&poly);
        assert_eq!(
            result,
            Polynomial::<1, 3>::from_column_slice([2.0, 0.0, 0.0].as_ref())
        );
    }

    #[derive(Debug, Clone, Copy)]
    struct ReasonableFloat(f64);

    impl Arbitrary for ReasonableFloat {
        fn arbitrary(g: &mut Gen) -> ReasonableFloat {
            let mut rng = SmallRng::seed_from_u64(u64::arbitrary(g));
            let range = Uniform::new_inclusive(-100_000.0, 100_000.0);
            ReasonableFloat(range.sample(&mut rng))
        }
    }

    #[derive(Debug, Clone)]
    struct ArbVec3(Vec3);

    impl Arbitrary for ArbVec3 {
        fn arbitrary(g: &mut Gen) -> ArbVec3 {
            let x = ReasonableFloat::arbitrary(g);
            let y = ReasonableFloat::arbitrary(g);
            let z = ReasonableFloat::arbitrary(g);
            ArbVec3(Vec3::new(x.0, y.0, z.0))
        }
    }

    #[test]
    fn test_quartic_spline() {
        let start = ReasonableFloat(0.0);
        let end = ReasonableFloat(1.0);
        let start_vel = ReasonableFloat(-0.5);
        let end_vel = ReasonableFloat(0.5);
        let start_acc = ReasonableFloat(0.0);

        let spline = QuarticSpline::new(
            SVector::<f64, 1>::new(start.0),
            SVector::<f64, 1>::new(end.0),
            SVector::<f64, 1>::new(start_vel.0),
            SVector::<f64, 1>::new(end_vel.0),
            SVector::<f64, 1>::new(start_acc.0),
        );

        let result = evaluate_polynomial(&spline.curve, 0.0);
        assert_feq(result[0], start.0);

        let result = evaluate_polynomial(&spline.curve, 1.0);
        assert_feq(result[0], end.0);

        let result_deriv = evaluate_polynomial(&derivative_polynomial(&spline.curve), 0.0);
        assert_feq(result_deriv[0], start_vel.0);

        let result_deriv = evaluate_polynomial(&derivative_polynomial(&spline.curve), 1.0);
        assert_feq(result_deriv[0], end_vel.0);

        let result_deriv_deriv = evaluate_polynomial(&derivative2nd_polynomial(&spline.curve), 0.0);
        assert_feq(result_deriv_deriv[0], start_acc.0);
    }

    #[quickcheck]
    fn quickcheck_quartic_spline(
        start: ReasonableFloat,
        end: ReasonableFloat,
        start_vel: ReasonableFloat,
        end_vel: ReasonableFloat,
        start_acc: ReasonableFloat,
    ) {
        let spline = QuarticSpline::new(
            SVector::<f64, 1>::new(start.0),
            SVector::<f64, 1>::new(end.0),
            SVector::<f64, 1>::new(start_vel.0),
            SVector::<f64, 1>::new(end_vel.0),
            SVector::<f64, 1>::new(start_acc.0),
        );

        let result = evaluate_polynomial(&spline.curve, 0.0);
        assert_feq(result[0], start.0);

        let result = evaluate_polynomial(&spline.curve, 1.0);
        assert_feq(result[0], end.0);

        let result_deriv = evaluate_polynomial(&derivative_polynomial(&spline.curve), 0.0);
        assert_feq(result_deriv[0], start_vel.0);

        let result_deriv = evaluate_polynomial(&derivative_polynomial(&spline.curve), 1.0);
        assert_feq(result_deriv[0], end_vel.0);

        let result_deriv_deriv = evaluate_polynomial(&derivative2nd_polynomial(&spline.curve), 0.0);
        assert_feq(result_deriv_deriv[0], start_acc.0);
    }

    #[quickcheck]
    fn quickcheck_quartic_spline_nd(
        start: ArbVec3,
        end: ArbVec3,
        start_vel: ArbVec3,
        end_vel: ArbVec3,
        start_acc: ArbVec3,
    ) {
        let spline = QuarticSpline::new(start.0, end.0, start_vel.0, end_vel.0, start_acc.0);

        let result = evaluate_polynomial(&spline.curve, 0.0);
        assert_feq((result - start.0).norm(), 0.0);

        let result = evaluate_polynomial(&spline.curve, 1.0);
        assert_feq((result - end.0).norm(), 0.0);

        let result_deriv = evaluate_polynomial(&derivative_polynomial(&spline.curve), 0.0);
        assert_feq((result_deriv - start_vel.0).norm(), 0.0);

        let result_deriv = evaluate_polynomial(&derivative_polynomial(&spline.curve), 1.0);
        assert_feq((result_deriv - end_vel.0).norm(), 0.0);

        let result_deriv_deriv = evaluate_polynomial(&derivative2nd_polynomial(&spline.curve), 0.0);
        assert_feq((result_deriv_deriv - start_acc.0).norm(), 0.0);
    }

    #[test]
    fn test_cube_spline() {
        let start = ReasonableFloat(0.0);
        let end = ReasonableFloat(1.0);
        let start_vel = ReasonableFloat(-0.5);
        let end_vel = ReasonableFloat(0.5);

        let spline = CubicSpline::new(
            SVector::<f64, 3>::new(start.0, 0.0, 0.0),
            SVector::<f64, 3>::new(end.0, 0.0, 0.0),
            SVector::<f64, 3>::new(start_vel.0, 0.0, 0.0),
            SVector::<f64, 3>::new(end_vel.0, 0.0, 0.0),
        );

        let result = spline.evaluate(0.0);
        assert_feq(result[0], start.0);

        let result = spline.evaluate(1.0);
        assert_feq(result[0], end.0);

        let result_deriv = spline.derivative(0.0);
        assert_feq(result_deriv[0], start_vel.0);

        let result_deriv = spline.derivative(1.0);
        assert_feq(result_deriv[0], end_vel.0);
    }

    #[quickcheck]
    fn quickcheck_cube_spline(start: ArbVec3, end: ArbVec3, start_vel: ArbVec3, end_vel: ArbVec3) {
        let spline = CubicSpline::new(start.0, end.0, start_vel.0, end_vel.0);

        let result = spline.evaluate(0.0);
        assert_feq((result - start.0).norm(), 0.0);

        let result = spline.evaluate(1.0);
        assert_feq((result - end.0).norm(), 0.0);

        let result_deriv = spline.derivative(0.0);
        assert_feq((result_deriv - start_vel.0).norm(), 0.0);

        let result_deriv = spline.derivative(1.0);
        assert_feq((result_deriv - end_vel.0).norm(), 0.0);
    }

    #[quickcheck]
    fn quickcheck_shifted_spline_offset(
        start: ArbVec3,
        end: ArbVec3,
        start_vel: ArbVec3,
        end_vel: ArbVec3,
        start_acc: ArbVec3,
        offset: ReasonableFloat,
    ) {
        let spline = QuarticSpline::new(start.0, end.0, start_vel.0, end_vel.0, start_acc.0);
        let interval = Interval::new_closed(
            spline.range().infim.val() + offset.0,
            spline.range().suprem.val() + offset.0,
        );
        let shifted_spline = ShiftedSpline::new(spline.clone(), interval);

        for t in (0..=100).map(|i| i as f64 / 100.0) {
            let result = shifted_spline.evaluate(t + offset.0);
            let expected = spline.evaluate(t);

            assert_vec_eq!(result, expected);

            let result = shifted_spline.derivative(t + offset.0);
            let expected = spline.derivative(t);
            assert!(
                (result - expected).norm() < 1e-4,
                "Expected: {:?}, got: {:?}",
                expected,
                result
            );

            let result = shifted_spline.derivative2nd(t + offset.0);
            let expected = spline.derivative2nd(t);
            assert!(
                (result - expected).norm() < 1e-4,
                "Expected: {:?}, got: {:?}",
                expected,
                result
            );
        }
    }
    #[quickcheck]
    fn quickcheck_shifted_spline(
        start: ArbVec3,
        end: ArbVec3,
        start_vel: ArbVec3,
        end_vel: ArbVec3,
        start_acc: ArbVec3,
        a: ReasonableFloat,
        b: ReasonableFloat,
    ) {
        if (a.0 - b.0).abs() < 1e-4 {
            return;
        }
        let (a, b) = if a.0 < b.0 { (a, b) } else { (b, a) };
        let spline = QuarticSpline::new(start.0, end.0, start_vel.0, end_vel.0, start_acc.0);
        let interval = Interval::new_closed(a.0, b.0);
        let shifted_spline = ShiftedSpline::new(spline, interval);

        let result = shifted_spline.evaluate(a.0);
        assert_vec_eq!(result, start.0);

        let result = shifted_spline.evaluate(b.0);
        assert_vec_eq!(result, end.0);

        let result = shifted_spline.derivative(a.0);
        assert_vec_eq!(result, start_vel.0);

        let result = shifted_spline.derivative(b.0);
        assert_vec_eq!(result, end_vel.0);

        let result = shifted_spline.derivative2nd(a.0);
        assert_vec_eq!(result, start_acc.0);
    }
}
