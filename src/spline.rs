use std::usize;

use nalgebra::wrap;
use nalgebra::SMatrix;
use nalgebra::SVector;

// Polynomial of degree DEG as a function R -> R^D
type Polynomial<const D: usize, const DEG: usize> = SMatrix<f64, D, DEG>;

fn evaluate_polynomial<const D: usize, const DEG: usize>(
    poly: &Polynomial<D, DEG>,
    x: f64,
) -> SVector<f64, D> {
    let powers = SVector::from_iterator((0..DEG).map(|i| x.powi(i as i32)));
    poly * powers
}

fn polynom_deriv_matrix<const DEG: usize>() -> SMatrix<f64, DEG, DEG> {
    let mut matrix = SMatrix::<f64, DEG, DEG>::zeros();
    for i in 0..DEG - 1 {
        matrix[(i, i)] = (i + 1) as f64;
    }
    matrix
}

fn derivative_polynomial<const D: usize, const DEG: usize>(
    poly: &Polynomial<D, DEG>,
) -> Polynomial<D, DEG> {
    let mut result = SMatrix::<f64, D, DEG>::zeros();
    for i in 0..D {
        for j in 0..DEG - 1 {
            result[(i, j)] = poly[(i, j + 1)] * (j + 1) as f64;
        }
    }
    result
}

fn derivative2nd_polynomial<const D: usize, const DEG: usize>(
    poly: &Polynomial<D, DEG>,
) -> Polynomial<D, DEG> {
    let mut result = SMatrix::<f64, D, DEG>::zeros();
    for i in 0..D {
        for j in 0..DEG - 2 {
            result[(i, j)] = poly[(i, j + 2)] * (j + 1) as f64 * (j + 2) as f64;
        }
    }
    result
}

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

        #[rustfmt::skip]
        let a_inv = SMatrix::<f64, 5, 5>::from_row_slice(
            [
                 1.0,  0.0,  0.0,  0.0,  0.0,
                 0.0,  0.0,  1.0,  0.0,  0.0, 
                 0.0,  0.0,  0.0,  0.0,  0.5, 
                -4.0,  4.0, -3.0, -1.0, -1.0, 
                 3.0, -3.0,  2.0,  1.0,  0.5,
            ]
            .as_ref(),
        );

        let curve = a_inv * b;

        return QuarticSpline {
            curve: curve.transpose(),
        };
    }
}

#[cfg(test)]
mod tests {
    use quickcheck::{Arbitrary, Gen};
    use quickcheck_macros::quickcheck;
    use rand::distributions::Uniform;
    use rand::prelude::*;
    use rand::rngs::SmallRng;
    use splines::spline;

    use super::*;
    use crate::utils::test_utils::assert_feq;

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
}
