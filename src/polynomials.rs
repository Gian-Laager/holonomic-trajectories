use nalgebra::{SMatrix, SVector};

// Polynomial of degree DEG as a function R -> R^D
pub type Polynomial<const D: usize, const DEG: usize> = SMatrix<f64, D, DEG>;

pub fn evaluate_polynomial<const D: usize, const DEG: usize>(
    poly: &Polynomial<D, DEG>,
    x: f64,
) -> SVector<f64, D> {
    let powers = SVector::from_iterator((0..DEG).map(|i| x.powi(i as i32)));
    poly * powers
}

pub fn polynom_deriv_matrix<const DEG: usize>() -> SMatrix<f64, DEG, DEG> {
    let mut matrix = SMatrix::<f64, DEG, DEG>::zeros();
    for i in 0..DEG - 1 {
        matrix[(i, i)] = (i + 1) as f64;
    }
    matrix
}

pub fn derivative_polynomial<const D: usize, const DEG: usize>(
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

pub fn derivative2nd_polynomial<const D: usize, const DEG: usize>(
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
