use std::fmt::Display;

use nalgebra::{allocator::Allocator, DefaultAllocator, SVector, Vector, Vector2, Vector3, U1};

pub type Vec3 = Vector3<f64>;
pub type Vec2 = Vector2<f64>;
pub type SVec<const D: usize> = SVector<f64, D>;

pub fn fwrap(x: f64, in_min: f64, in_max: f64, out_min: f64, out_max: f64) -> f64 {
    (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min
}

#[derive(Debug, Clone, Copy)]
pub enum IntervalEnd {
    Open(f64),
    Closed(f64),
}

impl IntervalEnd {
    pub fn val(&self) -> f64 {
        match self {
            IntervalEnd::Open(x) => *x,
            IntervalEnd::Closed(x) => *x,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Interval {
    pub infim: IntervalEnd,
    pub suprem: IntervalEnd,
}

impl Interval {
    pub fn new(a: IntervalEnd, b: IntervalEnd) -> Self {
        if a.val() < b.val() {
            Self {
                infim: a,
                suprem: b,
            }
        } else {
            Self {
                infim: b,
                suprem: a,
            }
        }
    }

    pub fn new_closed(a: f64, b: f64) -> Self {
        Self::new(IntervalEnd::Closed(a), IntervalEnd::Closed(b))
    }

    pub fn new_open(a: f64, b: f64) -> Self {
        Self::new(IntervalEnd::Open(a), IntervalEnd::Open(b))
    }

    pub fn contains(&self, x: f64) -> bool {
        let greater_then_infim = match self.infim {
            IntervalEnd::Open(a) => x > a,
            IntervalEnd::Closed(a) => x >= a,
        };

        let less_then_suprem = match self.suprem {
            IntervalEnd::Open(a) => x < a,
            IntervalEnd::Closed(a) => x <= a,
        };

        greater_then_infim && less_then_suprem
    }

    pub fn is_lower_bound(&self, x: f64) -> bool {
        match self.infim {
            IntervalEnd::Open(a) => x <= a,
            IntervalEnd::Closed(a) => x < a,
        }
    }

    pub fn is_upper_bound(&self, x: f64) -> bool {
        match self.suprem {
            IntervalEnd::Open(a) => x >= a,
            IntervalEnd::Closed(a) => x > a,
        }
    }

    pub fn length(&self) -> f64 {
        self.suprem.val() - self.infim.val()
    }
}

impl Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.infim {
            IntervalEnd::Open(a) => write!(f, "]{}, ", a)?,
            IntervalEnd::Closed(a) => write!(f, "[{}, ", a)?,
        }

        write!(f, ",")?;

        match self.suprem {
            IntervalEnd::Open(a) => write!(f, "{}[", a),
            IntervalEnd::Closed(a) => write!(f, "{}]", a),
        }
    }
}

#[cfg(test)]
pub mod test_utils {
    use super::*;
    use std::io::Write;

    pub fn assert_feq(a: f64, b: f64) {
        if !float_compare(a, b, 1e-4) {
            assert_eq!(a, b);
        }
    }

    #[macro_export]
    macro_rules! assert_vec_eq {
        ($a:expr, $b:expr) => {
            assert!($a.norm() - $b.norm() < 1e-4, "Vectors are not equal: {:?} != {:?}", $a, $b);
        };
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
