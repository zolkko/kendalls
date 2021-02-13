//! [Kendall's tau rank correlation](https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient).
//! At this point this is basically a copy-paste
//! from [Apache Commons Math](http://commons.apache.org/proper/commons-math/) library with some
//! additions taken from [scipy](https://github.com/scipy/scipy).
//!
//! Example usage:
//! ```
//! let res = kendalls::tau_b(&[1, 2, 3], &[3, 4, 5]);
//! assert_eq!(res, Ok(1.0));
//! ```
//! If you want to compute correlation, let's say, for `f64` type, then you will have to
//! provide either a custom comparator function or declare `Ord` trait for your custom floating point
//! numbers type (see [float](https://crates.io/crates/float) crate).
//!
//! ```
//! use std::cmp::Ordering;
//!
//! let res = kendalls::tau_b_with_comparator(
//!            &[1.0, 2.0],
//!            &[3.0, 4.0],
//!            |a: &f64, b: &f64| a.partial_cmp(&b).unwrap_or(Ordering::Greater),
//!        );
//! assert_eq!(res, Ok(1.0));
//! ```
//!
//! The function will return an error if you pass empty arrays into it or `x` and `y` arrays'
//! dimensions are not equal.
use std::cmp::Ordering;
use std::error::Error as StdError;
use std::fmt::{Display, Error as FmtError, Formatter};
use std::result::Result;


#[derive(Debug)]
pub struct Correlation {
    corr: f64,
}

impl Correlation {
    pub fn correlation(&self) -> f64 {
        self.corr
    }

    pub fn pvalue(&self) -> f64 {
        // _kendall_p_exact(n, int(min(C, (n*(n-1))//2-C)))
        todo!()
    }
}

fn factorial(num: usize) -> f64 {
    (1..=num).into_iter().map(|x| x as f64).product()
}

fn cumsum<T>(input: &[T]) -> Vec<T>
where
    T: std::ops::AddAssign + Copy + Default,
{
    let zero = Default::default();
    input
        .iter()
        .scan(zero, |acc, &x| {
            *acc += x;
            Some(*acc)
        })
        .collect()
}

/// Exact p-value, see Maurice G. Kendall, "Rank Correlation Methods" (4th Edition), Charles Griffin & Co., 1970.
fn p_exact(n: usize, c: usize) -> Result<f64, Error> {
    let c4 = 4 * c;
    let n2 = n * (n - 1);

    let prob = if n == 0 {
        return Err(Error::NMustBePositive(n));
    } else if c < 0 || c4 > n2 {
        return Err(Error::CNotSatisfied(c));
    } else if n == 1 {
        1.0
    } else if n == 2 {
        1.0
    } else if c == 0 {
        if n < 171 {
            2.0 / factorial(n)
        } else {
            0.0
        }
    } else if c == 1 {
        if n < 172 {
            2.0 / factorial(n - 1)
        } else {
            0.0
        }
    } else if c4 == n2 {
        1.0
    } else if n < 171 {
        let mut new = Vec::with_capacity(c + 1);
        for i in 0..(c + 1) {
            if i < 2 {
                new.push(1.0);
            } else {
                new.push(0.0);
            }
        }

        for j in 3..(n + 1) {
            new = cumsum(&new[..]);
            if j <= c {
                // new[j:] -= new[:c + 1 - j]
            }
        }

        2.0 * new.iter().sum::<f64>() / factorial(n)
    } else {
        let mut new = Vec::with_capacity(c + 1);
        for i in 0..(c + 1) {
            if i < 2 {
                new.push(1.0);
            } else {
                new.push(0.0);
            }
        }

        for j in 3..(n + 1) {
            new = new
                .iter()
                .scan(0.0, |acc, &x| {
                    *acc += x;
                    Some(*acc / (j as f64))
                })
                .collect();
            if j <= c {
                // new[j:] -= new[:c+1-j]
            }
        }

        new.iter().sum()
    };

    Ok(prob.max(0.0).min(1.0))
}

#[derive(Debug, PartialEq)]
pub enum Error {
    DimensionMismatch { expected: usize, got: usize },
    InsufficientLength,
    NMustBePositive(usize),
    /// c must satisfy 0 <= 4c <= n(n-1)
    CNotSatisfied(usize),
}

impl Display for Error {
    fn fmt(&self, f: &mut Formatter) -> Result<(), FmtError> {
        match self {
            Error::InsufficientLength => write!(f, "insufficient array length"),
            Error::DimensionMismatch { expected, got } => write!(f, "dimension mismatch: {} != {}", expected, got),
            Error::NMustBePositive(n) => write!(f, "n ({}) must be positive", n),
            Error::CNotSatisfied(c) => write!(f, "c ({}) must satisfy 0 <= 4c <= n(n-1)", c),
        }
    }
}

impl StdError for Error {}

/// Implementation of Kendall's Tau-b rank correlation between two arrays.
///
/// The definition of Kendall’s tau that is used is:
///
/// `tau = (P - Q) / sqrt((P + Q + T) * (P + Q + U))`
///
/// where P is the number of concordant pairs, Q the number of discordant pairs, T the number of
/// ties only in x, and U the number of ties only in y. If a tie occurs for the same pair in
/// both x and y, it is not added to either T or U.
pub fn tau_b<T>(x: &[T], y: &[T]) -> Result<Correlation, Error>
where
    T: Ord + Clone + Default,
{
    tau_b_with_comparator(x, y, |a, b| a.cmp(b))
}

/// The same as `tau_b` but also allow to specify custom comparator for numbers for
/// which [Ord] trait is not defined.
#[allow(clippy::many_single_char_names)]
pub fn tau_b_with_comparator<T, F>(x: &[T], y: &[T], mut comparator: F) -> Result<Correlation, Error>
where
    T: PartialOrd + Clone + Default,
    F: FnMut(&T, &T) -> Ordering,
{
    if x.len() != y.len() {
        return Err(Error::DimensionMismatch {
            expected: x.len(),
            got: y.len(),
        });
    }

    if x.is_empty() {
        return Err(Error::InsufficientLength);
    }

    let n = x.len();

    let mut pairs: Vec<(T, T)> = Vec::with_capacity(n);
    for pair in x.iter().cloned().zip(y.iter().cloned()) {
        pairs.push(pair);
    }

    pairs.sort_by(|pair1, pair2| {
        let res = comparator(&pair1.0, &pair2.0);
        if res == Ordering::Equal {
            comparator(&pair1.1, &pair2.1)
        } else {
            res
        }
    });

    let mut tied_x_pairs = 0usize;
    let mut tied_xy_pairs = 0usize;
    let mut consecutive_x_ties = 1usize;
    let mut consecutive_xy_ties = 1usize;

    for i in 1..n {
        let prev = &pairs[i - 1];
        let curr = &pairs[i];
        if curr.0 == prev.0 {
            consecutive_x_ties += 1;
            if curr.1 == prev.1 {
                consecutive_xy_ties += 1;
            } else {
                tied_xy_pairs += sum(consecutive_xy_ties - 1);
                consecutive_xy_ties = 1;
            }
        } else {
            tied_x_pairs += sum(consecutive_x_ties - 1);
            consecutive_x_ties = 1;
            tied_xy_pairs += sum(consecutive_xy_ties - 1);
            consecutive_xy_ties = 1;
        }
    }

    tied_x_pairs += sum(consecutive_x_ties - 1);
    tied_xy_pairs += sum(consecutive_xy_ties - 1);

    let mut swaps = 0usize;
    let mut pairs_dest: Vec<(T, T)> = vec![(Default::default(), Default::default()); n];

    let mut segment_size = 1usize;
    while segment_size < n {
        for offset in (0..n).step_by(2 * segment_size) {
            let mut i = offset;
            let i_end = n.min(i + segment_size);
            let mut j = i_end;
            let j_end = n.min(j + segment_size);

            let mut copy_location = offset;
            while i < i_end || j < j_end {
                if i < i_end {
                    if j < j_end {
                        let a = &pairs[i].1;
                        let b = &pairs[j].1;
                        if comparator(a, b) == Ordering::Greater {
                            pairs_dest[copy_location] = pairs[j].clone();
                            j += 1;
                            swaps += i_end - i;
                        } else {
                            pairs_dest[copy_location] = pairs[i].clone();
                            i += 1;
                        }
                    } else {
                        pairs_dest[copy_location] = pairs[i].clone();
                        i += 1;
                    }
                } else {
                    pairs_dest[copy_location] = pairs[j].clone();
                    j += 1;
                }
                copy_location += 1;
            }
        }

        std::mem::swap(&mut pairs, &mut pairs_dest);

        segment_size <<= 1;
    }

    let mut tied_y_pairs = 0usize;
    let mut consecutive_y_ties = 1usize;

    for i in 1..n {
        let prev = &pairs[i - 1];
        let curr = &pairs[i];
        if curr.1 == prev.1 {
            consecutive_y_ties += 1;
        } else {
            tied_y_pairs += sum(consecutive_y_ties - 1);
            consecutive_y_ties = 1;
        }
    }

    tied_y_pairs += sum(consecutive_y_ties - 1);

    // to prevent overflow on subtraction
    let num_pairs_f: f64 = ((n * (n - 1)) as f64) / 2.0; // sum(n - 1).as_();
    let tied_x_pairs_f: f64 = tied_x_pairs as f64;
    let tied_y_pairs_f: f64 = tied_y_pairs as f64;
    let tied_xy_pairs_f: f64 = tied_xy_pairs as f64;
    let swaps_f: f64 = (2 * swaps) as f64;

    // Note that tot = con + dis + (xtie - ntie) + (ytie - ntie) + ntie
    //               = con + dis + xtie + ytie - ntie
    //
    //           C-D = tot - xtie - ytie + ntie - 2 * dis
    let concordant_minus_discordant =
        num_pairs_f - tied_x_pairs_f - tied_y_pairs_f + tied_xy_pairs_f - swaps_f;
    let non_tied_pairs_multiplied = (num_pairs_f - tied_x_pairs_f) * (num_pairs_f - tied_y_pairs_f);

    let tau = concordant_minus_discordant / non_tied_pairs_multiplied.sqrt();
    let corr = tau.clamp(-1.0, 1.0);

    // return result
    Ok(Correlation { corr, })
}

/// Calculate statistical significance: Z.
///
/// Typically any value greater than 1.96 is going to be statistically significant
/// against the Z-table with alpha set at 0.05.
pub fn significance(tau: f64, n: usize) -> f64 {

    let n_tmp: f64 = (n * (n - 1)) as f64;
    let deter: f64 = (2 * (2 * n + 5)) as f64;
    (3.0 * tau * n_tmp.sqrt()) / deter.sqrt()
}

#[inline]
fn sum(n: usize) -> usize {
    n * (n + 1usize) / 2usize
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn xy_consecutive_pair_test() {
        let x = vec![
            12.0, 14.0, 14.0, 17.0, 19.0, 19.0, 19.0, 19.0, 19.0, 20.0, 21.0, 21.0, 21.0, 21.0, 21.0,
            22.0, 23.0, 24.0, 24.0, 24.0, 26.0, 26.0, 27.0,
        ];
        let y = vec![
            11.0, 4.0, 4.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 4.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
        ];

        let tau = tau_b_with_comparator(&x, &y, |a: &f64, b: &f64| {
            a.partial_cmp(&b).unwrap_or(Ordering::Greater)
        }).unwrap();

        approx::assert_abs_diff_eq!(tau, -0.3762015410475098);
    }

    #[test]
    fn shifted_test() {
        let comparator = |a: &f64, b: &f64| a.partial_cmp(&b).unwrap_or(Ordering::Greater);

        let x = &[1.0, 1.0, 2.0, 2.0, 3.0, 3.0];
        let y = &[1.0, 2.0, 2.0, 3.0, 3.0, 4.0];
        let res = tau_b_with_comparator(&x[..], &y[..], comparator).unwrap();
        approx::assert_abs_diff_eq!(res, 0.8006407690254358);

        let x = &[12.0, 2.0, 1.0, 12.0, 2.0];
        let y = &[1.0, 4.0, 7.0, 1.0, 0.0];
        let res = tau_b_with_comparator(&x[..], &y[..], comparator).unwrap();
        approx::assert_abs_diff_eq!(res, -0.4714045207910316);
    }

    #[test]
    fn simple_correlated_data() {
        let res = tau_b(&[1, 2, 3], &[3, 4, 5]).unwrap();
        assert_eq!(res, 1.0);
    }

    #[test]
    fn simple_correlated_reversed() {
        let res = tau_b(&[1, 2, 3], &[5, 4, 3]).unwrap();
        assert_eq!(res, -1.0);
    }

    #[test]
    fn simple_jumble() {
        let x = &[1.0, 2.0, 3.0, 4.0];
        let y = &[1.0, 3.0, 2.0, 4.0];

        // 6 pairs: (A,B) (A,C) (A,D) (B,C) (B,D) (C,D)
        // (B,C) is discordant, the other 5 are concordant
        let expected = (5.0 - 1.0) / 6.0;

        assert_eq!(
            tau_b_with_comparator(x, y, |a: &f64, b: &f64| a
                .partial_cmp(&b)
                .unwrap_or(Ordering::Greater)),
            Ok(expected)
        );
    }

    #[test]
    fn balanced_jumble() {
        let x = [1.0, 2.0, 3.0, 4.0];
        let y = [1.0, 4.0, 3.0, 2.0];

        // 6 pairs: (A,B) (A,C) (A,D) (B,C) (B,D) (C,D)
        // (A,B) (A,C), (A,D) are concordant, the other 3 are discordant

        assert_eq!(
            tau_b_with_comparator(&x, &y, |a: &f64, b: &f64| a
                .partial_cmp(&b)
                .unwrap_or(Ordering::Greater)),
            Ok(0.0)
        );
    }

    #[test]
    fn fails_if_dimentions_does_not_match() {
        let res = tau_b(&[1, 2, 3], &[5, 4]);
        assert_eq!(
            res,
            Err(Error::DimensionMismatch {
                expected: 3,
                got: 2
            })
        );
    }

    #[test]
    fn fails_if_arrays_are_empty() {
        let res = tau_b::<i32>(&[], &[]);
        assert_eq!(res, Err(Error::InsufficientLength));
    }

    #[test]
    fn it_format_dimension_mismatch_error() {
        let error = Error::DimensionMismatch { expected: 2, got: 1 };
        assert_eq!("dimension mismatch: 2 != 1", format!("{}", error));
    }

    #[test]
    fn it_format_insufficient_length_error() {
        let error = Error::InsufficientLength {} ;
        assert_eq!("insufficient array length", format!("{}", error));
    }

    #[test]
    fn significance_computed_correctly_for_certain_values() {
        let res = significance(0.818, 12);
        assert!(res > 3.7009);
        assert!(res < 3.709);
    }
}
