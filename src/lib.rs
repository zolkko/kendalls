//! [Kendall's tau rank correlation](https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient).
//! At this point this is basically a copy-paste
//! from [Apache Commons Math](http://commons.apache.org/proper/commons-math/) library with some
//! additions taken from [scipy](https://github.com/scipy/scipy)
//! and R [cor.test](https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/cor.test.R) function
//!
//! Example usage:
//! ```
//! let res = kendalls::tau_b(&[1, 2, 3], &[3, 4, 5]);
//! assert_eq!(res, Ok((1.0, 1.5666989036012806)));
//! ```
//! If you want to compute correlation, let's say, for `f64` type, then you will have to
//! provide either a custom comparator function or declare `Ord` trait for your custom floating point
//! numbers type (see [float](https://crates.io/crates/float) crate).
//!
//! ```
//! use std::cmp::Ordering;
//!
//! let (tau_b, _) = kendalls::tau_b_with_comparator(
//!            &[1.0, 2.0],
//!            &[3.0, 4.0],
//!            |a: &f64, b: &f64| a.partial_cmp(&b).unwrap_or(Ordering::Greater),
//!        ).unwrap();
//! assert_eq!(tau_b, 1.0);
//! ```
//!
//! The function will return an error if you pass empty arrays into it or `x` and `y` arrays'
//! dimensions are not equal.
use std::cmp::Ordering;
use std::error::Error as StdError;
use std::fmt::{Display, Error as FmtError, Formatter};
use std::result::Result;

#[derive(Debug, PartialEq)]
pub enum Error {
    DimensionMismatch { expected: usize, got: usize },
    InsufficientLength,
}

impl Display for Error {
    fn fmt(&self, f: &mut Formatter) -> Result<(), FmtError> {
        match self {
            Error::InsufficientLength => write!(f, "insufficient array length"),
            Error::DimensionMismatch { expected, got } => {
                write!(f, "dimension mismatch: {} != {}", expected, got)
            }
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
pub fn tau_b<T>(x: &[T], y: &[T]) -> Result<(f64, f64), Error>
where
    T: Ord + Clone + Default,
{
    tau_b_with_comparator(x, y, |a, b| a.cmp(b))
}

/// The same as `tau_b` but also allow to specify custom comparator for numbers for
/// which [Ord] trait is not defined.
#[allow(clippy::many_single_char_names)]
pub fn tau_b_with_comparator<T, F>(x: &[T], y: &[T], mut comparator: F) -> Result<(f64, f64), Error>
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

    let mut pairs: Vec<(T, T)> = x.iter().cloned().zip(y.iter().cloned()).collect();

    pairs.sort_by(|pair1, pair2| {
        let res = comparator(&pair1.0, &pair2.0);
        if res == Ordering::Equal {
            comparator(&pair1.1, &pair2.1)
        } else {
            res
        }
    });

    let mut v1_part_1 = 0usize;
    let mut v2_part_1 = 0isize;

    let mut tied_x_pairs = 0usize;
    let mut tied_xy_pairs = 0usize;
    let mut vt = 0usize;
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
            // TODO: refactor
            vt += consecutive_x_ties * (consecutive_x_ties - 1) * (2 * consecutive_x_ties + 5);
            v1_part_1 += consecutive_x_ties * (consecutive_x_ties - 1);

            let consecutive_x_ties_i = consecutive_x_ties as isize;
            v2_part_1 += consecutive_x_ties_i * (consecutive_x_ties_i - 1) * (consecutive_x_ties_i - 2);

            tied_x_pairs += sum(consecutive_x_ties - 1);
            consecutive_x_ties = 1;
            tied_xy_pairs += sum(consecutive_xy_ties - 1);
            consecutive_xy_ties = 1;

        }
    }

    vt += consecutive_x_ties * (consecutive_x_ties - 1) * (2 * consecutive_x_ties + 5);
    v1_part_1 += consecutive_x_ties * (consecutive_x_ties - 1);

    let consecutive_x_ties_i = consecutive_x_ties as isize;
    v2_part_1 += consecutive_x_ties_i * (consecutive_x_ties_i - 1) * (consecutive_x_ties_i - 2);

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

    let mut v1_part_2 = 0usize;
    let mut v2_part_2 = 0isize;
    let mut tied_y_pairs = 0usize;
    let mut consecutive_y_ties = 1usize;
    let mut vu = 0usize;

    for j in 1..n {
        let prev = &pairs[j - 1];
        let curr = &pairs[j];
        if curr.1 == prev.1 {
            consecutive_y_ties += 1;
        } else {
            // TODO: refactor
            vu += consecutive_y_ties * (consecutive_y_ties - 1) * (2 * consecutive_y_ties + 5);
            v1_part_2 += consecutive_y_ties * (consecutive_y_ties - 1);

            let consecutive_y_ties_i = consecutive_y_ties as isize;
            v2_part_2 += consecutive_y_ties_i * (consecutive_y_ties_i - 1) * (consecutive_y_ties_i - 2);
            
            tied_y_pairs += sum(consecutive_y_ties - 1);
            consecutive_y_ties = 1;
        }
    }

    vu += consecutive_y_ties * (consecutive_y_ties - 1) * (2 * consecutive_y_ties + 5);
    v1_part_2 += consecutive_y_ties * (consecutive_y_ties - 1);

    let consecutive_y_ties_i = consecutive_y_ties as isize;
    v2_part_2 += consecutive_y_ties_i * (consecutive_y_ties_i - 1) * (consecutive_y_ties_i - 2);

    tied_y_pairs += sum(consecutive_y_ties - 1);

    let v1 = (v1_part_1 * v1_part_2) as f64;
    let v2 = (v2_part_1 * v2_part_2) as f64;

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

    // non_tied_pairs_multiplied = ((n0 - n1) * (n0 - n2)).sqrt()
    let non_tied_pairs_multiplied = (num_pairs_f - tied_x_pairs_f) * (num_pairs_f - tied_y_pairs_f);

    let tau_b = concordant_minus_discordant / non_tied_pairs_multiplied.sqrt();

    // Significance
    let v0 = (n * (n - 1)) * (2 * n + 5);
    let n_f = n as f64;
            
    let var_s = (v0 - vt - vu) as f64 / 18.0 + v1 / (2.0 * n_f * (n_f - 1.0)) + v2 / (9.0 * n_f * (n_f - 1.0) * (n_f - 2.0));

    let s = tau_b * non_tied_pairs_multiplied.sqrt();
    let z = s / var_s.sqrt();

    // Limit range to fix computational errors
    Ok((tau_b.max(-1.0).min(1.0), z))
}

#[inline]
fn sum(n: usize) -> usize {
    n * (n + 1 as usize) / 2 as usize
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

        let (tau_b, z) = tau_b_with_comparator(&x, &y, |a: &f64, b: &f64| {
            a.partial_cmp(&b).unwrap_or(Ordering::Greater)
        }).unwrap();

        approx::assert_abs_diff_eq!(tau_b, -0.3762015410475098);
        approx::assert_abs_diff_eq!(z, -2.09764910068664);
    }

    #[test]
    fn shifted_test() {
        let comparator = |a: &f64, b: &f64| a.partial_cmp(&b).unwrap_or(Ordering::Greater);

        let x = &[1.0, 1.0, 2.0, 2.0, 3.0, 3.0];
        let y = &[1.0, 2.0, 2.0, 3.0, 3.0, 4.0];
        let (tau_b, z) = tau_b_with_comparator(&x[..], &y[..], comparator).unwrap();
        approx::assert_abs_diff_eq!(tau_b, 0.8006407690254358);
        approx::assert_abs_diff_eq!(z, 2.0526, epsilon = 0.0001);

        let x = &[12.0, 2.0, 1.0, 12.0, 2.0];
        let y = &[1.0, 4.0, 7.0, 1.0, 0.0];
        let (tau_b, z) = tau_b_with_comparator(&x[..], &y[..], comparator).unwrap();
        approx::assert_abs_diff_eq!(tau_b, -0.4714045207910316);
        approx::assert_abs_diff_eq!(z, -1.0742, epsilon = 0.0001);
    }

    #[test]
    fn simple_correlated_data() {
        let (tau_b, z) = tau_b(&[1, 2, 3], &[3, 4, 5]).unwrap();
        assert_eq!(tau_b, 1.0);
        approx::assert_abs_diff_eq!(z, 1.5666989036012806);
    }

    #[test]
    fn simple_correlated_reversed() {
        let (tau_b, z) = tau_b(&[1, 2, 3], &[5, 4, 3]).unwrap();
        assert_eq!(tau_b, -1.0);
        approx::assert_abs_diff_eq!(z, -1.5666989036012806);
    }

    #[test]
    fn simple_jumble() {
        let x = &[1.0, 2.0, 3.0, 4.0];
        let y = &[1.0, 3.0, 2.0, 4.0];

        // 6 pairs: (A,B) (A,C) (A,D) (B,C) (B,D) (C,D)
        // (B,C) is discordant, the other 5 are concordant
        let expected_tau_b = (5.0 - 1.0) / 6.0;
        let expected_z = 1.3587324409735149;

        assert_eq!(
            tau_b_with_comparator(x, y, |a: &f64, b: &f64| a
                .partial_cmp(&b)
                .unwrap_or(Ordering::Greater)),
            Ok((expected_tau_b, expected_z))
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
            Ok((0.0, 0.0))
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
}
