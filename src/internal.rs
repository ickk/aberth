use super::StopReason;
use core::{fmt::Debug, iter::zip};
use num_complex::{Complex, ComplexFloat};
use num_traits::{
  cast,
  float::{Float, FloatConst},
  identities::{One, Zero},
  MulAdd,
};

/// Find all of the roots of a polynomial using Aberth's method.
pub fn aberth_raw<F: Float + MulAdd<Output = F>>(
  polynomial: &[Complex<F>],
  dydx: &[Complex<F>],
  initial_guesses: &mut [Complex<F>],
  out: &mut [Complex<F>],
  max_iterations: u32,
  epsilon: F,
) -> StopReason {
  out.copy_from_slice(initial_guesses);
  let mut zs = initial_guesses;
  let mut new_zs = out;

  for iteration in 1..=max_iterations {
    let mut converged = true;

    for i in 0..zs.len() {
      let p_of_z = sample_polynomial(polynomial, zs[i]);
      let dydx_of_z = sample_polynomial(dydx, zs[i]);
      let sum = (0..zs.len())
        .filter(|&k| k != i)
        .fold(Complex::<F>::zero(), |acc, k| {
          acc + Complex::<F>::one() / (zs[i] - zs[k])
        });

      let new_z = zs[i] + p_of_z / (p_of_z * sum - dydx_of_z);
      new_zs[i] = new_z;

      if new_z.re.is_nan()
        || new_z.im.is_nan()
        || new_z.re.is_infinite()
        || new_z.im.is_infinite()
      {
        return StopReason::Failed(iteration);
      }

      if !new_z.approx_eq(zs[i], epsilon) {
        converged = false;
      }
    }
    if converged {
      return StopReason::Converged(iteration);
    }
    core::mem::swap(&mut zs, &mut new_zs);
  }
  StopReason::MaxIteration(max_iterations)
}

/// Sample the polynomial at some value of `x` using Horner's method.
///
/// Polynomial of the form f(x) = a + b*x + c*x^2 + d*x^3 + ...
///
/// `coefficients` is a slice containing the coefficients [a, b, c, d, ...]
pub(crate) fn sample_polynomial<F: Float + MulAdd<Output = F>>(
  coefficients: &[Complex<F>],
  x: Complex<F>,
) -> Complex<F> {
  #![allow(clippy::len_zero)]
  debug_assert!(coefficients.len() != 0);
  let mut r = Complex::zero();
  for &c in coefficients.iter().rev() {
    r = r.mul_add(x, c)
  }
  r
}

/// Compute the derivative of a polynomial.
///
/// Polynomial of the form f(x) = a + b*x + c*x^2 + d*x^3 + ...
///
/// `coefficients` is a slice containing the coefficients [a, b, c, d, ...]
/// starting from the coefficient of the x with degree 0.
pub(crate) fn derivative<F: Float>(
  polynomial: &[Complex<F>],
  out: &mut [Complex<F>],
) {
  polynomial
    .iter()
    .enumerate()
    .skip(1)
    .for_each(|(index, coefficient)| {
      // SAFETY: it's possible to cast any usize to a float
      let p = unsafe { F::from(index).unwrap_unchecked() };
      out[index - 1] = coefficient * p;
    })
}

// Initial guesses using the method from "Iteration Methods for Finding all
// Zeros of a Polynomial Simultaneously" by Oliver Aberth.
pub(crate) fn initial_guesses<
  F: Float + FloatConst + MulAdd<Output = F> + Debug,
>(
  polynomial: &[Complex<F>],
  out: &mut [Complex<F>],
) {
  // the degree of the polynomial
  let n = polynomial.len() - 1;
  // SAFETY: it's possible to cast any usize to a float
  let n_f: F = unsafe { cast(n).unwrap_unchecked() };
  // convert polynomial to monic form
  let monic = out;
  for (i, c) in polynomial.iter().enumerate() {
    monic[i] = c / polynomial[n]; // TODO: check this divide by zero
  }
  // let a = - c_1 / n
  let a: Complex<F> = -monic[n - 1] / n_f;
  // let z = w + a,
  let p_of_w = {
    // we can recycle monic on the fly.
    for coefficient_index in 0..=n {
      let c = monic[coefficient_index];
      monic[coefficient_index] = Complex::zero();
      for ((index, power), pascal) in zip(
        zip(0..=coefficient_index, (0..=coefficient_index).rev()),
        PascalRowIter::new(coefficient_index as u32),
      ) {
        // SAFETY: it's possible to cast any f64 to a float (including f32,
        // which will result in infinity if the value exceeds f32::MAX)
        let pascal: Complex<F> = unsafe { cast(pascal).unwrap_unchecked() };
        monic[index] =
          MulAdd::mul_add(c, pascal * a.powi(power as i32), monic[index]);
      }
    }
    monic
  };
  // convert P(w) into S(w)
  let s_of_w = {
    // skip the last coefficient
    p_of_w.iter_mut().take(n).for_each(|coefficient| {
      *coefficient = Complex::from(-coefficient.abs())
    });
    p_of_w
  };
  // find r_0
  let mut int = F::one();
  let r_0 = loop {
    let s_at_r0 = sample_polynomial(s_of_w, int.into());
    if s_at_r0.re > F::zero() {
      break int;
    }
    int = int + F::one();
  };

  {
    let guesses = s_of_w; // output

    let frac_2pi_n = F::TAU() / n_f;
    let frac_pi_2n = F::FRAC_PI_2() / n_f;

    for (k, guess) in guesses.iter_mut().enumerate().take(n) {
      // SAFETY: it's possible to cast any usize to a float
      let k_f = unsafe { cast(k).unwrap_unchecked() };
      let theta = MulAdd::mul_add(frac_2pi_n, k_f, frac_pi_2n);

      let real = r_0 * theta.cos();
      let imaginary = r_0 * theta.sin();

      let val = Complex::new(real, imaginary) + a;
      *guess = val;
    }
  }
}

/// An iterator over the numbers in a row of Pascal's Triangle.
pub(crate) struct PascalRowIter {
  n: u32,
  k: u32,
  previous: f64,
}

impl PascalRowIter {
  /// Create an iterator yielding the numbers in the nth row of Pascal's
  /// triangle.
  pub fn new(n: u32) -> Self {
    Self {
      n,
      k: 0,
      previous: 1.0,
    }
  }
}

impl Iterator for PascalRowIter {
  type Item = f64;

  fn next(&mut self) -> Option<Self::Item> {
    if self.k == 0 {
      self.k = 1;
      self.previous = 1.0;
      return Some(1.0);
    }
    if self.k > self.n {
      return None;
    }
    let new = self.previous * (self.n + 1 - self.k) as f64 / self.k as f64;
    self.k += 1;
    self.previous = new;
    Some(new)
  }
}

/// Some extra methods for Complex numbers
pub(crate) trait ComplexExt<F: Float> {
  fn approx_eq(self, w: Self, epsilon: F) -> bool;
}

impl<F: Float> ComplexExt<F> for Complex<F> {
  /// Cheap comparison of complex numbers to within some margin, epsilon.
  #[inline]
  fn approx_eq(self, w: Complex<F>, epsilon: F) -> bool {
    (self.re - w.re).abs() < epsilon && (self.im - w.im).abs() < epsilon
  }
}

pub(crate) use private::ComplexCoefficient;
mod private {
  use super::*;
  /// A trait to group real & complex float types into a single generic type
  pub trait ComplexCoefficient<F: Float>: Copy + Into<Complex<F>> {}
  impl ComplexCoefficient<f32> for f32 {}
  impl ComplexCoefficient<f64> for f64 {}
  impl ComplexCoefficient<f32> for Complex<f32> {}
  impl ComplexCoefficient<f64> for Complex<f64> {}
}
