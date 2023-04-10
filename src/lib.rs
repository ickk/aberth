#![no_std]
#![doc = include_str!("../README.md")]

use arrayvec::ArrayVec;
use core::iter::zip;
use num_complex::Complex;
use num_traits::{
  cast,
  float::{Float, FloatConst},
  MulAdd,
};

/// Find all of the roots of a polynomial using Aberth's method.
///
/// Polynomial of the form f(x) = a + b*x + c*x^2 + d*x^3 + ...
///
/// `polynomial` is a slice containing the coefficients [a, b, c, d, ...]
///
/// When two successive iterations produce roots with less than `epsilon`
/// delta, the roots are returned.
pub fn aberth<
  const TERMS: usize,
  F: Float + FloatConst + MulAdd<Output = F>,
>(
  polynomial: &[F; TERMS],
  epsilon: F,
) -> Result<ArrayVec<Complex<F>, TERMS>, &'static str> {
  let dydx = &derivative(polynomial);
  let mut zs: ArrayVec<Complex<F>, TERMS> = initial_guesses(polynomial);
  let mut new_zs = zs.clone();

  'iteration: for _ in 0..100 {
    let mut converged = true;

    for i in 0..zs.len() {
      let p_of_z = sample_polynomial(polynomial, zs[i]);
      let dydx_of_z = sample_polynomial(dydx, zs[i]);
      let sum = (0..zs.len())
        .filter(|&k| k != i)
        .fold(Complex::<F>::ZERO(), |acc, k| {
          acc + Complex::<F>::ONE() / (zs[i] - zs[k])
        });

      let new_z = zs[i] + p_of_z / (p_of_z * sum - dydx_of_z);
      new_zs[i] = new_z;

      if new_z.re.is_nan()
        || new_z.im.is_nan()
        || new_z.re.is_infinite()
        || new_z.im.is_infinite()
      {
        break 'iteration;
      }

      if !new_z.approx_eq(zs[i], epsilon) {
        converged = false;
      }
    }
    if converged {
      return Ok(new_zs);
    }
    core::mem::swap(&mut zs, &mut new_zs);
  }
  Err("Failed to converge.")
}

// Initial guesses using the method from "Iteration Methods for Finding all
// Zeros of a Polynomial Simultaneously" by Oliver Aberth.
fn initial_guesses<
  const TERMS: usize,
  F: Float + FloatConst + MulAdd<Output = F>,
>(
  polynomial: &[F; TERMS],
) -> ArrayVec<Complex<F>, TERMS> {
  // the degree of the polynomial
  let n = polynomial.len() - 1;
  let n_f = unsafe { cast(n).unwrap_unchecked() };
  // convert polynomial to monic form
  let mut monic: ArrayVec<F, TERMS> = ArrayVec::new();
  for c in polynomial {
    // SAFETY: we push only as many values as there are terms.
    unsafe { monic.push_unchecked(*c / polynomial[n]) };
  }
  // let a = - c_1 / n
  let a = -monic[n - 1] / n_f;
  // let z = w + a,
  let p_of_w = {
    // we can recycle monic on the fly.
    for coefficient_index in 0..=n {
      let c = monic[coefficient_index];
      monic[coefficient_index] = F::zero();
      for ((index, power), pascal) in zip(
        zip(0..=coefficient_index, (0..=coefficient_index).rev()),
        PascalRowIter::new(coefficient_index as u32),
      ) {
        let pascal: F = unsafe { cast(pascal).unwrap_unchecked() };
        monic[index] =
          MulAdd::mul_add(c, pascal * a.powi(power as i32), monic[index]);
      }
    }
    monic
  };
  // convert P(w) into S(w)
  let s_of_w = {
    let mut p = p_of_w;
    // skip the last coefficient
    for i in 0..n {
      p[i] = -p[i].abs()
    }
    p
  };
  // find r_0
  let mut int = F::one();
  let r_0 = loop {
    let s_at_r0 = sample_polynomial(&s_of_w, int.into());
    if s_at_r0.re > F::zero() {
      break int;
    }
    int = int + F::one();
  };
  drop(s_of_w);

  {
    let mut guesses: ArrayVec<Complex<F>, TERMS> = ArrayVec::new();

    let frac_2pi_n = F::TAU() / n_f;
    let frac_pi_2n = F::FRAC_PI_2() / n_f;

    for k in 0..n {
      let k_f = unsafe { cast(k).unwrap_unchecked() };
      let theta = MulAdd::mul_add(frac_2pi_n, k_f, frac_pi_2n);

      let real = MulAdd::mul_add(r_0, theta.cos(), a);
      let imaginary = r_0 * theta.sin();

      let val = Complex::new(real, imaginary);
      // SAFETY: we push 1 less values than there are terms.
      unsafe { guesses.push_unchecked(val) };
    }

    guesses
  }
}

/// An iterator over the numbers in a row of Pascal's Triangle.
pub struct PascalRowIter {
  n: u32,
  k: u32,
  previous: u32,
}

impl PascalRowIter {
  /// Create an iterator yielding the numbers in the nth row of Pascal's
  /// triangle.
  pub fn new(n: u32) -> Self {
    Self {
      n,
      k: 0,
      previous: 1,
    }
  }
}

impl Iterator for PascalRowIter {
  type Item = u32;

  fn next(&mut self) -> Option<Self::Item> {
    if self.k == 0 {
      self.k = 1;
      self.previous = 1;
      return Some(1);
    }
    if self.k > self.n {
      return None;
    }
    let new = self.previous * (self.n + 1 - self.k) / self.k;
    self.k += 1;
    self.previous = new;
    Some(new)
  }
}

/// Sample the polynomial at some value of `x` using Horner's method.
///
/// Polynomial of the form f(x) = a + b*x + c*x^2 + d*x^3 + ...
///
/// `coefficients` is a slice containing the coefficients [a, b, c, d, ...]
pub fn sample_polynomial<F: Float + MulAdd<Output = F>>(
  coefficients: &[F],
  x: Complex<F>,
) -> Complex<F> {
  let mut r = Complex::ZERO();
  for c in coefficients.iter().rev() {
    r = r.mul_add(x, c.into())
  }
  r
}

/// Compute the derivative of a polynomial.
///
/// Polynomial of the form f(x) = a + b*x + c*x^2 + d*x^3 + ...
///
/// `coefficients` is a slice containing the coefficients [a, b, c, d, ...]
/// starting from the coefficient of the x with degree 0.
pub fn derivative<const TERMS: usize, F: Float>(
  coefficients: &[F; TERMS],
) -> ArrayVec<F, TERMS> {
  coefficients
    .iter()
    .enumerate()
    .skip(1)
    .map(|(power, &coefficient)| {
      let p = unsafe { F::from(power).unwrap_unchecked() };
      p * coefficient
    })
    .collect()
}

/// Some extra methods for Complex numbers
trait ComplexExt<F: Float> {
  /// Cheap comparison of complex numbers to within some margin, epsilon.
  fn approx_eq(self, w: Self, epsilon: F) -> bool;
  /// The additive identity
  #[allow(non_snake_case)]
  fn ZERO() -> Self;
  /// The multiplicative identity
  #[allow(non_snake_case)]
  fn ONE() -> Self;
}

impl<F: Float> ComplexExt<F> for Complex<F> {
  #[inline]
  fn approx_eq(self, w: Complex<F>, epsilon: F) -> bool {
    (self.re - w.re).abs() < epsilon && (self.im - w.im).abs() < epsilon
  }

  #[inline]
  fn ZERO() -> Self {
    Self {
      re: F::zero(),
      im: F::zero(),
    }
  }

  #[inline]
  fn ONE() -> Self {
    Complex {
      re: F::one(),
      im: F::zero(),
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  const EPSILON: f32 = 0.000_05;
  const EPSILON_64: f64 = 0.000_000_000_005;

  fn unsorted_compare<F: Float>(
    zs: &[Complex<F>],
    ws: &[Complex<F>],
    epsilon: F,
  ) -> bool {
    zs.iter().fold(true, |acc, &z| {
      let w = ws.iter().find(|&&w| z.approx_eq(w, epsilon));
      acc && w.is_some()
    })
  }

  fn array_approx_eq<F: Float>(lhs: &[F], rhs: &[F], epsilon: F) -> bool {
    if lhs.len() != rhs.len() {
      return false;
    }
    for i in 0..lhs.len() {
      if (lhs[i] - rhs[i]).abs() > epsilon {
        return false;
      }
    }
    true
  }

  #[test]
  fn derivative() {
    use super::derivative;

    {
      let y = [0., 1., 2., 3., 4.];
      let dydx = derivative(&y);

      let expected = [1., 4., 9., 16.];
      assert!(array_approx_eq(&dydx, &expected, EPSILON));
    }

    {
      let y = [19., 2.3, 0., 8.3, 69.420];
      let dydx = derivative(&y);

      let expected = [2.3, 0., 24.9, 277.68];
      assert!(array_approx_eq(&dydx, &expected, EPSILON));
    }
  }

  #[test]
  fn sample_polynomial() {
    use super::sample_polynomial;

    {
      let y = [0., 1., 2., 3., 4.];

      let x_0 = 0.0.into();
      let y_0 = sample_polynomial(&y, x_0);
      let expected_0 = 0.0.into();
      assert!(y_0.approx_eq(expected_0, EPSILON));

      let x_1 = 1.0.into();
      let y_1 = sample_polynomial(&y, x_1);
      let expected_1 = 10.0.into();
      assert!(y_1.approx_eq(expected_1, EPSILON));

      let x_2 = (-1.0).into();
      let y_2 = sample_polynomial(&y, x_2);
      let expected_2 = 2.0.into();
      assert!(y_2.approx_eq(expected_2, EPSILON));

      let x_3 = 2.5.into();
      let y_3 = sample_polynomial(&y, x_3);
      let expected_3 = 218.125.into();
      assert!(y_3.approx_eq(expected_3, EPSILON));
    }

    {
      let y = [19., 2.3, 0., 8.3, 69.420];

      let x_0 = 0.0.into();
      let y_0 = sample_polynomial(&y, x_0);
      let expected_0 = 19.0.into();
      assert!(y_0.approx_eq(expected_0, EPSILON));

      let x_1 = 1.0.into();
      let y_1 = sample_polynomial(&y, x_1);
      let expected_1 = 99.02.into();
      assert!(y_1.approx_eq(expected_1, EPSILON));

      let x_2 = (-1.0).into();
      let y_2 = sample_polynomial(&y, x_2);
      let expected_2 = 77.82.into();
      assert!(y_2.approx_eq(expected_2, EPSILON));
    }
  }

  #[test]
  fn aberth() {
    use super::*;

    {
      let polynomial = [0., 1.];
      let roots = aberth(&polynomial, EPSILON).unwrap();
      assert!(roots[0].approx_eq(Complex::ZERO(), EPSILON));
    }

    {
      let polynomial = [1., 0., -1.];
      let roots = aberth(&polynomial, EPSILON).unwrap();
      let expected = [1.0.into(), (-1.0).into()];
      assert!(unsorted_compare(&roots, &expected, EPSILON));
    }

    {
      // x^3 -12x^2 + 39x - 28 = 0
      let polynomial = [-28., 39., -12., 1.];

      let roots = aberth(&polynomial, EPSILON).unwrap();
      let expected = [7.0.into(), 4.0.into(), 1.0.into()];
      assert!(unsorted_compare(&roots, &expected, EPSILON));
    }
    {
      // 2x^3 - 38x^2 + 228x - 432 = 0
      let polynomial = [-432., 228., -38., 2.];

      let roots = aberth(&polynomial, EPSILON).unwrap();
      let expected = [9.0.into(), 6.0.into(), 4.0.into()];
      assert!(unsorted_compare(&roots, &expected, EPSILON));
    }
    {
      // x^3 + 8 = 0
      let polynomial = [8., 0., 0., 1.];

      let roots = aberth(&polynomial, EPSILON).unwrap();
      let expected = [
        (-2.).into(),
        Complex::new(1., -3f32.sqrt()),
        Complex::new(1., 3f32.sqrt()),
      ];
      assert!(unsorted_compare(&roots, &expected, EPSILON));
    }
    {
      // 11x^9 + 4x^4 + 2x - 1 = 0
      let polynomial = [-1., 2., 0., 0., 4., 0., 0., 0., 0., 11.];

      let roots = aberth(&polynomial, EPSILON).unwrap();
      let expected = [
        (0.429326).into(),
        Complex::new(-0.802811, -0.229634),
        Complex::new(-0.802811, 0.229634),
        Complex::new(-0.344895, -0.842594),
        Complex::new(-0.344895, 0.842594),
        Complex::new(0.206720, -0.675070),
        Complex::new(0.206720, 0.675070),
        Complex::new(0.726324, -0.455503),
        Complex::new(0.726324, 0.455503),
      ];
      assert!(unsorted_compare(&roots, &expected, EPSILON));
    }
    {
      // 0 = - 20x^19 + 19x^18 - 18x^17 + 17x^16 - 16x^15
      //     + 15x^14 - 14x^13 + 13x^12 - 12x^11 + 11x^10
      //     - 10x^9  +  9x^8  -  8x^7  +  7x^6  -  6x^5
      //     +  5x^4  -  4x^3  +  3x^2  -  2x    +  1
      let polynomial = [
        1., -2., 3., -4., 5., -6., 7., -8., 9., -10., 11., -12., 13., -14.,
        15., -16., 17., -18., 19., -20.,
      ];

      let roots = aberth(&polynomial, EPSILON).unwrap();
      // found using wolfram alpha
      let expected = [
        0.834053.into(),
        Complex::new(-0.844_061, -0.321_794),
        Complex::new(-0.844_061, 0.321_794),
        Complex::new(-0.684_734, -0.550_992),
        Complex::new(-0.684_734, 0.550_992),
        Complex::new(-0.476_151, -0.721_437),
        Complex::new(-0.476_151, 0.721_437),
        Complex::new(-0.231_844, -0.822_470),
        Complex::new(-0.231_844, 0.822_470),
        Complex::new(0.028_207, -0.846_944),
        Complex::new(0.028_207, 0.846_944),
        Complex::new(0.281_692, -0.793_720),
        Complex::new(0.281_692, 0.793_720),
        Complex::new(0.506_511, -0.668_231),
        Complex::new(0.506_511, 0.668_231),
        Complex::new(0.682_933, -0.482_160),
        Complex::new(0.682_933, 0.482_160),
        Complex::new(0.795_421, -0.252_482),
        Complex::new(0.795_421, 0.252_482),
      ];
      assert!(unsorted_compare(&roots, &expected, EPSILON));
    }
  }

  #[test]
  fn aberth_f64() {
    use super::aberth;
    {
      // 0 = - 20x^19 + 19x^18 - 18x^17 + 17x^16 - 16x^15
      //     + 15x^14 - 14x^13 + 13x^12 - 12x^11 + 11x^10
      //     - 10x^9  +  9x^8  -  8x^7  +  7x^6  -  6x^5
      //     +  5x^4  -  4x^3  +  3x^2  -  2x    +  1
      let polynomial: [f64; 20] = [
        1., -2., 3., -4., 5., -6., 7., -8., 9., -10., 11., -12., 13., -14.,
        15., -16., 17., -18., 19., -20.,
      ];

      let roots = aberth(&polynomial, EPSILON_64).unwrap();
      let expected = [
        0.834_053_367_550.into(),
        Complex::new(-0.844_060_952_037, -0.321_793_977_746),
        Complex::new(-0.844_060_952_037, 0.321_793_977_746),
        Complex::new(-0.684_734_480_334, -0.550_992_054_369),
        Complex::new(-0.684_734_480_334, 0.550_992_054_369),
        Complex::new(-0.476_151_406_058, -0.721_436_901_065),
        Complex::new(-0.476_151_406_058, 0.721_436_901_065),
        Complex::new(-0.231_843_928_891, -0.822_470_497_825),
        Complex::new(-0.231_843_928_891, 0.822_470_497_825),
        Complex::new(0.028_207_047_127, -0.846_944_061_134),
        Complex::new(0.028_207_047_127, 0.846_944_061_134),
        Complex::new(0.281_691_706_643, -0.793_720_289_127),
        Complex::new(0.281_691_706_643, 0.793_720_289_127),
        Complex::new(0.506_511_447_570, -0.668_230_679_428),
        Complex::new(0.506_511_447_570, 0.668_230_679_428),
        Complex::new(0.682_933_030_868, -0.482_159_501_324),
        Complex::new(0.682_933_030_868, 0.482_159_501_324),
        Complex::new(0.795_420_851_336, -0.252_482_354_484),
        Complex::new(0.795_420_851_336, 0.252_482_354_484),
      ];
      assert!(unsorted_compare(&roots, &expected, EPSILON_64));
    }
  }

  #[test]
  fn pascal_triangle() {
    {
      let row = PascalRowIter::new(0)
        .collect::<ArrayVec<_, 1>>()
        .into_inner()
        .unwrap();
      let expected = [1];
      assert_eq!(row, expected);
    }
    {
      let row = PascalRowIter::new(1)
        .collect::<ArrayVec<_, 2>>()
        .into_inner()
        .unwrap();
      let expected = [1, 1];
      assert_eq!(row, expected);
    }
    {
      let row = PascalRowIter::new(2)
        .collect::<ArrayVec<_, 3>>()
        .into_inner()
        .unwrap();
      let expected = [1, 2, 1];
      assert_eq!(row, expected);
    }
    {
      let row = PascalRowIter::new(3)
        .collect::<ArrayVec<_, 4>>()
        .into_inner()
        .unwrap();
      let expected = [1, 3, 3, 1];
      assert_eq!(row, expected);
    }
    {
      let row = PascalRowIter::new(4)
        .collect::<ArrayVec<_, 5>>()
        .into_inner()
        .unwrap();
      let expected = [1, 4, 6, 4, 1];
      assert_eq!(row, expected);
    }
    {
      let row = PascalRowIter::new(5)
        .collect::<ArrayVec<_, 6>>()
        .into_inner()
        .unwrap();
      let expected = [1, 5, 10, 10, 5, 1];
      assert_eq!(row, expected);
    }
    {
      let row = PascalRowIter::new(6)
        .collect::<ArrayVec<_, 7>>()
        .into_inner()
        .unwrap();
      let expected = [1, 6, 15, 20, 15, 6, 1];
      assert_eq!(row, expected);
    }
    {
      let row = PascalRowIter::new(9)
        .collect::<ArrayVec<_, 10>>()
        .into_inner()
        .unwrap();
      let expected = [1, 9, 36, 84, 126, 126, 84, 36, 9, 1];
      assert_eq!(row, expected);
    }
  }
}
