#![doc = include_str!("../README.md")]
#![cfg_attr(not(any(feature = "std", test, doctest)), no_std)]

pub mod internal;
use internal::*;
pub use num_complex::Complex;

#[cfg(any(test, doctest))]
mod tests;

use arrayvec::ArrayVec;
use core::{fmt::Debug, ops::Deref};
use num_traits::{
  float::{Float, FloatConst},
  identities::Zero,
  MulAdd,
};

/// Find all of the roots of a polynomial using Aberth's method
///
/// Polynomial of the form `f(x) = a + b*x + c*x^2 + d*x^3 + ...`
///
/// `polynomial` is a slice containing the coefficients `[a, b, c, d, ...]`
///
/// When two successive iterations produce roots with less than `epsilon`
/// delta, the roots are returned.
pub fn aberth<
  const TERMS: usize,
  F: Float + FloatConst + MulAdd<Output = F> + Debug,
  C: ComplexCoefficient<F> + Into<Complex<F>>,
>(
  polynomial: &[C; TERMS],
  max_iterations: u32,
  epsilon: F,
) -> Roots<ArrayVec<Complex<F>, TERMS>> {
  let degree = TERMS - 1;

  let polynomial: &[Complex<F>; TERMS] = &polynomial.map(|v| v.into());

  let mut dydx: ArrayVec<_, TERMS> = ArrayVec::new_const();
  // SAFETY: we immediately populate every entry in dydx.
  unsafe {
    dydx.set_len(degree);
    derivative::<F>(polynomial, dydx.as_mut());
  }

  let mut guesses: ArrayVec<_, TERMS> = ArrayVec::new_const();
  // SAFETY: we immediately populate every entry in guesses.
  unsafe {
    guesses.set_len(TERMS);
    initial_guesses(polynomial, guesses.as_mut());
    guesses.set_len(degree);
  }

  let mut output: ArrayVec<_, TERMS> = ArrayVec::new_const();
  // SAFETY: we push 1 less elements than there are terms.
  unsafe {
    for _ in 0..degree {
      output.push_unchecked(Complex::zero());
    }
  }

  let stop_reason = aberth_raw(
    polynomial,
    dydx.as_ref(),
    guesses.as_mut(),
    output.as_mut(),
    max_iterations,
    epsilon,
  );

  Roots {
    roots: output,
    stop_reason,
  }
}

/// The roots of a polynomial
///
/// Dereferences to an array-slice containing `roots`.
///
/// `stop_reason` contains information for how the solver terminated and how
/// many iterations it took.
#[derive(Clone, Debug, PartialEq)]
pub struct Roots<Arr> {
  pub roots: Arr,
  pub stop_reason: StopReason,
}

impl<Arr> Deref for Roots<Arr> {
  type Target = Arr;

  fn deref(&self) -> &Arr {
    &self.roots
  }
}

/// The reason the solver terminated and the number of iterations it took.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum StopReason {
  /// converged to within the required precision
  Converged(/* iterations */ u32),
  /// reached the iteration limit
  MaxIteration(/* iterations */ u32),
  /// detected a NaN or Inf while iterating
  Failed(/* iterations */ u32),
}

#[cfg(feature = "std")]
pub use feature_std::AberthSolver;

#[cfg(feature = "std")]
mod feature_std {
  use super::Roots;
  use crate::internal::*;
  use core::fmt::Debug;
  use num_complex::Complex;
  use num_traits::{
    cast,
    float::{Float, FloatConst},
    identities::Zero,
    MulAdd,
  };

  impl<F: Clone> Roots<&[Complex<F>]> {
    /// Create an owned duplicate of `Roots` by allocating a `Vec` to hold the
    /// values
    pub fn to_owned(&self) -> Roots<Vec<Complex<F>>> {
      Roots {
        roots: self.roots.to_vec(),
        stop_reason: self.stop_reason,
      }
    }
  }

  /// A solver for polynomials with Float or ComplexFloat coefficients. Will
  /// find all complex-roots, using the Aberth-Ehrlich method.
  ///
  /// The solver allocates some memory, and will reuse this allocation for
  /// subsequent calls. This is good to use for polynomials of varying lengths,
  /// polynomials with many terms, and for use in hot-loops where you want to
  /// avoid repeated allocations.
  ///
  /// Note the returned solutions are not sorted in any particular order.
  ///
  /// usage example:
  ///
  /// ```rust
  /// use aberth::AberthSolver;
  ///
  /// let mut solver = AberthSolver::new();
  /// solver.epsilon = 0.001;
  /// solver.max_iterations = 10;
  ///
  /// // 11x^4 + 4x^3 + 2x - 1 = 0
  /// let polynomial_a = [-1., 2., 0., 4., 11.];
  /// // x^4 -12x^3 + 39x^2 - 28 = 0
  /// let polynomial_b = [-28., 0., 39., -12., 1.];
  ///
  /// for polynomial in [polynomial_a, polynomial_b] {
  ///   let roots = solver.find_roots(&polynomial);
  ///   // ...
  /// }
  /// ```
  ///
  /// If you want to hold onto the roots you previously found while reusing the
  /// solver, then you can create an owned version:
  /// ```rust
  /// use aberth::AberthSolver;
  ///
  /// let mut solver = AberthSolver::new();
  /// let roots_a = solver.find_roots(&[-1., 2., 0., 4., 11.]).to_owned();
  /// let roots_b = solver.find_roots(&[-28., 39., -12., 1.]);
  /// roots_a[0];
  /// ```
  /// or alternatively just copy the `.roots` field into a vec
  /// ```rust
  /// use aberth::{AberthSolver, Complex};
  ///
  /// let mut solver = AberthSolver::new();
  /// let roots_a: Vec<Complex<f32>> =
  ///   solver.find_roots(&[-1., 2., 0., 4., 11.]).to_vec();
  /// let roots_b = solver.find_roots(&[-28., 39., -12., 1.]);
  /// roots_a[0];
  /// ```
  #[derive(Debug, Clone)]
  pub struct AberthSolver<F>
  where
    F: Float,
  {
    pub max_iterations: u32,
    pub epsilon: F,
    data: Vec<Complex<F>>,
  }

  impl<F: Float + FloatConst + MulAdd<Output = F> + Default + Debug> Default
    for AberthSolver<F>
  {
    fn default() -> Self {
      AberthSolver::new()
    }
  }

  impl<F: Float + FloatConst + MulAdd<Output = F> + Default + Debug>
    AberthSolver<F>
  {
    pub fn new() -> Self {
      AberthSolver {
        max_iterations: 100,
        data: Vec::new(),
        epsilon: cast(0.001).unwrap(),
      }
    }

    /// Find all the complex roots of the polynomial
    ///
    /// Polynomial is given in the form `f(x) = a + b*x + c*x^2 + d*x^3 + ...`
    ///
    /// `polynomial` is a slice containing the coefficients `[a, b, c, d, ...]`
    pub fn find_roots<C: ComplexCoefficient<F>>(
      &mut self,
      polynomial: &[C],
    ) -> Roots<&[Complex<F>]> {
      let len = polynomial.len();
      let degree = len - 1;
      // ensure we have enough space allocated
      self
        .data
        .resize_with(len + degree + len + degree, Complex::zero);
      // get mutable slices to our data
      let (complex_poly, tail) = self.data.split_at_mut(len);
      let (dydx, tail) = tail.split_at_mut(degree);
      let (guesses, output) = tail.split_at_mut(len);

      // convert the polynomial to a complex type
      polynomial
        .iter()
        .enumerate()
        .for_each(|(i, &coefficient)| complex_poly[i] = coefficient.into());

      initial_guesses(complex_poly, guesses);
      let guesses = &mut guesses[0..degree];
      derivative(complex_poly, dydx);

      let stop_reason = aberth_raw(
        complex_poly,
        dydx,
        guesses,
        output,
        self.max_iterations,
        self.epsilon,
      );

      Roots {
        roots: output,
        stop_reason,
      }
    }
  }
}
