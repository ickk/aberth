use crate::*;
use num_complex::ComplexFloat;

const EPSILON: f32 = 0.000_05;
const EPSILON_64: f64 = 0.000_000_000_005;

fn array_approx_eq<F: Float + FloatConst>(
  lhs: &[Complex<F>],
  rhs: &[Complex<F>],
  epsilon: F,
) -> bool {
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

#[test]
fn derivative() {
  use super::derivative;

  {
    let y = [0., 1., 2., 3., 4.];
    let y = y.map(|v| Complex::from(v));
    let mut dydx = [Complex::zero(); 4];

    derivative(&y, &mut dydx);

    let expected = [1., 4., 9., 16.];
    let expected = expected.map(|v| Complex::from(v));
    assert!(array_approx_eq(&dydx, &expected, EPSILON));
  }

  {
    let y = [19., 2.3, 0., 8.3, 69.420];
    let y = y.map(|v| Complex::from(v));
    let mut dydx = [Complex::zero(); 4];

    derivative(&y, &mut dydx);

    let expected = [2.3, 0., 24.9, 277.68];
    let expected = expected.map(|v| Complex::from(v));
    assert!(array_approx_eq(&dydx, &expected, EPSILON));
  }
}

#[test]
fn sample_polynomial() {
  use super::sample_polynomial;

  {
    let y = [0., 1., 2., 3., 4.];
    let y = y.map(|v| Complex::from(v));

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
    let y = y.map(|v| Complex::from(v));

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
fn pascal_triangle() {
  {
    let row: Vec<_> = PascalRowIter::new(0).collect();
    let expected = [1];
    assert_eq!(row, expected);
  }
  {
    let row: Vec<_> = PascalRowIter::new(1).collect();
    let expected = [1, 1];
    assert_eq!(row, expected);
  }
  {
    let row: Vec<_> = PascalRowIter::new(2).collect();
    let expected = [1, 2, 1];
    assert_eq!(row, expected);
  }
  {
    let row: Vec<_> = PascalRowIter::new(3).collect();
    let expected = [1, 3, 3, 1];
    assert_eq!(row, expected);
  }
  {
    let row: Vec<_> = PascalRowIter::new(4).collect();
    let expected = [1, 4, 6, 4, 1];
    assert_eq!(row, expected);
  }
  {
    let row: Vec<_> = PascalRowIter::new(5).collect();
    let expected = [1, 5, 10, 10, 5, 1];
    assert_eq!(row, expected);
  }
  {
    let row: Vec<_> = PascalRowIter::new(6).collect();
    let expected = [1, 6, 15, 20, 15, 6, 1];
    assert_eq!(row, expected);
  }
  {
    let row: Vec<_> = PascalRowIter::new(9).collect();
    let expected = [1, 9, 36, 84, 126, 126, 84, 36, 9, 1];
    assert_eq!(row, expected);
  }
}

/// ```should_panic
/// use aberth::aberth;
///
/// let y: [f32; 0] = [];
/// let dydx = aberth(&y, 100, 0.1);
/// ```
fn _aberth_empty_array() {}

#[test]
fn aberth() {
  use super::*;

  {
    let polynomial = [0., 1.];

    let roots = aberth(&polynomial, 100, EPSILON);
    assert!(roots[0].approx_eq(Complex::zero(), EPSILON));
  }

  {
    let polynomial = [1., 0., -1.];

    let roots = aberth(&polynomial, 100, EPSILON);
    let expected = [1.0.into(), (-1.0).into()];
    assert!(unsorted_compare(&roots, &expected, EPSILON));
  }

  {
    // x^3 -12x^2 + 39x - 28 = 0
    let polynomial = [-28., 39., -12., 1.];

    let roots = aberth(&polynomial, 100, EPSILON);
    let expected = [7.0.into(), 4.0.into(), 1.0.into()];
    assert!(unsorted_compare(&roots, &expected, EPSILON));
  }
  {
    // 2x^3 - 38x^2 + 228x - 432 = 0
    let polynomial = [-432., 228., -38., 2.];

    let roots = aberth(&polynomial, 100, EPSILON);
    let expected = [9.0.into(), 6.0.into(), 4.0.into()];
    assert!(unsorted_compare(&roots, &expected, EPSILON));
  }
  {
    // x^3 + 8 = 0
    let polynomial = [8., 0., 0., 1.];

    let roots = aberth(&polynomial, 100, EPSILON);
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

    let roots = aberth(&polynomial, 100, EPSILON);
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
      1., -2., 3., -4., 5., -6., 7., -8., 9., -10., 11., -12., 13., -14., 15.,
      -16., 17., -18., 19., -20.,
    ];

    let roots = aberth(&polynomial, 100, EPSILON);
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
      1., -2., 3., -4., 5., -6., 7., -8., 9., -10., 11., -12., 13., -14., 15.,
      -16., 17., -18., 19., -20.,
    ];

    let roots = aberth(&polynomial, 100, EPSILON_64);
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

#[cfg(feature = "std")]
mod feature_std {
  use crate::*;

  #[test]
  fn aberth_solver() {
    use super::*;
    let mut solver = AberthSolver::new();
    solver.epsilon = EPSILON;

    {
      let polynomial = [0., 1.];

      let roots = solver.find_roots(&polynomial);
      let expected = [0.0.into()];
      assert!(unsorted_compare(&roots, &expected, EPSILON));
    }

    {
      let polynomial = [1., 0., -1.];

      let roots = solver.find_roots(&polynomial);
      let expected = [1.0.into(), (-1.0).into()];
      assert!(unsorted_compare(&roots, &expected, EPSILON));
    }

    {
      // x^3 -12x^2 + 39x - 28 = 0
      let polynomial = [-28., 39., -12., 1.];

      let roots = solver.find_roots(&polynomial);
      let expected = [7.0.into(), 4.0.into(), 1.0.into()];
      assert!(unsorted_compare(&roots, &expected, EPSILON));
    }
    {
      // 2x^3 - 38x^2 + 228x - 432 = 0
      let polynomial = [-432., 228., -38., 2.];

      let roots = solver.find_roots(&polynomial);
      let expected = [9.0.into(), 6.0.into(), 4.0.into()];
      assert!(unsorted_compare(&roots, &expected, EPSILON));
    }
    {
      // x^3 + 8 = 0
      let polynomial = [8., 0., 0., 1.];

      let roots = solver.find_roots(&polynomial);
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

      let roots = solver.find_roots(&polynomial);
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

      let roots = solver.find_roots(&polynomial);
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
  fn reuse_solver() {
    let mut solver = AberthSolver::new();
    let roots_a = solver.find_roots(&[-1., 2., 0., 4., 11.]).to_owned();
    let roots_b = solver.find_roots(&[-28., 39., -12., 1.]);
    roots_a.roots;
    roots_b.roots;
  }

  #[test]
  fn many_coefficients() {
    let polynomial = vec![
      1., -2., 3., -4., 5., -6., 7., -8., 9., -10., 11., -12., 13., -14., 15.,
      -16., 17., -18., 19., -20., 21., -22., 23., -24., 25., -26., 27., -28.,
      29., -30., 31., -32.,
    ];
    let mut solver = AberthSolver::new();
    // should complete without integer overflow
    solver.find_roots(&polynomial);
  }
}
