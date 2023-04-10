use num_traits::{cast, float::Float};

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Complex<F: Float> {
  pub real: F,
  pub imaginary: F,
}

impl<F: Float> Complex<F> {
  #[inline]
  #[allow(non_snake_case)]
  pub fn ZERO() -> Self {
    Complex {
      real: F::zero(),
      imaginary: F::zero(),
    }
  }

  #[inline]
  #[allow(non_snake_case)]
  pub fn ONE() -> Self {
    Complex {
      real: F::one(),
      imaginary: F::zero(),
    }
  }

  #[inline]
  #[allow(non_snake_case)]
  pub fn J() -> Self {
    Complex {
      real: F::zero(),
      imaginary: F::one(),
    }
  }

  #[inline]
  pub const fn new(real: F, imaginary: F) -> Self {
    Complex { real, imaginary }
  }

  #[inline]
  pub fn from_polar(modulus: F, argument: F) -> Self {
    Complex {
      real: modulus * argument.cos(),
      imaginary: modulus * argument.sin(),
    }
  }

  #[inline]
  pub fn conjugate(self) -> Self {
    Complex {
      real: self.real,
      imaginary: -self.imaginary,
    }
  }

  #[inline]
  pub fn square(self) -> Self {
    Complex {
      real: self.real * self.real - self.imaginary * self.imaginary,
      imaginary: F::from(2.).unwrap() * self.real * self.imaginary,
    }
  }

  #[inline]
  pub fn cube(self) -> Self {
    let re_2 = self.real * self.real;
    let im_2 = self.imaginary * self.imaginary;
    Complex {
      real: self.real * re_2 - F::from(3.).unwrap() * self.real * im_2,
      imaginary: F::from(3.).unwrap() * re_2 * self.imaginary
        - im_2 * self.imaginary,
    }
  }

  #[inline]
  pub fn reciprocal(self) -> Self {
    let denominator = self.real * self.real + self.imaginary * self.imaginary;
    Complex {
      real: self.real / denominator,
      imaginary: -self.imaginary / denominator,
    }
  }

  #[inline]
  pub const fn real(self) -> F {
    self.real
  }

  #[inline]
  pub const fn imaginary(self) -> F {
    self.imaginary
  }

  #[inline]
  pub fn abs(self) -> F {
    self.modulus()
  }

  #[inline]
  pub fn modulus(self) -> F {
    (self.real * self.real + self.imaginary * self.imaginary).sqrt()
  }

  #[inline]
  pub fn arg(self) -> F {
    F::atan2(self.imaginary, self.real)
  }

  #[inline]
  pub fn sqrt(self) -> Self {
    let half_modulus = F::from(0.5).unwrap() * self.modulus();
    let half_real = F::from(0.5).unwrap() * self.real;

    Complex {
      real: (half_real + half_modulus).sqrt(),
      imaginary: (half_modulus - half_real).sqrt().copysign(self.imaginary),
    }
  }

  #[inline]
  pub fn cbrt(self) -> Self {
    let mod_cbrt = self.modulus().cbrt();
    let arg_div_n = self.arg() / F::from(3.).unwrap();

    Complex {
      real: mod_cbrt * arg_div_n.cos(),
      imaginary: mod_cbrt * arg_div_n.sin(),
    }
  }

  #[inline]
  pub fn powi(self, power: i32) -> Self {
    let mod_pow_n = self.modulus().powi(power);
    let arg_mul_n = self.arg() * cast(power).unwrap();
    Complex {
      real: mod_pow_n * arg_mul_n.cos(),
      imaginary: mod_pow_n * arg_mul_n.sin(),
    }
  }

  #[inline]
  pub fn approx_eq(self, w: Complex<F>, epsilon: F) -> bool {
    (self.real - w.real).abs() < epsilon
      && (self.imaginary - w.imaginary).abs() < epsilon
  }
}

impl<F: Float> From<F> for Complex<F> {
  #[inline]
  fn from(real: F) -> Self {
    Complex {
      real,
      imaginary: F::zero(),
    }
  }
}

impl<F: Float> core::ops::Add<Complex<F>> for Complex<F> {
  type Output = Self;

  #[inline]
  fn add(self, rhs: Self) -> Complex<F> {
    Complex {
      real: self.real + rhs.real,
      imaginary: self.imaginary + rhs.imaginary,
    }
  }
}

impl<F: Float> core::ops::Add<F> for Complex<F> {
  type Output = Self;

  #[inline]
  fn add(self, rhs: F) -> Complex<F> {
    Complex {
      real: self.real + rhs,
      imaginary: self.imaginary,
    }
  }
}

impl<F: Float> core::ops::Sub<Complex<F>> for Complex<F> {
  type Output = Self;

  #[inline]
  fn sub(self, rhs: Self) -> Complex<F> {
    Complex {
      real: self.real - rhs.real,
      imaginary: self.imaginary - rhs.imaginary,
    }
  }
}

impl<F: Float> core::ops::Sub<F> for Complex<F> {
  type Output = Self;

  #[inline]
  fn sub(self, rhs: F) -> Complex<F> {
    Complex {
      real: self.real - rhs,
      imaginary: self.imaginary,
    }
  }
}

impl<F: Float> core::ops::Mul<Complex<F>> for Complex<F> {
  type Output = Self;

  #[inline]
  fn mul(self, rhs: Self) -> Complex<F> {
    Complex {
      real: self.real * rhs.real - self.imaginary * rhs.imaginary,
      imaginary: self.real * rhs.imaginary + self.imaginary * rhs.real,
    }
  }
}

impl<F: Float> core::ops::Mul<F> for Complex<F> {
  type Output = Self;

  #[inline]
  fn mul(self, rhs: F) -> Complex<F> {
    Complex {
      real: self.real * rhs,
      imaginary: self.imaginary * rhs,
    }
  }
}

impl<F: Float> core::ops::Mul<&F> for Complex<F> {
  type Output = Self;

  #[inline]
  fn mul(self, rhs: &F) -> Complex<F> {
    self * *rhs
  }
}

impl<F: Float> core::ops::Div<Complex<F>> for Complex<F> {
  type Output = Self;

  #[inline]
  fn div(self, rhs: Self) -> Complex<F> {
    let denominator = rhs.real * rhs.real + rhs.imaginary * rhs.imaginary;
    let real =
      (self.real * rhs.real + self.imaginary * rhs.imaginary) / denominator;
    let imaginary =
      (self.imaginary * rhs.real - self.real * rhs.imaginary) / denominator;
    Complex { real, imaginary }
  }
}

impl<F: Float> core::ops::Div<F> for Complex<F> {
  type Output = Self;

  #[inline]
  fn div(self, rhs: F) -> Complex<F> {
    Complex {
      real: self.real / rhs,
      imaginary: self.imaginary / rhs,
    }
  }
}

// E0210 means we have to implement these for each float type
impl core::ops::Add<Complex<f32>> for f32 {
  type Output = Complex<f32>;

  #[inline]
  fn add(self, rhs: Complex<f32>) -> Complex<f32> {
    Complex {
      real: self + rhs.real,
      imaginary: rhs.imaginary,
    }
  }
}
impl core::ops::Add<Complex<f64>> for f64 {
  type Output = Complex<f64>;

  #[inline]
  fn add(self, rhs: Complex<f64>) -> Complex<f64> {
    Complex {
      real: self + rhs.real,
      imaginary: rhs.imaginary,
    }
  }
}

impl core::ops::Sub<Complex<f32>> for f32 {
  type Output = Complex<f32>;

  #[inline]
  fn sub(self, rhs: Complex<f32>) -> Complex<f32> {
    Complex {
      real: self - rhs.real,
      imaginary: -rhs.imaginary,
    }
  }
}
impl core::ops::Sub<Complex<f64>> for f64 {
  type Output = Complex<f64>;

  #[inline]
  fn sub(self, rhs: Complex<f64>) -> Complex<f64> {
    Complex {
      real: self - rhs.real,
      imaginary: -rhs.imaginary,
    }
  }
}

impl core::ops::Mul<Complex<f32>> for f32 {
  type Output = Complex<f32>;

  #[inline]
  fn mul(self, rhs: Complex<f32>) -> Complex<f32> {
    Complex {
      real: self * rhs.real,
      imaginary: self * rhs.imaginary,
    }
  }
}
impl core::ops::Mul<Complex<f64>> for f64 {
  type Output = Complex<f64>;

  #[inline]
  fn mul(self, rhs: Complex<f64>) -> Complex<f64> {
    Complex {
      real: self * rhs.real,
      imaginary: self * rhs.imaginary,
    }
  }
}

impl core::ops::Mul<Complex<f32>> for &f32 {
  type Output = Complex<f32>;

  #[inline]
  fn mul(self, rhs: Complex<f32>) -> Complex<f32> {
    (*self).mul(rhs)
  }
}
impl core::ops::Mul<Complex<f64>> for &f64 {
  type Output = Complex<f64>;

  #[inline]
  fn mul(self, rhs: Complex<f64>) -> Complex<f64> {
    (*self).mul(rhs)
  }
}

impl core::ops::Div<Complex<f32>> for f32 {
  type Output = Complex<f32>;

  #[inline]
  fn div(self, rhs: Complex<f32>) -> Complex<f32> {
    Complex::from(self).div(rhs)
  }
}
impl core::ops::Div<Complex<f64>> for f64 {
  type Output = Complex<f64>;

  #[inline]
  fn div(self, rhs: Complex<f64>) -> Complex<f64> {
    Complex::from(self).div(rhs)
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use float_cmp::assert_approx_eq;
  use core::f32::consts::{SQRT_2, PI, FRAC_PI_2, FRAC_PI_4};

  const EPSILON: f32 = 0.0001;

  #[test]
  fn from_polar() {
    let z = Complex::from_polar(1., FRAC_PI_4);
    let expected = Complex::new(SQRT_2 / 2., SQRT_2 / 2.);

    assert!(z.approx_eq(expected, EPSILON));
  }

  #[test]
  fn j() {
    let z = Complex::J();
    let expected = Complex::new(0., 1.);

    assert!(z.approx_eq(expected, EPSILON));
  }

  #[test]
  fn conjugate() {
    let z = Complex::new(1.3, 3.4);
    let w = z.conjugate();
    let expected = Complex::new(1.3, -3.4);

    assert!(w.approx_eq(expected, EPSILON));
  }

  #[test]
  fn square() {
    let z = Complex::from_polar(3., FRAC_PI_4);
    let z_squared = z.square();
    let expected = Complex::new(0., 9.);

    assert!(z_squared.approx_eq(expected, EPSILON));
  }

  #[test]
  fn cube() {
    let z = Complex::from_polar(3., FRAC_PI_2);
    let z_cubed = z.cube();
    let expected = Complex::new(0., -27.);

    assert!(z_cubed.approx_eq(expected, EPSILON));
  }

  #[test]
  fn reciprocal() {
    let z = Complex::from_polar(5., 3. * PI / 4.);
    let z_reciprocal = z.reciprocal();
    let expected = Complex::from_polar(0.2, -3. * PI / 4.);

    assert!(z_reciprocal.approx_eq(expected, EPSILON));
  }

  #[test]
  fn real() {
    let z = Complex::new(1., 2.);
    let real = z.real();
    let expected = 1.;

    assert_eq!(real, expected);
  }

  #[test]
  fn imaginary() {
    let z = Complex::new(1., 2.);
    let im = z.imaginary();
    let expected = 2.;

    assert_eq!(im, expected);
  }

  #[test]
  fn modulus() {
    let z = Complex::new(1., 1.);
    let modulus = z.modulus();
    let expected = SQRT_2;

    assert_approx_eq!(f32, modulus, expected);
  }

  #[test]
  fn arg() {
    let z = Complex::new(0., 1.);
    let arg = z.arg();
    let expected = FRAC_PI_2;

    assert_approx_eq!(f32, arg, expected);
  }

  #[test]
  fn sqrt() {
    let z = Complex::new(0., 1.);
    let sqrt = z.sqrt();
    let expected = Complex::new(SQRT_2 / 2., SQRT_2 / 2.);

    assert!(sqrt.approx_eq(expected, EPSILON));
  }

  #[test]
  fn cbrt() {
    let z = Complex::new(0., 1.);
    let cbrt = z.cbrt();
    let expected = Complex::new(3f32.sqrt() / 2., 0.5);

    assert!(cbrt.approx_eq(expected, EPSILON));
  }

  #[test]
  fn powi() {
    let z = Complex::new(SQRT_2 / 2., SQRT_2 / 2.);

    {
      let pow2 = z.powi(2);
      let expected = Complex::new(0., 1.);
      assert!(pow2.approx_eq(expected, EPSILON));
    }

    {
      let pow3 = z.powi(3);
      let expected = Complex::new(-SQRT_2 / 2., SQRT_2 / 2.);
      assert!(pow3.approx_eq(expected, EPSILON));
    }

    {
      let pow5 = z.powi(5);
      let expected = Complex::new(-SQRT_2 / 2., -SQRT_2 / 2.);
      assert!(pow5.approx_eq(expected, EPSILON));
    }

    {
      let pow8 = z.powi(8);
      let expected = Complex::new(1., 0.);
      assert!(pow8.approx_eq(expected, EPSILON));
    }
  }

  #[test]
  fn from_f32() {
    let z = Complex::from(13.);
    let expected = Complex::new(13., 0.);

    assert!(z.approx_eq(expected, EPSILON));
  }

  #[test]
  fn add() {
    let z = Complex::new(1.5, 4.0);
    let w = Complex::new(0.2, 0.3);
    let result = z + w;
    let expected = Complex::new(1.7, 4.3);

    assert!(result.approx_eq(expected, EPSILON));
  }

  #[test]
  fn add_f32() {
    {
      let z = Complex::new(1.5, 4.0);
      let result = z + 13.;
      let expected = Complex::new(14.5, 4.);
      assert!(result.approx_eq(expected, EPSILON));
    }
    {
      let z = Complex::new(1.5, 4.0);
      let result = 13_f32 + z;
      let expected = Complex::new(14.5, 4.);
      assert!(result.approx_eq(expected, EPSILON));
    }
  }

  #[test]
  fn sub() {
    let z = Complex::new(1.5, 4.0);
    let w = Complex::new(0.2, 0.3);
    let result = z - w;
    let expected = Complex::new(1.3, 3.7);

    assert!(result.approx_eq(expected, EPSILON));
  }

  #[test]
  fn sub_f32() {
    {
      let z = Complex::new(1.5, 4.0);
      let result = z - 13.;
      let expected = Complex::new(-11.5, 4.0);
      assert!(result.approx_eq(expected, EPSILON));
    }

    {
      let z = Complex::new(1.5, 4.0);
      let result = 13_f32 - z;
      let expected = Complex::new(11.5, -4.0);
      assert!(result.approx_eq(expected, EPSILON));
    }
  }

  #[test]
  fn mul() {
    let z = Complex::new(3., 3.);
    let w = Complex::new(SQRT_2 / 2., SQRT_2 / 2.);
    let result = z * w;
    let expected = Complex::new(0., 3. * SQRT_2);

    assert!(result.approx_eq(expected, EPSILON));
  }

  #[test]
  fn mul_f32() {
    {
      let z = Complex::new(3., 1.);
      let result = z * 1.5;
      let expected = Complex::new(4.5, 1.5);
      assert!(result.approx_eq(expected, EPSILON));
    }

    {
      let z = Complex::new(3., 1.);
      let result = 1.5_f32 * z;
      let expected = Complex::new(4.5, 1.5);
      assert!(result.approx_eq(expected, EPSILON));
    }
  }

  #[test]
  fn div() {
    let z = Complex::new(-SQRT_2, SQRT_2);
    let w = Complex::from_polar(2., FRAC_PI_4);
    let result = z / w;
    let expected = Complex::new(0., 1.);

    assert!(result.approx_eq(expected, EPSILON));
  }

  #[test]
  fn div_f32() {
    {
      let z = Complex::new(3., 1.);
      let result = z / 1.5;
      let expected = Complex::new(2., 2. / 3.);
      assert!(result.approx_eq(expected, EPSILON));
    }

    {
      let z = Complex::from_polar(1., FRAC_PI_4);
      let result = 2. / z;
      let expected = Complex::new(SQRT_2, -SQRT_2);
      assert!(result.approx_eq(expected, EPSILON));
    }
  }
}
