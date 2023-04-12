aberth
======

An implementation of the
[Aberth-Ehrlich method](https://en.wikipedia.org/wiki/Aberth_method)
for finding the zeros of a polynomial.

Aberth's method uses an electrostatics analogy to model the approximations as
negative charges and the true zeros as positive charges. This enables
finding all complex roots simultaneously, converging cubically (worst-case it
converges linearly for zeros of multiplicity).

This crate is `#![no_std]` and tries to have minimal dependencies:
- [num-complex](https://crates.io/crates/num-complex) for Complex number types
- [num-traits](https://crates.io/crates/num-traits) to be generic over floating
point types.


Usage
-----

Add it to your project:
```sh
cargo add aberth
```

Specify the coefficients of your polynomial in an array in ascending order
and then call the `aberth` method on your polynomial.
```rust
use aberth::aberth;
const EPSILON: f32 = 0.001;

// 0 = -1 + 2x + 4x^4 + 11x^9
let polynomial = [-1., 2., 0., 0., 4., 0., 0., 0., 0., 11.];

let roots = aberth(&polynomial, EPSILON).unwrap();
// [
//   Complex { re:  0.4293261, im:  1.084202e-19 },
//   Complex { re:  0.7263235, im:  0.4555030 },
//   Complex { re:  0.2067199, im:  0.6750696 },
//   Complex { re: -0.3448952, im:  0.8425941 },
//   Complex { re: -0.8028113, im:  0.2296336 },
//   Complex { re: -0.8028113, im: -0.2296334 },
//   Complex { re: -0.3448952, im: -0.8425941 },
//   Complex { re:  0.2067200, im: -0.6750695 },
//   Complex { re:  0.7263235, im: -0.4555030 },
// ]
```

Note that the returned values are not sorted in any particular order.

The function returns an `Err` if it fails to converge after 100 cycles. This
is excessive. Even polynomials of degree 20 will converge after <20 iterations.
Using a larger epsilon can usually avoid these errors.

The coefficient of the highest degree term should not be zero or you will get
nonsense extra roots (probably at 0 + 0j).


`#![no_std]`
--------

To use in a `no_std` environment you must disable `default-features` and enable
the `libm` feature:
```toml
[dependencies]
aberth = { version = "0.0.4", default-features = false, features = ["libm"] }
```


License
-------

This crate is dual-licensed under either of
[Apache license, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0)
or the
[MIT license](http://opensource.org/licenses/MIT)
at your option.
