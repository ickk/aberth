aberth
======
[crates.io](https://crates.io/crates/aberth) | [docs.rs](https://docs.rs/aberth)

An implementation of the
[Aberth-Ehrlich method](https://en.wikipedia.org/wiki/Aberth_method)
for finding the zeros of a polynomial.

Aberth's method uses an electrostatics analogy to model the approximations as
negative charges and the true zeros as positive charges. This enables
finding all complex roots simultaneously, converging cubically (worst-case it
converges linearly for zeros of multiplicity).

This crate is `#![no_std]` and tries to have minimal dependencies. It uses
[arrayvec](https://crates.io/crates/arrayvec) to avoid allocations, which will
be removed when rust stabilises support for const-generics.


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
const MAX_ITERATIONS: u32 = 10;

// 0 = -1 + 2x + 4x^4 + 11x^9
let polynomial = [-1., 2., 0., 0., 4., 0., 0., 0., 0., 11.];

let roots = aberth(&polynomial, MAX_ITERATIONS, EPSILON);
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

The above method does not require any allocation, instead doing all the
computation on the stack. It is generic over any size of polynomial, but the
size of the polynomial must be known at compile time.

If `std` is available then there is also an `AberthSolver` struct which
allocates some memory to support dynamically sized polynomials at run time.
This may also be good to use when you are dealing with polynomials with many
terms, as it uses the heap instead of blowing up the stack.

```rust
use aberth::AberthSolver;

let mut solver = AberthSolver::new();
solver.epsilon = 0.001;
solver.max_iterations = 10;

// 0 = -1 + 2x + 4x^3 + 11x^4
let a = [-1., 2., 0., 4., 11.];
// 0 = -28 + 39x^2 - 12x^3 + x^4
let b = [-28., 0., 39., -12., 1.];

for polynomial in [a, b] {
  let roots = solver.find_roots(&polynomial);
  // ...
}
```

Note that the returned values are not sorted in any particular order.

The coefficient of the highest degree term should not be zero.


`#![no_std]`
------------

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
