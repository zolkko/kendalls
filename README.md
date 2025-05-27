# kendalls

[![crates.io](https://img.shields.io/crates/v/kendalls.svg)](https://crates.io/crates/kendalls)
[![docs.rs](https://docs.rs/kendalls/badge.svg)](https://docs.rs/kendalls)
[![codecov](https://codecov.io/gh/zolkko/kendalls/branch/master/graph/badge.svg)](https://codecov.io/gh/zolkko/kendalls)

[Kendall's rank correlation coefficient](https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient)

# Usage

Add this to your Cargo.toml:

```toml
[dependencies]
kendalls = "1.0.0"
```

and this to your crate root:
```rust
extern crate kendalls;
```

Example:

```rust
fn main() -> Result<(), kendalls::Error> {
    let (tau_b, significance) = kendalls::tau_b(&[1, 2, 3], &[3, 4, 5])?;
    assert_eq!(tau_b, 1.0);
    assert_eq!(significance, 1.5666989036012806);

    Ok(())
}
```
