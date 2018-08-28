# kendalls

[![Crates.io](https://img.shields.io/crates/d/kendalls.svg)](https://crates.io/crates/kendalls)
[![Travis CI Build Status](https://api.travis-ci.org/zolkko/kendalls.svg?branch=master)](https://travis-ci.org/zolkko/kendalls)
[![docs](https://img.shields.io/badge/docs-online-5023dd.svg)](https://docs.rs/kendalls/)

[Kendall rank correlation coefficient](https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient)

# Usage

Add this to your Cargo.toml:
```toml
[dependencies]
kendalls = "0.1.0"
```

and this to your crate root:
```rust
extern crate kendalls;
```

Example:
```rust
extern crate kendalls;

use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let correlation = kendalls::tau_b(&[1,2,3], &[3,2,1])?;
    println!("correlation = {:?} == -1.0", correlation);
    Ok(())
}
```
