extern crate criterion;
extern crate kendalls;
extern crate rand;
extern crate rand_distr;

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use kendalls::tau_b_with_comparator;
use rand::prelude::*;
use rand_distr::StandardNormal;
use std::time::Duration;


fn generate_data(len: usize) -> (Vec<f64>, Vec<f64>) {
    let mut rng = rand::rng();
    let x: Vec<f64> = (0..len).map(|_| StandardNormal.sample(&mut rng)).collect();

    let y: Vec<f64> = (0..len).map(|_| StandardNormal.sample(&mut rng)).collect();
    (x, y)
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("kendalls");
    group.measurement_time(Duration::from_secs(10));
    group.warm_up_time(Duration::from_secs(3));

    for n in [100, 1000, 10_000] {
        let (x, y) = generate_data(n);
        group.bench_with_input(BenchmarkId::new("taub", n), &n, |b, &n| {
            b.iter(|| {
                tau_b_with_comparator(
                    black_box(&x),
                    black_box(&y),
                    black_box(|a: &f64, b: &f64| {
                        a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Greater)
                    }),
                )
            });
        });
    }
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
