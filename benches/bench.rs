use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_vcf_merge(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-vcf-merge");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let a = manifest.join("tests/golden/a.vcf");
    let b = manifest.join("tests/golden/b.vcf");
    c.bench_function("rsomics-vcf-merge golden", |b_| {
        b_.iter(|| {
            let out = Command::new(black_box(bin))
                .args([a.to_str().unwrap(), b.to_str().unwrap()])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_vcf_merge);
criterion_main!(benches);
