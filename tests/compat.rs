use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-vcf-merge"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn have(tool: &str) -> bool {
    Command::new(tool)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// Data records only (drop ## header), for field-level comparison.
fn records(vcf: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(vcf)
        .lines()
        .filter(|l| !l.starts_with("##") && !l.is_empty())
        .map(str::to_owned)
        .collect()
}

#[test]
fn runs_with_fixture() {
    let out = ours().arg(golden("a.vcf")).output().expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

// Multi-sample merge must match `bcftools merge`: all inputs' samples as
// columns, union of sites, ./. fill for missing genotypes.
#[test]
fn merge_matches_bcftools() {
    if !have("bcftools") || !have("bgzip") || !have("tabix") {
        eprintln!("skipping: bcftools/bgzip/tabix not found");
        return;
    }
    let dir = std::env::temp_dir().join("rsomics-vcf-merge-compat");
    let _ = std::fs::create_dir_all(&dir);
    let prep = |name: &str| -> PathBuf {
        let plain = dir.join(name);
        std::fs::copy(golden(name), &plain).unwrap();
        let gz = dir.join(format!("{name}.gz"));
        let g = std::fs::File::create(&gz).unwrap();
        assert!(
            Command::new("bgzip")
                .arg("-c")
                .arg(&plain)
                .stdout(g)
                .status()
                .unwrap()
                .success()
        );
        assert!(
            Command::new("tabix")
                .args(["-fp", "vcf"])
                .arg(&gz)
                .status()
                .unwrap()
                .success()
        );
        gz
    };
    let a_gz = prep("s1.vcf");
    let b_gz = prep("s2.vcf");

    let ours_out = ours()
        .arg(golden("s1.vcf"))
        .arg(golden("s2.vcf"))
        .output()
        .unwrap();
    let bcf_out = Command::new("bcftools")
        .arg("merge")
        .arg(&a_gz)
        .arg(&b_gz)
        .output()
        .unwrap();
    assert!(bcf_out.status.success());

    assert_eq!(records(&ours_out.stdout), records(&bcf_out.stdout));
}
