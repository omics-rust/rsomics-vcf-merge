use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

struct Record {
    chrom: String,
    pos: u64,
    ref_allele: String,
    alts: Vec<String>,
    qual: String,
    filter: String,
    format: String,
    sample_vals: Vec<String>,
}

/// Merge QUAL across records: max numeric, "." if none numeric (matches bcftools merge).
fn merge_qual(quals: &[&str]) -> String {
    let mut best: Option<f64> = None;
    for q in quals {
        if let Ok(v) = q.parse::<f64>() {
            best = Some(best.map_or(v, |b: f64| b.max(v)));
        }
    }
    best.map_or_else(
        || ".".to_string(),
        |v| {
            if v.fract() == 0.0 {
                format!("{}", v as i64)
            } else {
                format!("{v}")
            }
        },
    )
}

/// Merge FILTER: union of non-PASS/non-"." filters; else PASS if any PASS; else ".".
fn merge_filter(filters: &[&str]) -> String {
    let mut named: Vec<&str> = Vec::new();
    let mut any_pass = false;
    for f in filters {
        match *f {
            "." => {}
            "PASS" => any_pass = true,
            other => {
                for part in other.split(';') {
                    if !named.contains(&part) {
                        named.push(part);
                    }
                }
            }
        }
    }
    if !named.is_empty() {
        named.join(";")
    } else if any_pass {
        "PASS".to_string()
    } else {
        ".".to_string()
    }
}

/// Missing per-sample value for a FORMAT string: GT -> ./., other keys -> .
fn missing_for(format: &str) -> String {
    format
        .split(':')
        .map(|k| if k == "GT" { "./." } else { "." })
        .collect::<Vec<_>>()
        .join(":")
}

struct Source {
    lines: std::io::Lines<BufReader<File>>,
    n_samples: usize,
    current: Option<Record>,
}

impl Source {
    fn open(path: &Path) -> Result<(Self, Vec<String>, Vec<String>)> {
        let file = File::open(path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
        let mut lines = BufReader::new(file).lines();
        let mut meta: Vec<String> = Vec::new();
        let mut samples: Vec<String> = Vec::new();
        for line in lines.by_ref() {
            let line = line.map_err(RsomicsError::Io)?;
            if line.starts_with("##") {
                meta.push(line);
            } else if line.starts_with('#') {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() > 9 {
                    samples = f[9..].iter().map(|s| (*s).to_string()).collect();
                }
                break;
            }
        }
        let n = samples.len();
        let mut src = Source {
            lines,
            n_samples: n,
            current: None,
        };
        src.advance()?;
        Ok((src, meta, samples))
    }

    fn advance(&mut self) -> Result<()> {
        self.current = None;
        for line in self.lines.by_ref() {
            let line = line.map_err(RsomicsError::Io)?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8 {
                continue;
            }
            self.current = Some(Record {
                chrom: f[0].to_string(),
                pos: f[1].parse().unwrap_or(0),
                ref_allele: f[3].to_string(),
                alts: f[4].split(',').map(str::to_owned).collect(),
                qual: f.get(5).copied().unwrap_or(".").to_string(),
                filter: f.get(6).copied().unwrap_or(".").to_string(),
                format: f.get(8).copied().unwrap_or("GT").to_string(),
                sample_vals: f
                    .get(9..)
                    .map_or_else(Vec::new, |s| s.iter().map(|x| (*x).to_string()).collect()),
            });
            break;
        }
        Ok(())
    }
}

/// Remap one allele index from a source's local ALT numbering to the merged one.
fn remap_allele(a: &str, map: &[usize]) -> String {
    if a == "." {
        return ".".to_string();
    }
    match a.parse::<usize>() {
        Ok(0) => "0".to_string(),
        Ok(local) => map
            .get(local - 1)
            .map_or_else(|| a.to_string(), ToString::to_string),
        Err(_) => a.to_string(),
    }
}

/// Remap the GT subfield (first FORMAT field) of a sample value; pass the rest through.
fn remap_sample(val: &str, map: &[usize]) -> String {
    let mut parts = val.splitn(2, ':');
    let gt = parts.next().unwrap_or(".");
    let rest = parts.next();
    let mut out = String::new();
    let mut num = String::new();
    for ch in gt.chars() {
        if ch == '/' || ch == '|' {
            out.push_str(&remap_allele(&num, map));
            num.clear();
            out.push(ch);
        } else {
            num.push(ch);
        }
    }
    out.push_str(&remap_allele(&num, map));
    if let Some(r) = rest {
        out.push(':');
        out.push_str(r);
    }
    out
}

pub fn merge(inputs: &[&Path], output: &mut dyn Write) -> Result<u64> {
    if inputs.is_empty() {
        return Err(RsomicsError::InvalidInput("no input files".into()));
    }
    let mut out = BufWriter::with_capacity(256 * 1024, output);

    let mut sources: Vec<Source> = Vec::with_capacity(inputs.len());
    let mut all_samples: Vec<String> = Vec::new();
    let mut header_meta: Vec<String> = Vec::new();
    let mut contig_order: HashMap<String, usize> = HashMap::new();
    for (i, path) in inputs.iter().enumerate() {
        let (src, meta, samples) = Source::open(path)?;
        if i == 0 {
            for m in &meta {
                if let Some(id) = m.strip_prefix("##contig=<ID=") {
                    let name = id.split([',', '>']).next().unwrap_or("");
                    if !name.is_empty() {
                        let next = contig_order.len();
                        contig_order.entry(name.to_string()).or_insert(next);
                    }
                }
            }
            header_meta = meta;
        }
        all_samples.extend(samples);
        sources.push(src);
    }

    for m in &header_meta {
        writeln!(out, "{m}").map_err(RsomicsError::Io)?;
    }
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        .map_err(RsomicsError::Io)?;
    for s in &all_samples {
        write!(out, "\t{s}").map_err(RsomicsError::Io)?;
    }
    writeln!(out).map_err(RsomicsError::Io)?;

    let rank = |c: &str| contig_order.get(c).copied().unwrap_or(usize::MAX);
    let mut count: u64 = 0;

    loop {
        // find the minimum (contig, pos) across sources with a pending record
        let mut min: Option<(usize, u64)> = None;
        for src in &sources {
            if let Some(r) = &src.current {
                let key = (rank(&r.chrom), r.pos);
                if min.is_none_or(|m| key < m) {
                    min = Some(key);
                }
            }
        }
        let Some((min_rank, min_pos)) = min else {
            break;
        };

        // gather every source positioned at this site
        let ref_allele = sources
            .iter()
            .find(|s| {
                s.current
                    .as_ref()
                    .is_some_and(|r| rank(&r.chrom) == min_rank && r.pos == min_pos)
            })
            .and_then(|s| s.current.as_ref())
            .map(|r| (r.chrom.clone(), r.ref_allele.clone(), r.format.clone()))
            .unwrap();
        let (chrom, ref_str, out_format) = ref_allele;
        let missing = missing_for(&out_format);

        // union ALTs in first-appearance order; build per-source allele maps
        let mut merged_alts: Vec<String> = Vec::new();
        let mut maps: Vec<Option<Vec<usize>>> = Vec::with_capacity(sources.len());
        let mut quals: Vec<&str> = Vec::new();
        let mut filters: Vec<&str> = Vec::new();
        for src in &sources {
            match &src.current {
                Some(r) if rank(&r.chrom) == min_rank && r.pos == min_pos => {
                    quals.push(&r.qual);
                    filters.push(&r.filter);
                    let mut m = Vec::with_capacity(r.alts.len());
                    for a in &r.alts {
                        let idx = if let Some(p) = merged_alts.iter().position(|x| x == a) {
                            p + 1
                        } else {
                            merged_alts.push(a.clone());
                            merged_alts.len()
                        };
                        m.push(idx);
                    }
                    maps.push(Some(m));
                }
                _ => maps.push(None),
            }
        }

        let qual = merge_qual(&quals);
        let filter = merge_filter(&filters);
        write!(
            out,
            "{chrom}\t{min_pos}\t.\t{ref_str}\t{}\t{qual}\t{filter}\t.\t{out_format}",
            merged_alts.join(",")
        )
        .map_err(RsomicsError::Io)?;

        for (i, src) in sources.iter().enumerate() {
            match (&maps[i], &src.current) {
                (Some(map), Some(r)) => {
                    for v in &r.sample_vals {
                        write!(out, "\t{}", remap_sample(v, map)).map_err(RsomicsError::Io)?;
                    }
                    for _ in r.sample_vals.len()..src.n_samples {
                        write!(out, "\t{missing}").map_err(RsomicsError::Io)?;
                    }
                }
                _ => {
                    for _ in 0..src.n_samples {
                        write!(out, "\t{missing}").map_err(RsomicsError::Io)?;
                    }
                }
            }
        }
        writeln!(out).map_err(RsomicsError::Io)?;
        count += 1;

        // advance every source consumed at this site
        for src in &mut sources {
            let at = src
                .current
                .as_ref()
                .is_some_and(|r| rank(&r.chrom) == min_rank && r.pos == min_pos);
            if at {
                src.advance()?;
            }
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
