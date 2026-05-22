use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

struct Site {
    // CHROM POS ID REF ALT QUAL FILTER INFO from the first input carrying the site
    template: [String; 8],
    // file index -> (FORMAT, that file's sample value columns)
    per_file: HashMap<usize, (String, Vec<String>)>,
}

/// Merge single- or multi-sample VCFs into one multi-sample VCF, matching
/// `bcftools merge`: union of sites, every input's samples as columns, missing
/// genotypes filled with `./.` (other FORMAT keys with `.`). Output is ordered
/// by the `##contig` header order then position.
pub fn merge(inputs: &[&Path], output: &mut dyn Write) -> Result<u64> {
    let mut out = BufWriter::with_capacity(64 * 1024, output);
    let mut header_meta: Vec<String> = Vec::new();
    let mut contig_order: HashMap<String, usize> = HashMap::new();
    let mut all_samples: Vec<String> = Vec::new();
    let mut file_sample_counts: Vec<usize> = Vec::new();
    let mut sites: HashMap<(String, u64, String, String), Site> = HashMap::new();
    let mut order: Vec<(String, u64, String, String)> = Vec::new();

    for (fi, input) in inputs.iter().enumerate() {
        let file = File::open(input)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
        let mut nsamp = 0;
        for line in BufReader::new(file).lines() {
            let line = line.map_err(RsomicsError::Io)?;
            if let Some(rest) = line.strip_prefix("##") {
                if fi == 0 {
                    header_meta.push(line.clone());
                    if let Some(id) = rest.strip_prefix("contig=<ID=") {
                        let name = id.split([',', '>']).next().unwrap_or("");
                        if !name.is_empty() {
                            let next = contig_order.len();
                            contig_order.entry(name.to_string()).or_insert(next);
                        }
                    }
                }
                continue;
            }
            if line.starts_with('#') {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() > 9 {
                    for s in &f[9..] {
                        all_samples.push((*s).to_string());
                    }
                    nsamp = f.len() - 9;
                }
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                continue;
            }
            let key = (
                f[0].to_string(),
                f[1].parse::<u64>().unwrap_or(0),
                f[3].to_string(),
                f[4].to_string(),
            );
            let format = f.get(8).copied().unwrap_or("GT").to_string();
            let values: Vec<String> = if f.len() > 9 {
                f[9..].iter().map(|s| (*s).to_string()).collect()
            } else {
                Vec::new()
            };
            let g = |i: usize| f.get(i).copied().unwrap_or(".").to_string();
            let site = sites.entry(key.clone()).or_insert_with(|| {
                order.push(key.clone());
                Site {
                    template: [
                        g(0),
                        g(1),
                        g(2),
                        f[3].to_string(),
                        f[4].to_string(),
                        g(5),
                        g(6),
                        g(7),
                    ],
                    per_file: HashMap::new(),
                }
            });
            site.per_file.insert(fi, (format, values));
        }
        file_sample_counts.push(nsamp);
    }

    let contig_rank = |chrom: &str| contig_order.get(chrom).copied().unwrap_or(usize::MAX);
    order.sort_by(|a, b| {
        contig_rank(&a.0)
            .cmp(&contig_rank(&b.0))
            .then(a.1.cmp(&b.1))
            .then(a.2.cmp(&b.2))
            .then(a.3.cmp(&b.3))
    });

    for h in &header_meta {
        writeln!(out, "{h}").map_err(RsomicsError::Io)?;
    }
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        .map_err(RsomicsError::Io)?;
    for s in &all_samples {
        write!(out, "\t{s}").map_err(RsomicsError::Io)?;
    }
    writeln!(out).map_err(RsomicsError::Io)?;

    let mut count: u64 = 0;
    for key in &order {
        let site = &sites[key];
        let mut present: Vec<usize> = site.per_file.keys().copied().collect();
        present.sort_unstable();
        let out_format = &site.per_file[&present[0]].0;
        let missing: String = out_format
            .split(':')
            .map(|k| if k == "GT" { "./." } else { "." })
            .collect::<Vec<_>>()
            .join(":");

        write!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            site.template[0],
            site.template[1],
            site.template[2],
            site.template[3],
            site.template[4],
            site.template[5],
            site.template[6],
            site.template[7],
            out_format
        )
        .map_err(RsomicsError::Io)?;
        for (fi, &cnt) in file_sample_counts.iter().enumerate() {
            if let Some((_fmt, vals)) = site.per_file.get(&fi) {
                for v in vals {
                    write!(out, "\t{v}").map_err(RsomicsError::Io)?;
                }
                for _ in vals.len()..cnt {
                    write!(out, "\t{missing}").map_err(RsomicsError::Io)?;
                }
            } else {
                for _ in 0..cnt {
                    write!(out, "\t{missing}").map_err(RsomicsError::Io)?;
                }
            }
        }
        writeln!(out).map_err(RsomicsError::Io)?;
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
