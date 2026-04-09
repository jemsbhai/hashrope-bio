//! E-G1 Gene Annotation: Region queries at real gene/exon boundaries (Rust)
//!
//! Extends E-G1 by querying at biologically meaningful positions from
//! UCSC refGene annotations instead of random positions.
//!
//! Usage:
//!     cargo run --release --bin bench_region_query_g1_genes -- \
//!         --fasta ../data/chr22.fa --genes ../data/chr22_genes.tsv --output ../results/

use std::env;
use std::fs;
use std::io::{BufRead, BufReader};
use std::process::Command;
use std::time::Instant;

use hashrope::rope::{Arena, Node, NodeInner};

// ---------------------------------------------------------------------------
// Gene annotation
// ---------------------------------------------------------------------------

#[derive(Debug)]
struct Gene {
    name: String,
    tx_start: u64,
    tx_end: u64,
    tx_len: u64,
    exon_starts: Vec<u64>,
    exon_ends: Vec<u64>,
}

fn load_genes(path: &str) -> Vec<Gene> {
    let file = fs::File::open(path).unwrap_or_else(|e| panic!("Cannot open {}: {}", path, e));
    let reader = BufReader::new(file);
    let mut genes = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line.expect("Failed to read line");
        if i == 0 { continue; } // skip header
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 { continue; }

        let name = fields[0].to_string();
        let tx_start: u64 = fields[3].parse().unwrap();
        let tx_end: u64 = fields[4].parse().unwrap();
        let tx_len: u64 = fields[5].parse().unwrap();

        let exon_starts: Vec<u64> = fields[9].split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().unwrap())
            .collect();
        let exon_ends: Vec<u64> = fields[10].split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().unwrap())
            .collect();

        genes.push(Gene { name, tx_start, tx_end, tx_len, exon_starts, exon_ends });
    }
    genes
}

// ---------------------------------------------------------------------------
// FASTA loading + rope building (same as E-G1)
// ---------------------------------------------------------------------------

fn load_fasta(path: &str) -> (Vec<u8>, String) {
    let file = fs::File::open(path).unwrap_or_else(|e| panic!("Cannot open {}: {}", path, e));
    let reader = BufReader::new(file);
    let mut seq_name = String::new();
    let mut seq = Vec::new();
    let mut found_header = false;
    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        let line = line.trim_end();
        if line.starts_with('>') {
            if found_header { break; }
            seq_name = line[1..].split_whitespace().next().unwrap_or("").to_string();
            found_header = true;
            continue;
        }
        if found_header {
            for &b in line.as_bytes() {
                if b != b' ' && b != b'\t' {
                    seq.push(b.to_ascii_uppercase());
                }
            }
        }
    }
    (seq, seq_name)
}

fn build_rope_chunked(arena: &mut Arena, seq: &[u8], chunk_size: usize) -> Node {
    let mut rope: Node = None;
    for chunk in seq.chunks(chunk_size) {
        let leaf = arena.from_bytes(chunk);
        rope = arena.concat(rope, leaf);
    }
    rope
}

fn count_leaves(arena: &Arena, node: Node) -> u64 {
    match node {
        None => 0,
        Some(id) => match arena.node(id) {
            NodeInner::Leaf { .. } => 1,
            NodeInner::Internal { left, right, .. } =>
                count_leaves(arena, Some(*left)) + count_leaves(arena, Some(*right)),
            NodeInner::Repeat { child, reps, .. } =>
                count_leaves(arena, Some(*child)) * reps,
        }
    }
}

// ---------------------------------------------------------------------------
// Statistics
// ---------------------------------------------------------------------------

fn median(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n == 0 { return 0.0; }
    if n % 2 == 0 { (sorted[n/2 - 1] + sorted[n/2]) / 2.0 } else { sorted[n/2] }
}

#[allow(dead_code)]
fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() { return 0.0; }
    let idx = (p / 100.0 * (sorted.len() - 1) as f64).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

// ---------------------------------------------------------------------------
// Core benchmark
// ---------------------------------------------------------------------------

#[derive(Debug)]
#[allow(dead_code)]
struct QueryResult {
    kind: String,       // "gene_body" or "exon"
    name: String,
    start: u64,
    len: u64,
    hashrope_ns: f64,
    baseline_ns: f64,
    speedup: f64,
    hash_correct: bool,
}

fn bench_genes(
    seq: &[u8],
    arena: &mut Arena,
    rope: Node,
    genes: &[Gene],
) -> Vec<QueryResult> {
    let seq_len = seq.len() as u64;
    let mut results = Vec::new();

    for gene in genes {
        // Skip genes outside sequence bounds
        if gene.tx_end > seq_len { continue; }

        // --- Gene body query ---
        let start = gene.tx_start;
        let len = gene.tx_len;

        // hashrope
        let t0 = Instant::now();
        let h_rope = arena.substr_hash(rope, start, len);
        let hr_ns = t0.elapsed().as_nanos() as f64;

        // baseline
        let t0 = Instant::now();
        let h_base = arena.hash_bytes(&seq[start as usize..gene.tx_end as usize]);
        let bl_ns = t0.elapsed().as_nanos() as f64;

        let speedup = if hr_ns > 0.0 { bl_ns / hr_ns } else { 0.0 };

        results.push(QueryResult {
            kind: "gene_body".into(),
            name: gene.name.clone(),
            start, len,
            hashrope_ns: hr_ns,
            baseline_ns: bl_ns,
            speedup,
            hash_correct: h_rope == h_base,
        });

        // --- Exon queries ---
        for i in 0..gene.exon_starts.len() {
            let ex_start = gene.exon_starts[i];
            let ex_end = gene.exon_ends[i];
            if ex_end > seq_len { continue; }
            let ex_len = ex_end - ex_start;
            if ex_len == 0 { continue; }

            let t0 = Instant::now();
            let h_rope = arena.substr_hash(rope, ex_start, ex_len);
            let hr_ns = t0.elapsed().as_nanos() as f64;

            let t0 = Instant::now();
            let h_base = arena.hash_bytes(&seq[ex_start as usize..ex_end as usize]);
            let bl_ns = t0.elapsed().as_nanos() as f64;

            let speedup = if hr_ns > 0.0 { bl_ns / hr_ns } else { 0.0 };

            results.push(QueryResult {
                kind: "exon".into(),
                name: format!("{}:exon{}", gene.name, i + 1),
                start: ex_start,
                len: ex_len,
                hashrope_ns: hr_ns,
                baseline_ns: bl_ns,
                speedup,
                hash_correct: h_rope == h_base,
            });
        }
    }
    results
}

// ---------------------------------------------------------------------------
// Size-bin analysis
// ---------------------------------------------------------------------------

#[allow(dead_code)]
struct SizeBin {
    label: String,
    lo: u64,
    hi: u64,
    count: usize,
    hashrope_median_ns: f64,
    baseline_median_ns: f64,
    speedup_median: f64,
    mismatches: usize,
}

fn analyze_by_size(results: &[QueryResult], kind_filter: &str) -> Vec<SizeBin> {
    let bins: &[(u64, u64, &str)] = &[
        (1, 100, "<100 bp"),
        (100, 500, "100-500 bp"),
        (500, 1000, "500-1K bp"),
        (1000, 5000, "1K-5K bp"),
        (5000, 10000, "5K-10K bp"),
        (10000, 50000, "10K-50K bp"),
        (50000, 100000, "50K-100K bp"),
        (100000, 500000, "100K-500K bp"),
        (500000, 5000000, ">500K bp"),
    ];

    let mut out = Vec::new();
    for &(lo, hi, label) in bins {
        let matching: Vec<&QueryResult> = results.iter()
            .filter(|r| r.kind == kind_filter && r.len >= lo && r.len < hi)
            .collect();
        if matching.is_empty() { continue; }

        let mut hr_times: Vec<f64> = matching.iter().map(|r| r.hashrope_ns).collect();
        let mut bl_times: Vec<f64> = matching.iter().map(|r| r.baseline_ns).collect();
        hr_times.sort_by(|a, b| a.partial_cmp(b).unwrap());
        bl_times.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let hr_med = median(&hr_times);
        let bl_med = median(&bl_times);
        let mismatches = matching.iter().filter(|r| !r.hash_correct).count();

        out.push(SizeBin {
            label: label.into(),
            lo, hi,
            count: matching.len(),
            hashrope_median_ns: hr_med,
            baseline_median_ns: bl_med,
            speedup_median: if hr_med > 0.0 { bl_med / hr_med } else { 0.0 },
            mismatches,
        });
    }
    out
}

// ---------------------------------------------------------------------------
// JSON output
// ---------------------------------------------------------------------------

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"").replace('\n', "\\n")
}

fn get_environment() -> Vec<(String, String)> {
    let mut env_data = Vec::new();
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH).unwrap().as_secs();
    env_data.push(("timestamp_unix".into(), now.to_string()));
    env_data.push(("os".into(), env::consts::OS.into()));
    env_data.push(("arch".into(), env::consts::ARCH.into()));
    if let Ok(output) = Command::new("wmic").args(["cpu", "get", "name"]).output() {
        let s = String::from_utf8_lossy(&output.stdout);
        for line in s.lines() {
            let line = line.trim();
            if !line.is_empty() && line != "Name" {
                env_data.push(("cpu_name".into(), line.into()));
                break;
            }
        }
    }
    if let Ok(output) = Command::new("rustc").arg("--version").output() {
        env_data.push(("rustc_version".into(), String::from_utf8_lossy(&output.stdout).trim().to_string()));
    }
    env_data.push(("language".into(), "rust".into()));
    env_data.push(("build_profile".into(), if cfg!(debug_assertions) { "debug" } else { "release" }.into()));
    env_data
}

fn days_to_ymd(mut days: u64) -> (u64, u64, u64) {
    let mut year = 1970;
    loop {
        let diy = if year % 4 == 0 && (year % 100 != 0 || year % 400 == 0) { 366 } else { 365 };
        if days < diy { break; }
        days -= diy;
        year += 1;
    }
    let md: &[u64] = if year % 4 == 0 && (year % 100 != 0 || year % 400 == 0) {
        &[31,29,31,30,31,30,31,31,30,31,30,31]
    } else {
        &[31,28,31,30,31,30,31,31,30,31,30,31]
    };
    let mut month = 1;
    for &d in md { if days < d { break; } days -= d; month += 1; }
    (year, month, days + 1)
}

fn write_json(
    path: &str,
    env_data: &[(String, String)],
    gene_bins: &[SizeBin],
    exon_bins: &[SizeBin],
    total_queries: usize,
    total_genes: usize,
    total_exons: usize,
    total_mismatches: usize,
    chunk_size: usize,
    seq_len: u64,
) {
    let mut o = String::new();
    o.push_str("{\n");
    o.push_str("  \"benchmark\": \"region_query_genes_rust\",\n");

    // Environment
    o.push_str("  \"environment\": {\n");
    for (i, (k, v)) in env_data.iter().enumerate() {
        let c = if i + 1 < env_data.len() { "," } else { "" };
        o.push_str(&format!("    \"{}\": \"{}\"{}\n", escape_json(k), escape_json(v), c));
    }
    o.push_str("  },\n");

    // Config
    o.push_str(&format!("  \"config\": {{\n\
        \x20   \"chunk_size\": {},\n\
        \x20   \"seq_len\": {},\n\
        \x20   \"genes_queried\": {},\n\
        \x20   \"exons_queried\": {},\n\
        \x20   \"total_queries\": {},\n\
        \x20   \"total_hash_mismatches\": {},\n\
        \x20   \"data_source\": \"UCSC refGene hg38, chr22\"\n\
        \x20 }},\n", chunk_size, seq_len, total_genes, total_exons, total_queries, total_mismatches));

    // Gene body bins
    o.push_str("  \"gene_body_results\": [\n");
    for (i, b) in gene_bins.iter().enumerate() {
        let c = if i + 1 < gene_bins.len() { "," } else { "" };
        o.push_str(&format!(
            "    {{\"size_bin\": \"{}\", \"count\": {}, \"hashrope_median_ns\": {:.1}, \"baseline_median_ns\": {:.1}, \"speedup\": {:.2}, \"mismatches\": {}}}{}\n",
            b.label, b.count, b.hashrope_median_ns, b.baseline_median_ns, b.speedup_median, b.mismatches, c));
    }
    o.push_str("  ],\n");

    // Exon bins
    o.push_str("  \"exon_results\": [\n");
    for (i, b) in exon_bins.iter().enumerate() {
        let c = if i + 1 < exon_bins.len() { "," } else { "" };
        o.push_str(&format!(
            "    {{\"size_bin\": \"{}\", \"count\": {}, \"hashrope_median_ns\": {:.1}, \"baseline_median_ns\": {:.1}, \"speedup\": {:.2}, \"mismatches\": {}}}{}\n",
            b.label, b.count, b.hashrope_median_ns, b.baseline_median_ns, b.speedup_median, b.mismatches, c));
    }
    o.push_str("  ]\n");

    o.push_str("}\n");
    fs::write(path, &o).unwrap_or_else(|e| panic!("Cannot write {}: {}", path, e));
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut fasta_path = String::new();
    let mut genes_path = String::new();
    let mut output_dir: Option<String> = None;
    let mut chunk_size: usize = 256;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--fasta" => { i += 1; fasta_path = args[i].clone(); }
            "--genes" => { i += 1; genes_path = args[i].clone(); }
            "--output" => { i += 1; output_dir = Some(args[i].clone()); }
            "--chunk-size" => { i += 1; chunk_size = args[i].parse().expect("Invalid chunk size"); }
            other => panic!("Unknown argument: {}", other),
        }
        i += 1;
    }

    if fasta_path.is_empty() || genes_path.is_empty() {
        eprintln!("Usage: bench_region_query_g1_genes --fasta <path> --genes <tsv> [--output <dir>] [--chunk-size 256]");
        std::process::exit(1);
    }

    eprintln!("hashrope-bio E-G1 Gene Annotation (Rust)");
    eprintln!("Build profile: {}", if cfg!(debug_assertions) { "debug" } else { "release" });

    // Load sequence
    eprintln!("\nLoading sequence...");
    let t0 = Instant::now();
    let (seq, seq_name) = load_fasta(&fasta_path);
    eprintln!("  {} — {} bp in {:.3}s", seq_name, seq.len(), t0.elapsed().as_secs_f64());

    // Load gene annotations
    eprintln!("Loading gene annotations...");
    let genes = load_genes(&genes_path);
    let total_exons_count: usize = genes.iter().map(|g| g.exon_starts.len()).sum();
    eprintln!("  {} genes, {} exons", genes.len(), total_exons_count);

    // Build rope
    eprintln!("Building rope (chunk_size={})...", chunk_size);
    let mut arena = Arena::new();
    let rope = build_rope_chunked(&mut arena, &seq, chunk_size);
    let leaves = count_leaves(&arena, rope);
    let height = arena.height(rope);
    eprintln!("  {} leaves, height {}", leaves, height);

    // Run benchmark
    eprintln!("\nRunning gene/exon queries...");
    let t0 = Instant::now();
    let results = bench_genes(&seq, &mut arena, rope, &genes);
    let wall = t0.elapsed().as_secs_f64();
    eprintln!("  {} queries in {:.3}s", results.len(), wall);

    let mismatches = results.iter().filter(|r| !r.hash_correct).count();
    if mismatches > 0 {
        eprintln!("  *** {} HASH MISMATCHES! ***", mismatches);
    } else {
        eprintln!("  Hash verification: {}/{} correct (100%)", results.len(), results.len());
    }

    // Analyze by size
    let gene_count = results.iter().filter(|r| r.kind == "gene_body").count();
    let exon_count = results.iter().filter(|r| r.kind == "exon").count();

    let gene_bins = analyze_by_size(&results, "gene_body");
    let exon_bins = analyze_by_size(&results, "exon");

    // Print summary
    eprintln!("\n{}", "=".repeat(90));
    eprintln!("Gene Body Queries ({} genes)", gene_count);
    eprintln!("{}", "=".repeat(90));
    eprintln!("{:<16} {:>6} {:>12} {:>12} {:>10}", "Size bin", "Count", "rope_ns", "base_ns", "Speedup");
    eprintln!("{}", "-".repeat(90));
    for b in &gene_bins {
        eprintln!("{:<16} {:>6} {:>12.0} {:>12.0} {:>9.1}×", b.label, b.count, b.hashrope_median_ns, b.baseline_median_ns, b.speedup_median);
    }

    eprintln!("\n{}", "=".repeat(90));
    eprintln!("Exon Queries ({} exons)", exon_count);
    eprintln!("{}", "=".repeat(90));
    eprintln!("{:<16} {:>6} {:>12} {:>12} {:>10}", "Size bin", "Count", "rope_ns", "base_ns", "Speedup");
    eprintln!("{}", "-".repeat(90));
    for b in &exon_bins {
        eprintln!("{:<16} {:>6} {:>12.0} {:>12.0} {:>9.1}×", b.label, b.count, b.hashrope_median_ns, b.baseline_median_ns, b.speedup_median);
    }

    // Top 10 genes by speedup
    let mut gene_results: Vec<&QueryResult> = results.iter().filter(|r| r.kind == "gene_body").collect();
    gene_results.sort_by(|a, b| b.speedup.partial_cmp(&a.speedup).unwrap());
    eprintln!("\n--- Top 10 Genes by Speedup ---");
    for r in gene_results.iter().take(10) {
        eprintln!("  {:>20} {:>9} bp  rope={:>8.0}ns  base={:>10.0}ns  {:>6.1}×",
                  r.name, r.len, r.hashrope_ns, r.baseline_ns, r.speedup);
    }

    // Save JSON
    if let Some(ref dir) = output_dir {
        fs::create_dir_all(dir).ok();
        let env_data = get_environment();
        let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs();
        let days = now / 86400;
        let tod = now % 86400;
        let (y, m, d) = days_to_ymd(days);
        let ts = format!("{:04}{:02}{:02}T{:02}{:02}{:02}Z", y, m, d, tod/3600, (tod%3600)/60, tod%60);

        let latest = format!("{}/region_query_genes_rust.json", dir);
        let archive = format!("{}/region_query_genes_rust_{}.json", dir, ts);

        write_json(&latest, &env_data, &gene_bins, &exon_bins,
                    results.len(), gene_count, exon_count, mismatches,
                    chunk_size, seq.len() as u64);
        write_json(&archive, &env_data, &gene_bins, &exon_bins,
                    results.len(), gene_count, exon_count, mismatches,
                    chunk_size, seq.len() as u64);

        eprintln!("\n  Results saved:");
        eprintln!("    Latest:  {}", latest);
        eprintln!("    Archive: {}", archive);
    }
}
