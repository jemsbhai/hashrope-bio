//! E-G1: Genome-Scale substr_hash Query Performance (Rust)
//!
//! Measures region identity query time at various region sizes L on a
//! pre-built genome rope. Compares substr_hash (O(log w)) vs baseline
//! slice+hash (O(L)).
//!
//! Usage:
//!     cargo run --release --bin bench_region_query_g1 -- --fasta ../data/chr22.fa --output ../../../results/

use std::env;
use std::fs;
use std::io::{BufRead, BufReader};
use std::process::Command;
use std::time::Instant;

use hashrope::rope::{Arena, Node, NodeInner};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

const QUERIES_PER_SIZE: usize = 1000;
const RNG_SEED: u64 = 42;

// ---------------------------------------------------------------------------
// Simple LCG RNG (deterministic, no deps)
// ---------------------------------------------------------------------------

struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self { Self(seed) }
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn next_range(&mut self, max: u64) -> u64 {
        self.next() % max
    }
}

// ---------------------------------------------------------------------------
// FASTA loading (same as E-G4)
// ---------------------------------------------------------------------------

fn load_fasta(path: &str) -> (Vec<u8>, String) {
    let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("Cannot open {}: {}", path, e));
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
    if seq.is_empty() { panic!("No sequence data in {}", path); }
    (seq, seq_name)
}

// ---------------------------------------------------------------------------
// Rope building
// ---------------------------------------------------------------------------

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
        Some(id) => count_leaves_inner(arena, id),
    }
}

fn count_leaves_inner(arena: &Arena, id: u32) -> u64 {
    match arena.node(id) {
        NodeInner::Leaf { .. } => 1,
        NodeInner::Internal { left, right, .. } => {
            count_leaves_inner(arena, *left) + count_leaves_inner(arena, *right)
        }
        NodeInner::Repeat { child, reps, .. } => count_leaves_inner(arena, *child) * reps,
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

fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() { return 0.0; }
    let idx = (p / 100.0 * (sorted.len() - 1) as f64).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

// ---------------------------------------------------------------------------
// Core benchmark
// ---------------------------------------------------------------------------

struct RegionResult {
    chunk_size: usize,
    region_len: u64,
    leaves: u64,
    height: u64,
    // hashrope
    hashrope_median_ns: f64,
    hashrope_p5_ns: f64,
    hashrope_p95_ns: f64,
    // baseline
    baseline_median_ns: f64,
    baseline_p5_ns: f64,
    baseline_p95_ns: f64,
    // derived
    speedup_median: f64,
    queries: usize,
}

fn bench_region_queries(
    seq: &[u8],
    chunk_size: usize,
    region_sizes: &[u64],
) -> Vec<RegionResult> {
    let mut results = Vec::new();
    let seq_len = seq.len() as u64;

    // Build rope once per chunk size
    let mut arena = Arena::new();
    let rope = build_rope_chunked(&mut arena, seq, chunk_size);
    let height = arena.height(rope);
    let leaves = count_leaves(&arena, rope);

    eprintln!("  Built rope: chunk_size={}, leaves={}, height={}", chunk_size, leaves, height);

    for &l in region_sizes {
        if l > seq_len {
            eprintln!("    L={}: skipped (exceeds seq_len={})", l, seq_len);
            continue;
        }

        let max_start = seq_len - l;
        let mut rng = Lcg::new(RNG_SEED);

        // Generate random start positions
        let starts: Vec<u64> = (0..QUERIES_PER_SIZE).map(|_| rng.next_range(max_start)).collect();

        // --- hashrope: substr_hash ---
        // Warmup
        for &s in starts.iter().take(50) {
            std::hint::black_box(arena.substr_hash(rope, s, l));
        }

        let mut hashrope_times = Vec::with_capacity(QUERIES_PER_SIZE);
        for &s in &starts {
            let t0 = Instant::now();
            std::hint::black_box(arena.substr_hash(rope, s, l));
            let elapsed = t0.elapsed().as_nanos() as f64;
            hashrope_times.push(elapsed);
        }
        hashrope_times.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // --- baseline: slice + hash ---
        // Warmup
        for &s in starts.iter().take(50) {
            let slice = &seq[s as usize..(s + l) as usize];
            std::hint::black_box(arena.hash_bytes(slice));
        }

        let mut baseline_times = Vec::with_capacity(QUERIES_PER_SIZE);
        for &s in &starts {
            let slice = &seq[s as usize..(s + l) as usize];
            let t0 = Instant::now();
            std::hint::black_box(arena.hash_bytes(slice));
            let elapsed = t0.elapsed().as_nanos() as f64;
            baseline_times.push(elapsed);
        }
        baseline_times.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Verify correctness on a subset
        let mut mismatches = 0;
        for &s in starts.iter().take(100) {
            let h_rope = arena.substr_hash(rope, s, l);
            let h_base = arena.hash_bytes(&seq[s as usize..(s + l) as usize]);
            if h_rope != h_base { mismatches += 1; }
        }
        if mismatches > 0 {
            eprintln!("    L={}: {} HASH MISMATCHES!", l, mismatches);
        }

        let hr_med = median(&hashrope_times);
        let bl_med = median(&baseline_times);
        let speedup = if hr_med > 0.0 { bl_med / hr_med } else { 0.0 };

        eprintln!("    L={:>9}: rope={:>9.0}ns  base={:>9.0}ns  speedup={:>6.1}×  (p5={:.0}/{:.0}, p95={:.0}/{:.0})",
                  l, hr_med, bl_med, speedup,
                  percentile(&hashrope_times, 5.0), percentile(&baseline_times, 5.0),
                  percentile(&hashrope_times, 95.0), percentile(&baseline_times, 95.0));

        results.push(RegionResult {
            chunk_size,
            region_len: l,
            leaves,
            height,
            hashrope_median_ns: hr_med,
            hashrope_p5_ns: percentile(&hashrope_times, 5.0),
            hashrope_p95_ns: percentile(&hashrope_times, 95.0),
            baseline_median_ns: bl_med,
            baseline_p5_ns: percentile(&baseline_times, 5.0),
            baseline_p95_ns: percentile(&baseline_times, 95.0),
            speedup_median: speedup,
            queries: QUERIES_PER_SIZE,
        });
    }

    results
}

// ---------------------------------------------------------------------------
// JSON output
// ---------------------------------------------------------------------------

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"").replace('\n', "\\n")
}

fn write_json(
    path: &str,
    env_data: &[(String, String)],
    all_results: &[RegionResult],
    config: &[(String, String)],
) {
    let mut out = String::new();
    out.push_str("{\n");
    out.push_str("  \"benchmark\": \"region_query_rust\",\n");

    // Environment
    out.push_str("  \"environment\": {\n");
    for (i, (k, v)) in env_data.iter().enumerate() {
        let comma = if i + 1 < env_data.len() { "," } else { "" };
        out.push_str(&format!("    \"{}\": \"{}\"{}\n", escape_json(k), escape_json(v), comma));
    }
    out.push_str("  },\n");

    // Config
    out.push_str("  \"config\": {\n");
    for (i, (k, v)) in config.iter().enumerate() {
        let comma = if i + 1 < config.len() { "," } else { "" };
        if v.parse::<f64>().is_ok() || v == "true" || v == "false" || v.starts_with('[') {
            out.push_str(&format!("    \"{}\": {}{}\n", escape_json(k), v, comma));
        } else {
            out.push_str(&format!("    \"{}\": \"{}\"{}\n", escape_json(k), escape_json(v), comma));
        }
    }
    out.push_str("  },\n");

    // Results
    out.push_str("  \"results\": [\n");
    for (i, r) in all_results.iter().enumerate() {
        let comma = if i + 1 < all_results.len() { "," } else { "" };
        out.push_str(&format!(
            "    {{\n\
             \x20     \"chunk_size\": {},\n\
             \x20     \"region_len\": {},\n\
             \x20     \"leaves\": {},\n\
             \x20     \"height\": {},\n\
             \x20     \"hashrope_median_ns\": {:.1},\n\
             \x20     \"hashrope_p5_ns\": {:.1},\n\
             \x20     \"hashrope_p95_ns\": {:.1},\n\
             \x20     \"baseline_median_ns\": {:.1},\n\
             \x20     \"baseline_p5_ns\": {:.1},\n\
             \x20     \"baseline_p95_ns\": {:.1},\n\
             \x20     \"speedup_median\": {:.2},\n\
             \x20     \"queries\": {}\n\
             \x20   }}{}\n",
            r.chunk_size, r.region_len, r.leaves, r.height,
            r.hashrope_median_ns, r.hashrope_p5_ns, r.hashrope_p95_ns,
            r.baseline_median_ns, r.baseline_p5_ns, r.baseline_p95_ns,
            r.speedup_median, r.queries, comma
        ));
    }
    out.push_str("  ]\n");
    out.push_str("}\n");

    fs::write(path, &out).unwrap_or_else(|e| panic!("Cannot write {}: {}", path, e));
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
    if let Ok(output) = Command::new("wmic").args(["computersystem", "get", "totalphysicalmemory"]).output() {
        let s = String::from_utf8_lossy(&output.stdout);
        for line in s.lines() {
            let line = line.trim();
            if !line.is_empty() && line != "TotalPhysicalMemory" {
                if let Ok(bytes) = line.parse::<u64>() {
                    env_data.push(("ram_gb".into(), format!("{:.1}", bytes as f64 / (1024.0*1024.0*1024.0))));
                }
                break;
            }
        }
    }
    if let Ok(output) = Command::new("rustc").arg("--version").output() {
        env_data.push(("rustc_version".into(), String::from_utf8_lossy(&output.stdout).trim().to_string()));
    }
    env_data.push(("hashrope_crate_version".into(), "0.2.1".into()));
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

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut fasta_path = String::new();
    let mut output_dir: Option<String> = None;
    let mut chunk_sizes_str: Option<String> = None;
    let mut region_sizes_str: Option<String> = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--fasta" => { i += 1; fasta_path = args[i].clone(); }
            "--output" => { i += 1; output_dir = Some(args[i].clone()); }
            "--chunk-sizes" => { i += 1; chunk_sizes_str = Some(args[i].clone()); }
            "--region-sizes" => { i += 1; region_sizes_str = Some(args[i].clone()); }
            other => panic!("Unknown argument: {}", other),
        }
        i += 1;
    }

    if fasta_path.is_empty() {
        eprintln!("Usage: bench_region_query_g1 --fasta <path> [--output <dir>] [--chunk-sizes 256,1024,4096] [--region-sizes 100,500,...]");
        std::process::exit(1);
    }

    let chunk_sizes: Vec<usize> = chunk_sizes_str
        .as_deref()
        .unwrap_or("256,1024,4096")
        .split(',')
        .map(|s| s.trim().parse().expect("Invalid chunk size"))
        .collect();

    let region_sizes: Vec<u64> = region_sizes_str
        .as_deref()
        .unwrap_or("100,500,1000,5000,10000,50000,100000,500000,1000000")
        .split(',')
        .map(|s| s.trim().parse().expect("Invalid region size"))
        .collect();

    eprintln!("hashrope-bio E-G1 (Rust): Region Query Scaling");
    eprintln!("Build profile: {}", if cfg!(debug_assertions) { "debug" } else { "release" });
    eprintln!("FASTA: {}", fasta_path);
    eprintln!("Chunk sizes: {:?}", chunk_sizes);
    eprintln!("Region sizes: {:?}", region_sizes);
    eprintln!("Queries per size: {}", QUERIES_PER_SIZE);
    eprintln!("RNG seed: {}", RNG_SEED);

    // Load sequence
    eprintln!("\nLoading sequence ...");
    let t0 = Instant::now();
    let (seq, seq_name) = load_fasta(&fasta_path);
    eprintln!("  Loaded: {}, {} bp ({:.1} MB) in {:.3}s",
              seq_name, seq.len(), seq.len() as f64 / (1024.0*1024.0), t0.elapsed().as_secs_f64());

    // Run benchmarks for each chunk size
    let mut all_results = Vec::new();

    for &cs in &chunk_sizes {
        eprintln!("\n=== chunk_size={} ===", cs);
        let results = bench_region_queries(&seq, cs, &region_sizes);
        all_results.extend(results);
    }

    // Summary table
    eprintln!("\n{}", "=".repeat(100));
    eprintln!("E-G1: Region Query Scaling Summary");
    eprintln!("{}", "=".repeat(100));
    eprintln!("{:>6} {:>9} {:>10} {:>10} {:>8} {:>14} {:>14}",
              "chunk", "L (bp)", "rope_ns", "base_ns", "speedup", "rope_p5-p95", "base_p5-p95");
    eprintln!("{}", "-".repeat(100));

    for r in &all_results {
        eprintln!("{:>6} {:>9} {:>10.0} {:>10.0} {:>7.1}× {:>6.0}-{:<6.0} {:>6.0}-{:<6.0}",
                  r.chunk_size, r.region_len,
                  r.hashrope_median_ns, r.baseline_median_ns,
                  r.speedup_median,
                  r.hashrope_p5_ns, r.hashrope_p95_ns,
                  r.baseline_p5_ns, r.baseline_p95_ns);
    }

    // Hash cross-validation
    eprintln!("\n--- Hash Cross-Validation ---");
    {
        let arena = Arena::new();
        eprintln!("  Full sequence hash: {}", arena.hash_bytes(&seq));
    }

    // Save JSON
    if let Some(ref dir) = output_dir {
        fs::create_dir_all(dir).ok();
        let env_data = get_environment();
        let config: Vec<(String, String)> = vec![
            ("fasta".into(), fasta_path.clone()),
            ("seq_name".into(), seq_name.clone()),
            ("seq_len".into(), seq.len().to_string()),
            ("chunk_sizes".into(), format!("[{}]", chunk_sizes.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(", "))),
            ("region_sizes".into(), format!("[{}]", region_sizes.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(", "))),
            ("queries_per_size".into(), QUERIES_PER_SIZE.to_string()),
            ("rng_seed".into(), RNG_SEED.to_string()),
            ("hash_base".into(), "131".into()),
            ("hash_prime".into(), "2305843009213693951".into()),
        ];

        let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs();
        let secs_in_day = 86400u64;
        let days = now / secs_in_day;
        let tod = now % secs_in_day;
        let (y, m, d) = days_to_ymd(days);
        let ts = format!("{:04}{:02}{:02}T{:02}{:02}{:02}Z", y, m, d, tod/3600, (tod%3600)/60, tod%60);

        let latest = format!("{}/region_query_rust.json", dir);
        let archive = format!("{}/region_query_rust_{}.json", dir, ts);

        write_json(&latest, &env_data, &all_results, &config);
        write_json(&archive, &env_data, &all_results, &config);

        eprintln!("\n  Results saved:");
        eprintln!("    Latest:  {}", latest);
        eprintln!("    Archive: {}", archive);
    }
}
