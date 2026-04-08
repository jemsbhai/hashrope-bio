//! E-G3: Mutation Localization via Binary Search (Rust)
//!
//! Introduces synthetic SNPs into chr22, then uses binary search via
//! substr_hash to pinpoint mutation positions. Measures comparison count
//! and wall-clock time vs linear scan baseline.
//!
//! Note: chr22 contains large N-runs (centromere/telomere). Mutations are
//! placed only on valid nucleotides (A/C/G/T) to ensure a real difference.
//! Regions that are entirely N are skipped and re-sampled.
//!
//! Usage:
//!     cargo run --release --bin bench_mutation_g3 -- --fasta ../data/chr22.fa --output ../../../results/

use std::env;
use std::fs;
use std::io::{BufRead, BufReader};
use std::process::Command;
use std::time::Instant;

use hashrope::rope::{Arena, Node};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

const REGION_SIZES: &[u64] = &[1_000, 10_000, 100_000, 1_000_000, 10_000_000, 50_818_468];
const TRIALS_PER_SIZE: usize = 100;
const CHUNK_SIZE: usize = 256; // Best query performance from E-G1
const RNG_SEED: u64 = 42;

// ---------------------------------------------------------------------------
// LCG RNG
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
// FASTA loading
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

// ---------------------------------------------------------------------------
// Nucleotide helpers
// ---------------------------------------------------------------------------

fn is_valid_nucleotide(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

fn flip_nucleotide(b: u8) -> u8 {
    match b {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        b'T' => b'A',
        _ => panic!("flip_nucleotide called on non-ACGT byte: {}", b),
    }
}

/// Find a valid region start AND mutation position.
/// Retries both region and position until a region with valid nucleotides is found.
/// Returns (region_start, mut_pos_absolute).
fn find_valid_region_and_mutation(
    seq: &[u8], seq_len: u64, n: u64, rng: &mut Lcg,
) -> (u64, usize) {
    for _ in 0..10_000 {
        let region_start = if n < seq_len { rng.next_range(seq_len - n) } else { 0 };
        // Try to find a valid nucleotide within this region
        for _ in 0..100 {
            let offset = rng.next_range(n);
            let pos = (region_start + offset) as usize;
            if is_valid_nucleotide(seq[pos]) {
                return (region_start, pos);
            }
        }
        // Linear scan fallback for this region
        let start = region_start as usize;
        for i in start..start + n as usize {
            if is_valid_nucleotide(seq[i]) {
                return (region_start, i);
            }
        }
        // Region is all-N — try a different region
    }
    panic!("Could not find any valid region after 10,000 attempts");
}

// ---------------------------------------------------------------------------
// Binary search localization
// ---------------------------------------------------------------------------

fn binary_search_mutation(
    arena: &mut Arena,
    ref_rope: Node,
    sample_rope: Node,
    start: u64,
    length: u64,
) -> (u64, u64) {
    // Returns (found_position, comparison_count)
    let mut lo = start;
    let mut hi = start + length;
    let mut comparisons = 0u64;

    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        let h_ref = arena.substr_hash(ref_rope, lo, mid - lo);
        let h_sam = arena.substr_hash(sample_rope, lo, mid - lo);
        comparisons += 1;
        if h_ref != h_sam {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    (lo, comparisons)
}

// ---------------------------------------------------------------------------
// Linear scan baseline
// ---------------------------------------------------------------------------

fn linear_scan(ref_bytes: &[u8], sample_bytes: &[u8], start: usize, length: usize) -> Option<usize> {
    for i in start..start + length {
        if ref_bytes[i] != sample_bytes[i] {
            return Some(i);
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Core benchmark
// ---------------------------------------------------------------------------

struct MutationResult {
    region_size: u64,
    trials: usize,
    expected_comparisons: u64,
    median_comparisons: u64,
    median_search_ns: f64,
    p5_search_ns: f64,
    p95_search_ns: f64,
    median_linear_ns: f64,
    p5_linear_ns: f64,
    p95_linear_ns: f64,
    speedup_median: f64,
    all_correct: bool,
}

fn median_u64(sorted: &[u64]) -> u64 {
    let n = sorted.len();
    if n % 2 == 0 { (sorted[n/2 - 1] + sorted[n/2]) / 2 } else { sorted[n/2] }
}

fn median_f64(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n % 2 == 0 { (sorted[n/2 - 1] + sorted[n/2]) / 2.0 } else { sorted[n/2] }
}

fn percentile_f64(sorted: &[f64], p: f64) -> f64 {
    let idx = (p / 100.0 * (sorted.len() - 1) as f64).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

fn ceil_log2(n: u64) -> u64 {
    if n <= 1 { return 0; }
    64 - (n - 1).leading_zeros() as u64
}

fn bench_mutation(seq: &[u8]) -> Vec<MutationResult> {
    let mut results = Vec::new();
    let seq_len = seq.len() as u64;

    for &n in REGION_SIZES {
        let n = n.min(seq_len);
        let expected_comps = ceil_log2(n);
        let actual_trials = if n == seq_len { TRIALS_PER_SIZE.min(20) } else { TRIALS_PER_SIZE };

        eprintln!("  N={:>12}: expected ⌈log₂(N)⌉={} comparisons, {} trials ...",
                  n, expected_comps, actual_trials);

        let mut rng = Lcg::new(RNG_SEED);
        let mut comp_counts = Vec::with_capacity(actual_trials);
        let mut search_times = Vec::with_capacity(actual_trials);
        let mut linear_times = Vec::with_capacity(actual_trials);
        let mut all_correct = true;

        for _ in 0..actual_trials {
            // Pick a random region with at least one valid nucleotide, and a mutation position
            let (region_start, mut_pos) = find_valid_region_and_mutation(seq, seq_len, n, &mut rng);

            // Create mutated sequence — assert the flip actually changes the byte
            let mut sample = seq.to_vec();
            let original = seq[mut_pos];
            let flipped = flip_nucleotide(original);
            assert_ne!(original, flipped, "flip_nucleotide must change the byte");
            sample[mut_pos] = flipped;

            // Build ropes (fresh arena per trial to avoid arena bloat)
            let mut arena = Arena::new();
            let ref_rope = build_rope_chunked(&mut arena, &seq[region_start as usize..(region_start + n) as usize], CHUNK_SIZE);
            let sample_rope = build_rope_chunked(&mut arena, &sample[region_start as usize..(region_start + n) as usize], CHUNK_SIZE);

            // --- Binary search ---
            let t0 = Instant::now();
            let (found, comps) = binary_search_mutation(&mut arena, ref_rope, sample_rope, 0, n);
            let search_ns = t0.elapsed().as_nanos() as f64;

            let found_abs = region_start + found;
            let correct = found_abs as usize == mut_pos;
            if !correct {
                eprintln!("    MISMATCH: expected={}, found={}, region_start={}, relative_expected={}, relative_found={}",
                          mut_pos, found_abs, region_start, mut_pos - region_start as usize, found);
                all_correct = false;
            }

            comp_counts.push(comps);
            search_times.push(search_ns);

            // --- Linear scan baseline ---
            let t0 = Instant::now();
            let lin_pos = linear_scan(&seq, &sample, region_start as usize, n as usize);
            let linear_ns = t0.elapsed().as_nanos() as f64;

            // Verify linear scan agrees
            if let Some(lp) = lin_pos {
                if lp != mut_pos {
                    eprintln!("    LINEAR MISMATCH: expected={}, found={}", mut_pos, lp);
                    all_correct = false;
                }
            } else {
                eprintln!("    LINEAR SCAN: no difference found (impossible — mutation was verified)");
                all_correct = false;
            }

            linear_times.push(linear_ns);
        }

        comp_counts.sort();
        search_times.sort_by(|a, b| a.partial_cmp(b).unwrap());
        linear_times.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let med_search = median_f64(&search_times);
        let med_linear = median_f64(&linear_times);
        let speedup = if med_search > 0.0 { med_linear / med_search } else { 0.0 };

        eprintln!("    comps={}, search={:.0}ns, linear={:.0}ns, speedup={:.0}×, correct={}",
                  median_u64(&comp_counts), med_search, med_linear, speedup, all_correct);

        results.push(MutationResult {
            region_size: n,
            trials: actual_trials,
            expected_comparisons: expected_comps,
            median_comparisons: median_u64(&comp_counts),
            median_search_ns: med_search,
            p5_search_ns: percentile_f64(&search_times, 5.0),
            p95_search_ns: percentile_f64(&search_times, 95.0),
            median_linear_ns: med_linear,
            p5_linear_ns: percentile_f64(&linear_times, 5.0),
            p95_linear_ns: percentile_f64(&linear_times, 95.0),
            speedup_median: speedup,
            all_correct,
        });
    }

    results
}

// ---------------------------------------------------------------------------
// JSON output
// ---------------------------------------------------------------------------

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

fn write_json(
    path: &str,
    env_data: &[(String, String)],
    results: &[MutationResult],
    config: &[(String, String)],
) {
    let mut out = String::new();
    out.push_str("{\n");
    out.push_str("  \"benchmark\": \"mutation_localization_rust\",\n");

    out.push_str("  \"environment\": {\n");
    for (i, (k, v)) in env_data.iter().enumerate() {
        let comma = if i + 1 < env_data.len() { "," } else { "" };
        out.push_str(&format!("    \"{}\": \"{}\"{}\n", escape_json(k), escape_json(v), comma));
    }
    out.push_str("  },\n");

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

    out.push_str("  \"results\": [\n");
    for (i, r) in results.iter().enumerate() {
        let comma = if i + 1 < results.len() { "," } else { "" };
        out.push_str(&format!(
            "    {{\n\
             \x20     \"region_size\": {},\n\
             \x20     \"trials\": {},\n\
             \x20     \"expected_comparisons\": {},\n\
             \x20     \"median_comparisons\": {},\n\
             \x20     \"median_search_ns\": {:.1},\n\
             \x20     \"p5_search_ns\": {:.1},\n\
             \x20     \"p95_search_ns\": {:.1},\n\
             \x20     \"median_linear_ns\": {:.1},\n\
             \x20     \"p5_linear_ns\": {:.1},\n\
             \x20     \"p95_linear_ns\": {:.1},\n\
             \x20     \"speedup_median\": {:.2},\n\
             \x20     \"all_correct\": {}\n\
             \x20   }}{}\n",
            r.region_size, r.trials, r.expected_comparisons, r.median_comparisons,
            r.median_search_ns, r.p5_search_ns, r.p95_search_ns,
            r.median_linear_ns, r.p5_linear_ns, r.p95_linear_ns,
            r.speedup_median, r.all_correct, comma
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

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--fasta" => { i += 1; fasta_path = args[i].clone(); }
            "--output" => { i += 1; output_dir = Some(args[i].clone()); }
            other => panic!("Unknown argument: {}", other),
        }
        i += 1;
    }

    if fasta_path.is_empty() {
        eprintln!("Usage: bench_mutation_g3 --fasta <path> [--output <dir>]");
        std::process::exit(1);
    }

    eprintln!("hashrope-bio E-G3 (Rust): Mutation Localization via Binary Search");
    eprintln!("Build profile: {}", if cfg!(debug_assertions) { "debug" } else { "release" });
    eprintln!("FASTA: {}", fasta_path);
    eprintln!("Chunk size: {}", CHUNK_SIZE);
    eprintln!("Trials per size: {}", TRIALS_PER_SIZE);
    eprintln!("RNG seed: {}", RNG_SEED);

    // Load sequence and report N-content
    eprintln!("\nLoading sequence ...");
    let t0 = Instant::now();
    let (seq, seq_name) = load_fasta(&fasta_path);
    let n_count = seq.iter().filter(|&&b| !is_valid_nucleotide(b)).count();
    let valid_count = seq.len() - n_count;
    eprintln!("  Loaded: {}, {} bp ({:.1} MB) in {:.3}s",
              seq_name, seq.len(), seq.len() as f64 / (1024.0*1024.0), t0.elapsed().as_secs_f64());
    eprintln!("  Valid nucleotides (ACGT): {} ({:.1}%)", valid_count, 100.0 * valid_count as f64 / seq.len() as f64);
    eprintln!("  Non-ACGT (N etc.): {} ({:.1}%)", n_count, 100.0 * n_count as f64 / seq.len() as f64);

    // Run benchmark
    eprintln!("\n--- Mutation Localization Benchmark (chunk_size={}) ---", CHUNK_SIZE);
    let results = bench_mutation(&seq);

    // Summary table
    eprintln!("\n{}", "=".repeat(110));
    eprintln!("E-G3: Mutation Localization Summary");
    eprintln!("{}", "=".repeat(110));
    eprintln!("{:>12} {:>8} {:>8} {:>12} {:>12} {:>10} {:>8}",
              "N", "⌈log₂N⌉", "comps", "search_ns", "linear_ns", "speedup", "correct");
    eprintln!("{}", "-".repeat(110));
    for r in &results {
        eprintln!("{:>12} {:>8} {:>8} {:>12.0} {:>12.0} {:>9.0}× {:>8}",
                  r.region_size, r.expected_comparisons, r.median_comparisons,
                  r.median_search_ns, r.median_linear_ns,
                  r.speedup_median, r.all_correct);
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
            ("valid_nucleotides".into(), valid_count.to_string()),
            ("non_acgt_count".into(), n_count.to_string()),
            ("chunk_size".into(), CHUNK_SIZE.to_string()),
            ("region_sizes".into(), format!("[{}]", REGION_SIZES.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(", "))),
            ("trials_per_size".into(), TRIALS_PER_SIZE.to_string()),
            ("rng_seed".into(), RNG_SEED.to_string()),
            ("hash_base".into(), "131".into()),
            ("hash_prime".into(), "2305843009213693951".into()),
        ];

        let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs();
        let days = now / 86400;
        let tod = now % 86400;
        let (y, m, d) = days_to_ymd(days);
        let ts = format!("{:04}{:02}{:02}T{:02}{:02}{:02}Z", y, m, d, tod/3600, (tod%3600)/60, tod%60);

        let latest = format!("{}/mutation_localization_rust.json", dir);
        let archive = format!("{}/mutation_localization_rust_{}.json", dir, ts);

        write_json(&latest, &env_data, &results, &config);
        write_json(&archive, &env_data, &results, &config);

        eprintln!("\n  Results saved:");
        eprintln!("    Latest:  {}", latest);
        eprintln!("    Archive: {}", archive);
    }
}
