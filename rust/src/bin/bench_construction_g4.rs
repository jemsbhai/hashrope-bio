//! E-G4: Rope Construction Cost and Amortization (Rust)
//!
//! Mirrors the Python bench_construction.py for direct cross-language comparison.
//! Measures construction time, tree metrics, and amortization at genome scale.
//!
//! Usage:
//!     cargo run --release --bin bench_construction_g4 -- --fasta ../data/chr22.fa --output ../results/

use std::env;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::process::Command;
use std::time::Instant;

use hashrope::rope::{Arena, Node, NodeInner};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

const DEFAULT_CHUNK_SIZES: &[usize] = &[256, 1024, 4096, 16384];
const AMORT_QUERY_LEN: u64 = 10_000;
const AMORT_QUERY_ITERS: u64 = 100_000;
const AMORT_WARMUP: u64 = 10_000;

// ---------------------------------------------------------------------------
// FASTA loading (pure Rust, zero deps)
// ---------------------------------------------------------------------------

fn load_fasta(path: &str) -> (Vec<u8>, String) {
    let file = File::open(path).unwrap_or_else(|e| panic!("Cannot open {}: {}", path, e));
    let reader = BufReader::new(file);

    let mut seq_name = String::new();
    let mut seq = Vec::new();
    let mut found_header = false;

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        let line = line.trim_end();
        if line.starts_with('>') {
            if found_header {
                break; // second sequence — stop
            }
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

    if seq.is_empty() {
        panic!("No sequence data found in {}", path);
    }

    (seq, seq_name)
}

// ---------------------------------------------------------------------------
// Tree introspection (walks only reachable nodes)
// ---------------------------------------------------------------------------

struct NodeCounts {
    leaves: usize,
    internals: usize,
    repeats: usize,
}

impl NodeCounts {
    fn total(&self) -> usize {
        self.leaves + self.internals + self.repeats
    }
}

fn count_reachable(arena: &Arena, node: Node) -> NodeCounts {
    let mut counts = NodeCounts { leaves: 0, internals: 0, repeats: 0 };
    if let Some(id) = node {
        count_walk(arena, id, &mut counts);
    }
    counts
}

fn count_walk(arena: &Arena, id: u32, counts: &mut NodeCounts) {
    match arena.node(id) {
        NodeInner::Leaf { .. } => counts.leaves += 1,
        NodeInner::Internal { left, right, .. } => {
            counts.internals += 1;
            count_walk(arena, *left, counts);
            count_walk(arena, *right, counts);
        }
        NodeInner::Repeat { child, .. } => {
            counts.repeats += 1;
            count_walk(arena, *child, counts);
        }
    }
}

// ---------------------------------------------------------------------------
// Core benchmark
// ---------------------------------------------------------------------------

struct ConstructionResult {
    chunk_size: usize,
    seq_len: usize,
    leaves: usize,
    internals: usize,
    total_nodes: usize,
    arena_nodes: usize, // includes unreachable from rebalancing
    height: u64,
    construction_time_s: f64,
    throughput_mbs: f64,
    hash_verified: bool,
}

struct AmortizationResult {
    chunk_size: usize,
    construction_time_s: f64,
    query_len: u64,
    hashrope_query_ns: f64,
    baseline_hash_ns: f64,
    saving_per_query_ns: f64,
    amortization_queries: i64,
}

fn build_rope_chunked(arena: &mut Arena, seq: &[u8], chunk_size: usize) -> Node {
    let mut rope: Node = None;
    for chunk in seq.chunks(chunk_size) {
        let leaf = arena.from_bytes(chunk);
        rope = arena.concat(rope, leaf);
    }
    rope
}

fn bench_construction(seq: &[u8], chunk_sizes: &[usize]) -> Vec<ConstructionResult> {
    let mut results = Vec::new();

    // Compute ground-truth hash once
    let mut verify_arena = Arena::new();
    let full_hash = verify_arena.hash_bytes(seq);

    for &cs in chunk_sizes {
        eprint!("  chunk_size={} ... ", cs);

        let mut arena = Arena::new();
        let t0 = Instant::now();
        let rope = build_rope_chunked(&mut arena, seq, cs);
        let elapsed = t0.elapsed().as_secs_f64();

        let height = arena.height(rope);
        let counts = count_reachable(&arena, rope);
        let arena_total = arena.node_count();
        let rope_hash = arena.hash(rope);
        let hash_ok = rope_hash == full_hash;
        let throughput = (seq.len() as f64 / (1024.0 * 1024.0)) / elapsed;

        if hash_ok {
            eprintln!("OK ({:.4}s, height={}, {} leaves, arena={})",
                      elapsed, height, counts.leaves, arena_total);
        } else {
            eprintln!("HASH MISMATCH! expected={} got={}", full_hash, rope_hash);
        }

        results.push(ConstructionResult {
            chunk_size: cs,
            seq_len: seq.len(),
            leaves: counts.leaves,
            internals: counts.internals,
            total_nodes: counts.total(),
            arena_nodes: arena_total,
            height,
            construction_time_s: elapsed,
            throughput_mbs: throughput,
            hash_verified: hash_ok,
        });
    }

    results
}

fn bench_amortization(
    seq: &[u8],
    chunk_sizes: &[usize],
    construction_results: &[ConstructionResult],
) -> Vec<AmortizationResult> {
    let mut results = Vec::new();
    let query_len = std::cmp::min(AMORT_QUERY_LEN, seq.len() as u64 / 2);
    let start = seq.len() as u64 / 4;

    for (i, &cs) in chunk_sizes.iter().enumerate() {
        eprint!("  Amortization chunk_size={} (query_len={}) ... ", cs, query_len);

        let mut arena = Arena::new();
        let rope = build_rope_chunked(&mut arena, seq, cs);

        // --- hashrope substr_hash timing ---
        for _ in 0..AMORT_WARMUP {
            std::hint::black_box(arena.substr_hash(rope, start, query_len));
        }
        let t0 = Instant::now();
        for _ in 0..AMORT_QUERY_ITERS {
            std::hint::black_box(arena.substr_hash(rope, start, query_len));
        }
        let hashrope_ns = t0.elapsed().as_nanos() as f64 / AMORT_QUERY_ITERS as f64;

        // --- baseline: slice + hash ---
        let slice = &seq[start as usize..(start + query_len) as usize];
        for _ in 0..AMORT_WARMUP {
            std::hint::black_box(arena.hash_bytes(slice));
        }
        let t0 = Instant::now();
        for _ in 0..AMORT_QUERY_ITERS {
            std::hint::black_box(arena.hash_bytes(slice));
        }
        let baseline_ns = t0.elapsed().as_nanos() as f64 / AMORT_QUERY_ITERS as f64;

        let saving_ns = baseline_ns - hashrope_ns;
        let amort_queries = if saving_ns > 0.0 {
            (construction_results[i].construction_time_s * 1e9 / saving_ns) as i64 + 1
        } else {
            -1 // hashrope is slower
        };

        eprintln!("hashrope={:.0}ns, baseline={:.0}ns, saving={:.0}ns/query, amort={}",
                  hashrope_ns, baseline_ns, saving_ns, amort_queries);

        results.push(AmortizationResult {
            chunk_size: cs,
            construction_time_s: construction_results[i].construction_time_s,
            query_len,
            hashrope_query_ns: hashrope_ns,
            baseline_hash_ns: baseline_ns,
            saving_per_query_ns: saving_ns,
            amortization_queries: amort_queries,
        });
    }

    results
}

// ---------------------------------------------------------------------------
// Environment metadata
// ---------------------------------------------------------------------------

fn get_environment() -> Vec<(String, String)> {
    let mut env_data = Vec::new();

    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_secs();
    env_data.push(("timestamp_unix".into(), now.to_string()));
    env_data.push(("os".into(), env::consts::OS.into()));
    env_data.push(("arch".into(), env::consts::ARCH.into()));

    // CPU name (Windows)
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

    // RAM (Windows)
    if let Ok(output) = Command::new("wmic").args(["computersystem", "get", "totalphysicalmemory"]).output() {
        let s = String::from_utf8_lossy(&output.stdout);
        for line in s.lines() {
            let line = line.trim();
            if !line.is_empty() && line != "TotalPhysicalMemory" {
                if let Ok(bytes) = line.parse::<u64>() {
                    env_data.push(("ram_gb".into(), format!("{:.1}", bytes as f64 / (1024.0 * 1024.0 * 1024.0))));
                }
                break;
            }
        }
    }

    // Rust version
    if let Ok(output) = Command::new("rustc").arg("--version").output() {
        let s = String::from_utf8_lossy(&output.stdout).trim().to_string();
        env_data.push(("rustc_version".into(), s));
    }

    // hashrope crate version
    env_data.push(("hashrope_crate_version".into(), "0.2.1".into()));
    env_data.push(("hashrope_bio_crate_version".into(), "0.1.0".into()));
    env_data.push(("language".into(), "rust".into()));
    env_data.push(("build_profile".into(), if cfg!(debug_assertions) { "debug" } else { "release" }.into()));

    env_data
}

// ---------------------------------------------------------------------------
// JSON output (manual — no serde dependency)
// ---------------------------------------------------------------------------

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"").replace('\n', "\\n")
}

fn write_json(
    path: &str,
    env_data: &[(String, String)],
    construction: &[ConstructionResult],
    amortization: &[AmortizationResult],
    config: &[(String, String)],
) {
    let mut out = String::new();
    out.push_str("{\n");
    out.push_str("  \"benchmark\": \"construction_rust\",\n");

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
        // Try to emit numbers without quotes
        if v.parse::<f64>().is_ok() || v == "true" || v == "false" || v.starts_with('[') {
            out.push_str(&format!("    \"{}\": {}{}\n", escape_json(k), v, comma));
        } else {
            out.push_str(&format!("    \"{}\": \"{}\"{}\n", escape_json(k), escape_json(v), comma));
        }
    }
    out.push_str("  },\n");

    // Construction results
    out.push_str("  \"construction_results\": [\n");
    for (i, r) in construction.iter().enumerate() {
        let comma = if i + 1 < construction.len() { "," } else { "" };
        out.push_str(&format!(
            "    {{\n\
             \x20     \"chunk_size\": {},\n\
             \x20     \"seq_len\": {},\n\
             \x20     \"leaves\": {},\n\
             \x20     \"internals\": {},\n\
             \x20     \"total_reachable_nodes\": {},\n\
             \x20     \"arena_total_nodes\": {},\n\
             \x20     \"height\": {},\n\
             \x20     \"construction_time_s\": {:.6},\n\
             \x20     \"throughput_mbs\": {:.2},\n\
             \x20     \"hash_verified\": {}\n\
             \x20   }}{}\n",
            r.chunk_size, r.seq_len, r.leaves, r.internals,
            r.total_nodes, r.arena_nodes, r.height,
            r.construction_time_s, r.throughput_mbs, r.hash_verified, comma
        ));
    }
    out.push_str("  ],\n");

    // Amortization results
    out.push_str("  \"amortization_results\": [\n");
    for (i, r) in amortization.iter().enumerate() {
        let comma = if i + 1 < amortization.len() { "," } else { "" };
        out.push_str(&format!(
            "    {{\n\
             \x20     \"chunk_size\": {},\n\
             \x20     \"construction_time_s\": {:.6},\n\
             \x20     \"query_len\": {},\n\
             \x20     \"hashrope_query_ns\": {:.2},\n\
             \x20     \"baseline_hash_ns\": {:.2},\n\
             \x20     \"saving_per_query_ns\": {:.2},\n\
             \x20     \"amortization_queries\": {}\n\
             \x20   }}{}\n",
            r.chunk_size, r.construction_time_s, r.query_len,
            r.hashrope_query_ns, r.baseline_hash_ns,
            r.saving_per_query_ns, r.amortization_queries, comma
        ));
    }
    out.push_str("  ]\n");
    out.push_str("}\n");

    fs::write(path, &out).unwrap_or_else(|e| panic!("Cannot write {}: {}", path, e));
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = env::args().collect();

    let mut fasta_path = String::new();
    let mut output_dir: Option<String> = None;
    let mut chunk_sizes_str: Option<String> = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--fasta" => { i += 1; fasta_path = args[i].clone(); }
            "--output" => { i += 1; output_dir = Some(args[i].clone()); }
            "--chunk-sizes" => { i += 1; chunk_sizes_str = Some(args[i].clone()); }
            other => panic!("Unknown argument: {}", other),
        }
        i += 1;
    }

    if fasta_path.is_empty() {
        eprintln!("Usage: bench_construction_g4 --fasta <path> [--output <dir>] [--chunk-sizes 256,1024,4096,16384]");
        std::process::exit(1);
    }

    let chunk_sizes: Vec<usize> = chunk_sizes_str
        .as_deref()
        .unwrap_or("256,1024,4096,16384")
        .split(',')
        .map(|s| s.trim().parse().expect("Invalid chunk size"))
        .collect();

    eprintln!("hashrope-bio E-G4 (Rust): Rope Construction Cost and Amortization");
    eprintln!("Build profile: {}", if cfg!(debug_assertions) { "debug" } else { "release" });
    eprintln!("FASTA: {}", fasta_path);
    eprintln!("Chunk sizes: {:?}", chunk_sizes);

    // Load sequence
    eprintln!("\nLoading sequence from {} ...", fasta_path);
    let t0 = Instant::now();
    let (seq, seq_name) = load_fasta(&fasta_path);
    let load_time = t0.elapsed().as_secs_f64();
    eprintln!("  Loaded: {}, {} bp ({:.1} MB) in {:.3}s",
              seq_name, seq.len(), seq.len() as f64 / (1024.0 * 1024.0), load_time);

    // Construction benchmark
    eprintln!("\n--- Construction Benchmark ---");
    let construction_results = bench_construction(&seq, &chunk_sizes);

    // Print table
    eprintln!("\n{:>8} {:>12} {:>8} {:>8} {:>8} {:>8} {:>10} {:>10}",
              "chunk", "seq_len", "leaves", "internal", "total", "height", "time_s", "MB/s");
    eprintln!("{}", "-".repeat(86));
    for r in &construction_results {
        eprintln!("{:>8} {:>12} {:>8} {:>8} {:>8} {:>8} {:>10.4} {:>10.1}",
                  r.chunk_size, r.seq_len, r.leaves, r.internals,
                  r.total_nodes, r.height, r.construction_time_s, r.throughput_mbs);
    }

    // Amortization benchmark
    eprintln!("\n--- Amortization Benchmark ---");
    let amortization_results = bench_amortization(&seq, &chunk_sizes, &construction_results);

    eprintln!("\n{:>8} {:>12} {:>10} {:>10} {:>10} {:>12}",
              "chunk", "construct_s", "rope_ns", "base_ns", "save_ns", "amort_q");
    eprintln!("{}", "-".repeat(66));
    for r in &amortization_results {
        let amort_str = if r.amortization_queries > 0 {
            format!("{}", r.amortization_queries)
        } else {
            "N/A".to_string()
        };
        eprintln!("{:>8} {:>12.4} {:>10.1} {:>10.1} {:>10.1} {:>12}",
                  r.chunk_size, r.construction_time_s,
                  r.hashrope_query_ns, r.baseline_hash_ns,
                  r.saving_per_query_ns, amort_str);
    }

    // Hash cross-validation with Python
    {
        let mut arena = Arena::new();
        let full_hash = arena.hash_bytes(&seq);
        eprintln!("\n--- Hash Cross-Validation ---");
        eprintln!("  Full sequence hash: {}", full_hash);
    }

    // Save JSON results
    if let Some(ref dir) = output_dir {
        fs::create_dir_all(dir).ok();

        let env_data = get_environment();

        let config: Vec<(String, String)> = vec![
            ("fasta".into(), fasta_path.clone()),
            ("seq_name".into(), seq_name.clone()),
            ("seq_len".into(), seq.len().to_string()),
            ("chunk_sizes".into(), format!("[{}]", chunk_sizes.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(", "))),
            ("amort_query_len".into(), AMORT_QUERY_LEN.to_string()),
            ("amort_query_iters".into(), AMORT_QUERY_ITERS.to_string()),
            ("amort_warmup".into(), AMORT_WARMUP.to_string()),
            ("hash_base".into(), "131".into()),
            ("hash_prime".into(), "2305843009213693951".into()),
        ];

        // Timestamp for archival
        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_secs();
        // Simple UTC timestamp: YYYYMMDDTHHMMSSZ
        let secs_in_day = 86400u64;
        let days = now / secs_in_day;
        let time_of_day = now % secs_in_day;
        let hours = time_of_day / 3600;
        let minutes = (time_of_day % 3600) / 60;
        let seconds = time_of_day % 60;
        // Approximate date calculation (good enough for filenames)
        let (year, month, day) = days_to_ymd(days);
        let timestamp = format!("{:04}{:02}{:02}T{:02}{:02}{:02}Z", year, month, day, hours, minutes, seconds);

        let latest = format!("{}/construction_rust.json", dir);
        let archive = format!("{}/construction_rust_{}.json", dir, timestamp);

        write_json(&latest, &env_data, &construction_results, &amortization_results, &config);
        write_json(&archive, &env_data, &construction_results, &amortization_results, &config);

        eprintln!("\n  Results saved:");
        eprintln!("    Latest:  {}", latest);
        eprintln!("    Archive: {}", archive);
    }
}

// Simple days-since-epoch to (year, month, day)
fn days_to_ymd(mut days: u64) -> (u64, u64, u64) {
    // Approximate — sufficient for timestamp filenames
    let mut year = 1970;
    loop {
        let days_in_year = if year % 4 == 0 && (year % 100 != 0 || year % 400 == 0) { 366 } else { 365 };
        if days < days_in_year {
            break;
        }
        days -= days_in_year;
        year += 1;
    }
    let month_days: &[u64] = if year % 4 == 0 && (year % 100 != 0 || year % 400 == 0) {
        &[31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    } else {
        &[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    };
    let mut month = 1;
    for &md in month_days {
        if days < md {
            break;
        }
        days -= md;
        month += 1;
    }
    (year, month, days + 1)
}
