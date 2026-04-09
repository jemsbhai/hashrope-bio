//! E-CG1: Cross-Species Gene Conservation Screening (Rust)
//!
//! Screens 643 chr22 genes for human-chimp divergence using hashrope
//! on UCSC LASTZ pairwise alignment blocks.
//!
//! For each gene:
//!   1. Find overlapping alignment blocks
//!   2. Build hashrope for both human and chimp aligned sequences
//!   3. Compare full-sequence hashes (O(log w) each)
//!   4. For diverged genes: binary search to localize first divergence
//!   5. Time hashrope vs linear byte baseline
//!
//! Usage:
//!     cargo run --release --bin bench_cross_species_cg1 -- \
//!       --axt ../data/hg38_panTro6_chr22.axt \
//!       --genes ../data/chr22_genes.tsv \
//!       --output ../results/

use std::env;
use std::fs;
use std::io::{BufRead, BufReader};
use std::process::Command;
use std::time::Instant;

use hashrope::rope::{Arena, Node};

const CHUNK_SIZE: usize = 256;

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

struct AxtBlock {
    t_start: u64,
    t_end: u64,
    t_seq: Vec<u8>,
    q_seq: Vec<u8>,
}

struct Gene {
    name: String,
    tx_start: u64,
    tx_end: u64,
    tx_len: u64,
}

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

fn load_axt(path: &str) -> Vec<AxtBlock> {
    let file = fs::File::open(path).unwrap_or_else(|e| panic!("Cannot open {}: {}", path, e));
    let reader = BufReader::new(file);
    let mut blocks = Vec::new();
    let mut state = 0u8; // 0=header, 1=target, 2=query, 3=blank
    let mut t_start = 0u64;
    let mut t_end = 0u64;
    let mut t_seq = Vec::new();

    for line in reader.lines() {
        let line = line.expect("read error");
        let line = line.trim_end();
        match state {
            0 => {
                if line.is_empty() || line.starts_with('#') { continue; }
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 9 {
                    t_start = parts[2].parse().unwrap();
                    t_end = parts[3].parse().unwrap();
                    state = 1;
                }
            }
            1 => {
                t_seq = line.as_bytes().to_vec();
                for b in t_seq.iter_mut() { *b = b.to_ascii_uppercase(); }
                state = 2;
            }
            2 => {
                let mut q_seq: Vec<u8> = line.as_bytes().to_vec();
                for b in q_seq.iter_mut() { *b = b.to_ascii_uppercase(); }
                blocks.push(AxtBlock { t_start, t_end, t_seq: t_seq.clone(), q_seq });
                state = 3;
            }
            3 => { state = 0; }
            _ => {}
        }
    }
    blocks.sort_by_key(|b| b.t_start);
    blocks
}

fn load_genes(path: &str) -> Vec<Gene> {
    let file = fs::File::open(path).unwrap_or_else(|e| panic!("Cannot open {}: {}", path, e));
    let reader = BufReader::new(file);
    let mut genes = Vec::new();
    let mut first = true;
    for line in reader.lines() {
        let line = line.expect("read error");
        if first { first = false; continue; } // skip header
        let parts: Vec<&str> = line.trim().split('\t').collect();
        if parts.len() < 6 { continue; }
        genes.push(Gene {
            name: parts[0].to_string(),
            tx_start: parts[3].parse().unwrap(),
            tx_end: parts[4].parse().unwrap(),
            tx_len: parts[5].parse().unwrap(),
        });
    }
    genes
}

// ---------------------------------------------------------------------------
// Alignment region extraction
// ---------------------------------------------------------------------------

fn extract_aligned_region(block: &AxtBlock, start: u64, end: u64) -> (Vec<u8>, Vec<u8>) {
    let region_start = start.max(block.t_start);
    let region_end = end.min(block.t_end);
    if region_start >= region_end { return (Vec::new(), Vec::new()); }

    let mut genome_pos = block.t_start;
    let mut align_start: Option<usize> = None;
    let mut align_end: Option<usize> = None;

    for (i, &t_byte) in block.t_seq.iter().enumerate() {
        if t_byte != b'-' {
            if genome_pos == region_start && align_start.is_none() {
                align_start = Some(i);
            }
            genome_pos += 1;
            if genome_pos == region_end {
                align_end = Some(i + 1);
                break;
            }
        }
    }

    match (align_start, align_end) {
        (Some(s), Some(e)) => (block.t_seq[s..e].to_vec(), block.q_seq[s..e].to_vec()),
        (Some(s), None) => (block.t_seq[s..].to_vec(), block.q_seq[s..].to_vec()),
        _ => (Vec::new(), Vec::new()),
    }
}

// ---------------------------------------------------------------------------
// Rope building
// ---------------------------------------------------------------------------

fn build_rope(arena: &mut Arena, seq: &[u8]) -> Node {
    let mut rope: Node = None;
    for chunk in seq.chunks(CHUNK_SIZE) {
        let leaf = arena.from_bytes(chunk);
        rope = arena.concat(rope, leaf);
    }
    rope
}

// ---------------------------------------------------------------------------
// Environment metadata
// ---------------------------------------------------------------------------

fn get_os_info() -> String {
    #[cfg(target_os = "windows")]
    { "Windows".to_string() }
    #[cfg(target_os = "linux")]
    { "Linux".to_string() }
    #[cfg(target_os = "macos")]
    { "macOS".to_string() }
    #[cfg(not(any(target_os = "windows", target_os = "linux", target_os = "macos")))]
    { "Unknown".to_string() }
}

fn get_rustc_version() -> String {
    Command::new("rustc").arg("--version").output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|_| "unknown".into())
}

fn days_to_ymd(days: u64) -> (u64, u64, u64) {
    // Simplified — good enough for timestamps
    let y = 1970 + days / 365;
    let rem = days % 365;
    let m = rem / 30 + 1;
    let d = rem % 30 + 1;
    (y, m.min(12), d.min(28))
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut axt_path = String::new();
    let mut genes_path = String::new();
    let mut output_dir: Option<String> = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--axt" => { i += 1; axt_path = args[i].clone(); }
            "--genes" => { i += 1; genes_path = args[i].clone(); }
            "--output" => { i += 1; output_dir = Some(args[i].clone()); }
            other => panic!("Unknown argument: {}", other),
        }
        i += 1;
    }

    if axt_path.is_empty() || genes_path.is_empty() {
        eprintln!("Usage: bench_cross_species_cg1 --axt <path> --genes <path> [--output <dir>]");
        std::process::exit(1);
    }

    eprintln!("hashrope-bio E-CG1 (Rust): Cross-Species Gene Conservation");
    eprintln!("Build: {}", if cfg!(debug_assertions) { "debug" } else { "release" });
    eprintln!("Chunk size: {}", CHUNK_SIZE);

    // Load data
    eprintln!("\nLoading alignment blocks...");
    let t0 = Instant::now();
    let blocks = load_axt(&axt_path);
    eprintln!("  {} blocks in {:.3}s", blocks.len(), t0.elapsed().as_secs_f64());

    eprintln!("Loading gene annotations...");
    let genes = load_genes(&genes_path);
    eprintln!("  {} genes", genes.len());

    // Screen each gene
    eprintln!("\nScreening {} genes...", genes.len());

    let mut conserved_count = 0u64;
    let mut diverged_count = 0u64;
    let mut no_alignment_count = 0u64;
    let mut hash_errors = 0u64;
    let mut gene_results: Vec<String> = Vec::new(); // JSON fragments

    // Timing accumulators
    let mut total_hashrope_screen_ns = 0u64;
    let mut total_baseline_screen_ns = 0u64;
    let mut total_hashrope_search_ns = 0u64;
    let mut total_linear_count_ns = 0u64;
    let mut total_search_genes = 0u64;

    // Per-gene detailed results for top-N reporting
    struct GeneResult {
        name: String,
        aligned_len: usize,
        is_conserved: bool,
        hashrope_screen_ns: u64,
        baseline_screen_ns: u64,
        // Divergence details (only if diverged)
        mismatches: usize,
        mismatch_rate: f64,
        search_comparisons: u64,
        search_ns: u64,
        linear_count_ns: u64,
    }

    let mut detailed: Vec<GeneResult> = Vec::new();

    for (gi, gene) in genes.iter().enumerate() {
        // Find overlapping blocks
        let mut human_parts: Vec<Vec<u8>> = Vec::new();
        let mut chimp_parts: Vec<Vec<u8>> = Vec::new();

        for block in &blocks {
            if block.t_start >= gene.tx_end || block.t_end <= gene.tx_start { continue; }
            let (h, c) = extract_aligned_region(block, gene.tx_start, gene.tx_end);
            if !h.is_empty() {
                human_parts.push(h);
                chimp_parts.push(c);
            }
        }

        if human_parts.is_empty() {
            no_alignment_count += 1;
            continue;
        }

        let human_seq: Vec<u8> = human_parts.into_iter().flatten().collect();
        let chimp_seq: Vec<u8> = chimp_parts.into_iter().flatten().collect();
        let aligned_len = human_seq.len();
        assert_eq!(aligned_len, chimp_seq.len());

        // --- Hashrope screen ---
        let t0 = Instant::now();
        let mut arena = Arena::new();
        let h_rope = build_rope(&mut arena, &human_seq);
        let c_rope = build_rope(&mut arena, &chimp_seq);
        let h_hash = arena.substr_hash(h_rope, 0, aligned_len as u64);
        let c_hash = arena.substr_hash(c_rope, 0, aligned_len as u64);
        let is_conserved = h_hash == c_hash;
        let hashrope_screen_ns = t0.elapsed().as_nanos() as u64;

        // --- Byte baseline ---
        let t0 = Instant::now();
        let byte_match = human_seq == chimp_seq;
        let baseline_screen_ns = t0.elapsed().as_nanos() as u64;

        // Correctness check
        if is_conserved != byte_match {
            hash_errors += 1;
            eprintln!("  HASH ERROR: {} hashrope={} bytes={}", gene.name, is_conserved, byte_match);
            continue;
        }

        total_hashrope_screen_ns += hashrope_screen_ns;
        total_baseline_screen_ns += baseline_screen_ns;

        let mut mismatches = 0usize;
        let mut search_comparisons = 0u64;
        let mut search_ns = 0u64;
        let mut linear_count_ns = 0u64;

        if is_conserved {
            conserved_count += 1;
        } else {
            diverged_count += 1;

            // Binary search for first divergence
            let t0 = Instant::now();
            let mut lo = 0u64;
            let mut hi = aligned_len as u64;
            let mut comps = 0u64;
            while hi - lo > 1 {
                let mid = (lo + hi) / 2;
                let hh = arena.substr_hash(h_rope, lo, mid - lo);
                let ch = arena.substr_hash(c_rope, lo, mid - lo);
                comps += 1;
                if hh != ch { hi = mid; } else { lo = mid; }
            }
            search_ns = t0.elapsed().as_nanos() as u64;
            search_comparisons = comps;

            // Linear mismatch count (ground truth)
            let t0 = Instant::now();
            for j in 0..aligned_len {
                if human_seq[j] != chimp_seq[j] { mismatches += 1; }
            }
            linear_count_ns = t0.elapsed().as_nanos() as u64;

            total_hashrope_search_ns += search_ns;
            total_linear_count_ns += linear_count_ns;
            total_search_genes += 1;
        }

        detailed.push(GeneResult {
            name: gene.name.clone(),
            aligned_len,
            is_conserved,
            hashrope_screen_ns,
            baseline_screen_ns,
            mismatches,
            mismatch_rate: if aligned_len > 0 { mismatches as f64 / aligned_len as f64 * 100.0 } else { 0.0 },
            search_comparisons,
            search_ns,
            linear_count_ns,
        });

        if (gi + 1) % 100 == 0 {
            eprintln!("  {}/{} genes...", gi + 1, genes.len());
        }
    }

    // -----------------------------------------------------------------------
    // Report
    // -----------------------------------------------------------------------
    let total_with_alignment = conserved_count + diverged_count;
    let pct_conserved = if total_with_alignment > 0 {
        conserved_count as f64 / total_with_alignment as f64 * 100.0
    } else { 0.0 };

    eprintln!("\n{}", "=".repeat(90));
    eprintln!("E-CG1 Results: Human vs Chimp chr22 Gene Conservation");
    eprintln!("{}", "=".repeat(90));
    eprintln!("  Genes screened:       {}", genes.len());
    eprintln!("  With alignment:       {}", total_with_alignment);
    eprintln!("  No alignment:         {}", no_alignment_count);
    eprintln!("  Hash errors:          {}", hash_errors);
    eprintln!("  ---");
    eprintln!("  Conserved (identical): {}", conserved_count);
    eprintln!("  Diverged:              {}", diverged_count);
    eprintln!("  Conservation rate:     {:.1}%", pct_conserved);
    eprintln!("  ---");
    eprintln!("  Screening time (hashrope): {:.3} ms", total_hashrope_screen_ns as f64 / 1e6);
    eprintln!("  Screening time (baseline): {:.3} ms", total_baseline_screen_ns as f64 / 1e6);
    let screen_speedup = if total_hashrope_screen_ns > 0 {
        total_baseline_screen_ns as f64 / total_hashrope_screen_ns as f64
    } else { 0.0 };
    eprintln!("  Screen speedup: {:.2}×", screen_speedup);

    // Top diverged genes
    let mut diverged: Vec<&GeneResult> = detailed.iter().filter(|g| !g.is_conserved).collect();
    diverged.sort_by(|a, b| b.mismatches.cmp(&a.mismatches));
    eprintln!("\n  Top 10 most diverged genes:");
    eprintln!("  {:<15} {:>8} {:>11} {:>7} {:>6} {:>12} {:>12} {:>8}",
              "Gene", "Aligned", "Mismatches", "Rate", "Comps", "Search(µs)", "Linear(µs)", "Speedup");
    eprintln!("  {}", "-".repeat(85));
    for g in diverged.iter().take(10) {
        let speedup = if g.search_ns > 0 { g.linear_count_ns as f64 / g.search_ns as f64 } else { 0.0 };
        eprintln!("  {:<15} {:>8} {:>11} {:>6.2}% {:>6} {:>11.1} {:>11.1} {:>7.1}×",
                  g.name, g.aligned_len, g.mismatches, g.mismatch_rate,
                  g.search_comparisons,
                  g.search_ns as f64 / 1000.0,
                  g.linear_count_ns as f64 / 1000.0,
                  speedup);
    }

    // Conserved genes
    let mut conserved: Vec<&GeneResult> = detailed.iter().filter(|g| g.is_conserved).collect();
    conserved.sort_by(|a, b| b.aligned_len.cmp(&a.aligned_len));
    if !conserved.is_empty() {
        eprintln!("\n  Conserved genes (identical, top 10 by size):");
        for g in conserved.iter().take(10) {
            eprintln!("    {:<15} {} bp", g.name, g.aligned_len);
        }
    }

    // Aggregate binary search stats
    if total_search_genes > 0 {
        eprintln!("\n  Binary search aggregate ({} diverged genes):", total_search_genes);
        eprintln!("    Total search time:  {:.3} ms", total_hashrope_search_ns as f64 / 1e6);
        eprintln!("    Total linear time:  {:.3} ms", total_linear_count_ns as f64 / 1e6);
        let agg_speedup = if total_hashrope_search_ns > 0 {
            total_linear_count_ns as f64 / total_hashrope_search_ns as f64
        } else { 0.0 };
        eprintln!("    Aggregate speedup:  {:.1}×", agg_speedup);
    }

    // Screen speedup by gene size bin
    eprintln!("\n  Screening speedup by gene size:");
    let bins: &[(usize, usize, &str)] = &[
        (0, 1000, "<1K"),
        (1000, 5000, "1K-5K"),
        (5000, 10000, "5K-10K"),
        (10000, 50000, "10K-50K"),
        (50000, 100000, "50K-100K"),
        (100000, 500000, "100K-500K"),
        (500000, usize::MAX, ">500K"),
    ];
    eprintln!("  {:<12} {:>6} {:>14} {:>14} {:>10}",
              "Size bin", "Genes", "Hashrope(µs)", "Baseline(µs)", "Speedup");
    eprintln!("  {}", "-".repeat(60));
    for &(lo, hi, label) in bins {
        let in_bin: Vec<&GeneResult> = detailed.iter()
            .filter(|g| g.aligned_len >= lo && g.aligned_len < hi)
            .collect();
        if in_bin.is_empty() { continue; }
        let n = in_bin.len();
        let hr_median = {
            let mut vals: Vec<u64> = in_bin.iter().map(|g| g.hashrope_screen_ns).collect();
            vals.sort();
            vals[vals.len() / 2]
        };
        let bl_median = {
            let mut vals: Vec<u64> = in_bin.iter().map(|g| g.baseline_screen_ns).collect();
            vals.sort();
            vals[vals.len() / 2]
        };
        let speedup = if hr_median > 0 { bl_median as f64 / hr_median as f64 } else { 0.0 };
        eprintln!("  {:<12} {:>6} {:>13.1} {:>13.1} {:>9.1}×",
                  label, n, hr_median as f64 / 1000.0, bl_median as f64 / 1000.0, speedup);
    }

    eprintln!("{}", "=".repeat(90));

    // Save JSON
    if let Some(ref dir) = output_dir {
        fs::create_dir_all(dir).ok();
        let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs();
        let days = now / 86400;
        let tod = now % 86400;
        let (y, m, d) = days_to_ymd(days);
        let ts = format!("{:04}{:02}{:02}T{:02}{:02}{:02}Z", y, m, d, tod/3600, (tod%3600)/60, tod%60);

        // Build JSON manually (no serde dependency)
        let mut json = String::from("{\n");
        json.push_str("  \"benchmark\": \"cross_species_conservation_rust\",\n");
        json.push_str(&format!("  \"os\": \"{}\",\n", get_os_info()));
        json.push_str(&format!("  \"rustc\": \"{}\",\n", get_rustc_version()));
        json.push_str(&format!("  \"chunk_size\": {},\n", CHUNK_SIZE));
        json.push_str(&format!("  \"alignment_blocks\": {},\n", blocks.len()));
        json.push_str(&format!("  \"genes_screened\": {},\n", genes.len()));
        json.push_str(&format!("  \"conserved\": {},\n", conserved_count));
        json.push_str(&format!("  \"diverged\": {},\n", diverged_count));
        json.push_str(&format!("  \"no_alignment\": {},\n", no_alignment_count));
        json.push_str(&format!("  \"hash_errors\": {},\n", hash_errors));
        json.push_str(&format!("  \"conservation_rate_pct\": {:.2},\n", pct_conserved));
        json.push_str(&format!("  \"total_hashrope_screen_ns\": {},\n", total_hashrope_screen_ns));
        json.push_str(&format!("  \"total_baseline_screen_ns\": {},\n", total_baseline_screen_ns));
        json.push_str(&format!("  \"total_hashrope_search_ns\": {},\n", total_hashrope_search_ns));
        json.push_str(&format!("  \"total_linear_count_ns\": {},\n", total_linear_count_ns));

        // Per-gene results
        json.push_str("  \"per_gene\": [\n");
        for (i, g) in detailed.iter().enumerate() {
            json.push_str(&format!(
                "    {{\"gene\":\"{}\",\"aligned_len\":{},\"conserved\":{},\"hashrope_screen_ns\":{},\"baseline_screen_ns\":{},\"mismatches\":{},\"mismatch_rate\":{:.4},\"search_comps\":{},\"search_ns\":{},\"linear_count_ns\":{}}}",
                g.name, g.aligned_len, g.is_conserved,
                g.hashrope_screen_ns, g.baseline_screen_ns,
                g.mismatches, g.mismatch_rate,
                g.search_comparisons, g.search_ns, g.linear_count_ns
            ));
            if i < detailed.len() - 1 { json.push_str(",\n"); } else { json.push('\n'); }
        }
        json.push_str("  ]\n");
        json.push_str("}\n");

        let latest = format!("{}/cross_species_rust.json", dir);
        let archive = format!("{}/cross_species_rust_{}.json", dir, ts);
        fs::write(&latest, &json).expect("write latest");
        fs::write(&archive, &json).expect("write archive");
        eprintln!("\nSaved: {}", latest);
        eprintln!("Saved: {}", archive);
    }
}
