//! E-G3 Real-Data Validation: ClinVar Pathogenic SNV Localization
//!
//! Uses real pathogenic/likely_pathogenic SNVs from ClinVar (NCBI) on chr22
//! to validate binary search mutation localization with gold-standard clinical data.
//!
//! Prerequisites:
//!     python scripts/extract_clinvar_chr22.py --fasta ../data/chr22.fa
//!
//! Usage:
//!     cargo run --release --bin bench_mutation_g3_clinvar -- \
//!         --fasta ../data/chr22.fa \
//!         --variants ../data/clinvar_chr22_pathogenic_snvs.tsv \
//!         --output ../../../results/

use std::env;
use std::fs;
use std::io::{BufRead, BufReader};
use std::process::Command;
use std::time::Instant;

use hashrope::rope::{Arena, Node};

const CHUNK_SIZE: usize = 256;

// ---------------------------------------------------------------------------
// FASTA / TSV loading
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
                if b != b' ' && b != b'\t' { seq.push(b.to_ascii_uppercase()); }
            }
        }
    }
    if seq.is_empty() { panic!("No sequence data in {}", path); }
    (seq, seq_name)
}

struct ClinVarSNV {
    pos_0based: usize,
    ref_allele: u8,
    alt_allele: u8,
    clinvar_id: String,
    clnsig: String,
    gene: String,
    disease: String,
}

fn load_clinvar_tsv(path: &str) -> Vec<ClinVarSNV> {
    let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("Cannot open {}: {}", path, e));
    let reader = BufReader::new(file);
    let mut variants = Vec::new();

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('#') || line.starts_with("pos_1based") { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 { continue; }

        let ref_validated: bool = fields[8].trim().eq_ignore_ascii_case("true");
        if !ref_validated { continue; } // Skip variants with REF mismatch

        variants.push(ClinVarSNV {
            pos_0based: fields[1].parse().expect("Invalid pos_0based"),
            ref_allele: fields[2].as_bytes()[0],
            alt_allele: fields[3].as_bytes()[0],
            clinvar_id: fields[4].to_string(),
            clnsig: fields[5].to_string(),
            gene: fields[6].to_string(),
            disease: fields[7].to_string(),
        });
    }

    variants
}

// ---------------------------------------------------------------------------
// Rope building + binary search (same as E-G3)
// ---------------------------------------------------------------------------

fn build_rope_chunked(arena: &mut Arena, seq: &[u8], chunk_size: usize) -> Node {
    let mut rope: Node = None;
    for chunk in seq.chunks(chunk_size) {
        let leaf = arena.from_bytes(chunk);
        rope = arena.concat(rope, leaf);
    }
    rope
}

fn binary_search_mutation(
    arena: &mut Arena, ref_rope: Node, sample_rope: Node,
    start: u64, length: u64,
) -> (u64, u64) {
    let mut lo = start;
    let mut hi = start + length;
    let mut comparisons = 0u64;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        let h_ref = arena.substr_hash(ref_rope, lo, mid - lo);
        let h_sam = arena.substr_hash(sample_rope, lo, mid - lo);
        comparisons += 1;
        if h_ref != h_sam { hi = mid; } else { lo = mid; }
    }
    (lo, comparisons)
}

// ---------------------------------------------------------------------------
// Core benchmark
// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut fasta_path = String::new();
    let mut variants_path = String::new();
    let mut output_dir: Option<String> = None;
    let mut max_variants: usize = 200; // Default: 200 for fast run (~3 min)

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--fasta" => { i += 1; fasta_path = args[i].clone(); }
            "--variants" => { i += 1; variants_path = args[i].clone(); }
            "--output" => { i += 1; output_dir = Some(args[i].clone()); }
            "--max-variants" => { i += 1; max_variants = args[i].parse().expect("Invalid --max-variants"); }
            "--all" => { max_variants = usize::MAX; }
            other => panic!("Unknown argument: {}", other),
        }
        i += 1;
    }

    if fasta_path.is_empty() || variants_path.is_empty() {
        eprintln!("Usage: bench_mutation_g3_clinvar --fasta <chr22.fa> --variants <clinvar.tsv> [--output <dir>] [--max-variants N] [--all]");
        std::process::exit(1);
    }

    eprintln!("hashrope-bio E-G3 ClinVar Validation (Rust)");
    eprintln!("Build profile: {}", if cfg!(debug_assertions) { "debug" } else { "release" });

    // Load sequence
    eprintln!("\nLoading chr22 ...");
    let t0 = Instant::now();
    let (seq, seq_name) = load_fasta(&fasta_path);
    eprintln!("  {} bp in {:.3}s", seq.len(), t0.elapsed().as_secs_f64());

    // Build reference rope ONCE (amortized, realistic usage)
    eprintln!("Building reference rope (chunk_size={}) ...", CHUNK_SIZE);
    let t0 = Instant::now();
    let mut arena = Arena::new();
    let ref_rope = build_rope_chunked(&mut arena, &seq, CHUNK_SIZE);
    let construction_time = t0.elapsed().as_secs_f64();
    let height = arena.height(ref_rope);
    let leaves = arena.weight(ref_rope);
    eprintln!("  Built in {:.3}s, {} leaves, height {}", construction_time, leaves, height);

    // Load ClinVar variants
    eprintln!("\nLoading ClinVar variants from {} ...", variants_path);
    let mut variants = load_clinvar_tsv(&variants_path);
    let total_available = variants.len();
    if variants.len() > max_variants {
        variants.truncate(max_variants);
    }
    eprintln!("  Loaded {} pathogenic SNVs, using {} (--max-variants {}{})",
              total_available, variants.len(), max_variants,
              if max_variants == usize::MAX { ", --all" } else { "" });
    eprintln!("  Estimated runtime: ~{:.0}s ({:.1} min) — 2 rope constructions per variant",
              variants.len() as f64 * 1.0, variants.len() as f64 * 1.0 / 60.0);

    // Validate each REF allele against our loaded sequence
    let mut ref_checks_passed = 0;
    let mut ref_checks_failed = 0;
    for v in &variants {
        if v.pos_0based < seq.len() && seq[v.pos_0based] == v.ref_allele {
            ref_checks_passed += 1;
        } else {
            ref_checks_failed += 1;
            eprintln!("  REF CHECK FAIL: pos={} expected={} actual={}",
                      v.pos_0based, v.ref_allele as char,
                      if v.pos_0based < seq.len() { seq[v.pos_0based] as char } else { '?' });
        }
    }
    eprintln!("  REF allele validation: {}/{} passed, {} failed",
              ref_checks_passed, variants.len(), ref_checks_failed);
    if ref_checks_failed > 0 {
        eprintln!("  FATAL: REF allele mismatches — aborting");
        std::process::exit(1);
    }

    // Run binary search on every ClinVar variant
    eprintln!("\n--- ClinVar Mutation Localization ({} variants) ---", variants.len());
    let seq_len = seq.len() as u64;
    let expected_comps = {
        let mut n = seq_len;
        let mut c = 0u64;
        while n > 1 { n = (n + 1) / 2; c += 1; }
        c
    };

    let mut correct = 0usize;
    let mut incorrect = 0usize;
    let mut total_comps = 0u64;
    let mut search_times = Vec::with_capacity(variants.len());
    let mut incorrect_details: Vec<String> = Vec::new();

    let overall_t0 = Instant::now();

    for v in &variants {
        // Create mutated sequence and build sample rope
        let mut sample = seq.clone();
        sample[v.pos_0based] = v.alt_allele;

        // Build sample rope in a SEPARATE arena to avoid bloating the reference arena
        let mut sample_arena = Arena::new();
        let sample_rope = build_rope_chunked(&mut sample_arena, &sample, CHUNK_SIZE);

        // We need both ropes in the same arena for binary search.
        // Rebuild ref in the sample arena too.
        let ref_rope_local = build_rope_chunked(&mut sample_arena, &seq, CHUNK_SIZE);

        let t0 = Instant::now();
        let (found, comps) = binary_search_mutation(
            &mut sample_arena, ref_rope_local, sample_rope, 0, seq_len,
        );
        let elapsed_ns = t0.elapsed().as_nanos() as f64;

        search_times.push(elapsed_ns);
        total_comps += comps;

        if found as usize == v.pos_0based {
            correct += 1;
        } else {
            incorrect += 1;
            if incorrect_details.len() < 20 {
                incorrect_details.push(format!(
                    "  ClinVar {} ({} {}): expected pos={}, found={}, gene={}, disease={}",
                    v.clinvar_id, v.clnsig, v.ref_allele as char,
                    v.pos_0based, found, v.gene, v.disease
                ));
            }
        }
    }

    let overall_elapsed = overall_t0.elapsed().as_secs_f64();

    // Statistics
    search_times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = search_times.len();
    let median_ns = if n % 2 == 0 { (search_times[n/2-1] + search_times[n/2]) / 2.0 } else { search_times[n/2] };
    let p5_ns = search_times[(0.05 * (n - 1) as f64).round() as usize];
    let p95_ns = search_times[(0.95 * (n - 1) as f64).round() as usize];
    let avg_comps = total_comps as f64 / variants.len() as f64;

    eprintln!("\n=== E-G3 ClinVar Validation Results ===");
    eprintln!("  Variants tested: {}", variants.len());
    eprintln!("  Correctly localized: {} ({:.1}%)", correct, 100.0 * correct as f64 / variants.len() as f64);
    eprintln!("  Incorrect: {}", incorrect);
    eprintln!("  Average comparisons: {:.1} (expected ⌈log₂({})⌉ = {})", avg_comps, seq_len, expected_comps);
    eprintln!("  Search time: median={:.0}ns, p5={:.0}ns, p95={:.0}ns", median_ns, p5_ns, p95_ns);
    eprintln!("  Total wall-clock: {:.2}s for {} variants ({:.1} variants/sec)",
              overall_elapsed, variants.len(), variants.len() as f64 / overall_elapsed);

    if !incorrect_details.is_empty() {
        eprintln!("\n  First {} mismatches:", incorrect_details.len());
        for d in &incorrect_details { eprintln!("{}", d); }
    }

    // Save JSON
    if let Some(ref dir) = output_dir {
        fs::create_dir_all(dir).ok();

        let mut env_data = Vec::new();
        let now = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs();
        env_data.push(("timestamp_unix".to_string(), now.to_string()));
        env_data.push(("os".into(), env::consts::OS.into()));
        env_data.push(("arch".into(), env::consts::ARCH.into()));
        if let Ok(output) = Command::new("wmic").args(["cpu", "get", "name"]).output() {
            let s = String::from_utf8_lossy(&output.stdout);
            for line in s.lines() {
                let line = line.trim();
                if !line.is_empty() && line != "Name" { env_data.push(("cpu_name".into(), line.into())); break; }
            }
        }
        if let Ok(output) = Command::new("rustc").arg("--version").output() {
            env_data.push(("rustc_version".into(), String::from_utf8_lossy(&output.stdout).trim().to_string()));
        }
        env_data.push(("language".into(), "rust".into()));
        env_data.push(("build_profile".into(), if cfg!(debug_assertions) { "debug" } else { "release" }.into()));

        // Build JSON manually
        let mut out = String::new();
        out.push_str("{\n");
        out.push_str("  \"benchmark\": \"mutation_localization_clinvar_rust\",\n");

        out.push_str("  \"environment\": {\n");
        for (i, (k, v)) in env_data.iter().enumerate() {
            let comma = if i + 1 < env_data.len() { "," } else { "" };
            out.push_str(&format!("    \"{}\": \"{}\"{}\n", k, v.replace('\\', "\\\\").replace('"', "\\\""), comma));
        }
        out.push_str("  },\n");

        out.push_str("  \"config\": {\n");
        out.push_str(&format!("    \"fasta\": \"{}\",\n", fasta_path.replace('\\', "\\\\")));
        out.push_str(&format!("    \"variants_file\": \"{}\",\n", variants_path.replace('\\', "\\\\")));
        out.push_str(&format!("    \"seq_name\": \"{}\",\n", seq_name));
        out.push_str(&format!("    \"seq_len\": {},\n", seq.len()));
        out.push_str(&format!("    \"chunk_size\": {},\n", CHUNK_SIZE));
        out.push_str(&format!("    \"total_clinvar_snvs\": {},\n", total_available));
        out.push_str(&format!("    \"variants_tested\": {},\n", variants.len()));
        out.push_str(&format!("    \"max_variants_setting\": {},\n", max_variants));
        out.push_str("    \"data_source\": \"ClinVar VCF (GRCh38), NCBI FTP\",\n");
        out.push_str("    \"citation\": \"Landrum MJ et al. Nucleic Acids Res. 2018;46(D1):D1062-D1067. PMID: 29165669\"\n");
        out.push_str("  },\n");

        out.push_str("  \"results\": {\n");
        out.push_str(&format!("    \"variants_tested\": {},\n", variants.len()));
        out.push_str(&format!("    \"correctly_localized\": {},\n", correct));
        out.push_str(&format!("    \"incorrect\": {},\n", incorrect));
        out.push_str(&format!("    \"accuracy_pct\": {:.4},\n", 100.0 * correct as f64 / variants.len() as f64));
        out.push_str(&format!("    \"avg_comparisons\": {:.2},\n", avg_comps));
        out.push_str(&format!("    \"expected_comparisons\": {},\n", expected_comps));
        out.push_str(&format!("    \"median_search_ns\": {:.1},\n", median_ns));
        out.push_str(&format!("    \"p5_search_ns\": {:.1},\n", p5_ns));
        out.push_str(&format!("    \"p95_search_ns\": {:.1},\n", p95_ns));
        out.push_str(&format!("    \"total_wall_clock_s\": {:.3},\n", overall_elapsed));
        out.push_str(&format!("    \"variants_per_sec\": {:.1},\n", variants.len() as f64 / overall_elapsed));
        out.push_str(&format!("    \"ref_allele_validations_passed\": {},\n", ref_checks_passed));
        out.push_str(&format!("    \"ref_allele_validations_failed\": {}\n", ref_checks_failed));
        out.push_str("  }\n");
        out.push_str("}\n");

        let days = now / 86400;
        let tod = now % 86400;
        let (y, m, d) = days_to_ymd(days);
        let ts = format!("{:04}{:02}{:02}T{:02}{:02}{:02}Z", y, m, d, tod/3600, (tod%3600)/60, tod%60);

        let latest = format!("{}/clinvar_validation_rust.json", dir);
        let archive = format!("{}/clinvar_validation_rust_{}.json", dir, ts);
        fs::write(&latest, &out).unwrap();
        fs::write(&archive, &out).unwrap();
        eprintln!("\n  Results saved:");
        eprintln!("    Latest:  {}", latest);
        eprintln!("    Archive: {}", archive);
    }
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
