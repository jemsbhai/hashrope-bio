//! Gene-level change detection via hash comparison.
//!
//! Given two ropes (e.g., reference and patient genome) and a list of
//! gene annotations, reports which genes and exons differ — in O(N · log w)
//! total time, where N is the number of regions and w is the rope weight.
//!
//! # Two-pass architecture
//!
//! - **Pass 1 (screen)**: O(log w) hash comparison per region — identifies changes
//! - **Pass 2 (detail)**: only changed regions need expensive byte-level analysis
//!
//! # Example
//!
//! ```no_run
//! use hashrope::Arena;
//! use hashrope_bio::gene_diff::*;
//!
//! let genes = load_gene_regions("chr22_genes.tsv").unwrap();
//! let mut arena = Arena::new();
//! let ref_rope = arena.from_bytes(b"ACGT...");
//! let sample_rope = arena.from_bytes(b"ACGT...");
//! let report = diff_genes(&mut arena, ref_rope, &mut arena, sample_rope, &genes, None);
//! ```

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use hashrope::Arena;
use hashrope::rope::Node;

// ---------------------------------------------------------------------------
// Data types
// ---------------------------------------------------------------------------

/// A single exon with genomic coordinates (0-based, half-open).
#[derive(Debug, Clone)]
pub struct Exon {
    /// 1-based exon number within the gene.
    pub index: usize,
    /// 0-based genomic start position.
    pub start: u64,
    /// 0-based genomic end position (exclusive).
    pub end: u64,
}

impl Exon {
    pub fn length(&self) -> u64 {
        self.end - self.start
    }
}

/// A gene with its transcript boundaries and exon structure.
#[derive(Debug, Clone)]
pub struct GeneRegion {
    pub name: String,
    pub transcript: String,
    pub chrom: String,
    pub strand: char,
    /// 0-based genomic start of transcript.
    pub tx_start: u64,
    /// 0-based genomic end of transcript (exclusive).
    pub tx_end: u64,
    pub exons: Vec<Exon>,
}

impl GeneRegion {
    pub fn tx_length(&self) -> u64 {
        self.tx_end - self.tx_start
    }

    pub fn exon_count(&self) -> usize {
        self.exons.len()
    }
}

/// Result of comparing one exon between two ropes.
#[derive(Debug)]
pub struct ExonDiff {
    pub exon_index: usize,
    pub start: u64,
    pub length: u64,
    pub ref_hash: u64,
    pub sample_hash: u64,
}

impl ExonDiff {
    pub fn changed(&self) -> bool {
        self.ref_hash != self.sample_hash
    }
}

/// Result of comparing one gene (body + exons) between two ropes.
#[derive(Debug)]
pub struct GeneDiff {
    pub gene_name: String,
    pub tx_start: u64,
    pub tx_length: u64,
    pub body_ref_hash: u64,
    pub body_sample_hash: u64,
    pub exon_diffs: Vec<ExonDiff>,
}

impl GeneDiff {
    pub fn body_changed(&self) -> bool {
        self.body_ref_hash != self.body_sample_hash
    }

    pub fn changed_exons(&self) -> Vec<&ExonDiff> {
        self.exon_diffs.iter().filter(|e| e.changed()).collect()
    }

    pub fn changed_exon_count(&self) -> usize {
        self.exon_diffs.iter().filter(|e| e.changed()).count()
    }
}

/// Full report from diffing a set of genes between two ropes.
#[derive(Debug)]
pub struct GeneDiffReport {
    pub gene_diffs: Vec<GeneDiff>,
}

impl GeneDiffReport {
    pub fn total_genes(&self) -> usize {
        self.gene_diffs.len()
    }

    pub fn changed_genes(&self) -> Vec<&GeneDiff> {
        self.gene_diffs.iter().filter(|g| g.body_changed()).collect()
    }

    pub fn unchanged_genes(&self) -> Vec<&GeneDiff> {
        self.gene_diffs.iter().filter(|g| !g.body_changed()).collect()
    }

    pub fn total_exons(&self) -> usize {
        self.gene_diffs.iter().map(|g| g.exon_diffs.len()).sum()
    }

    pub fn changed_exon_count(&self) -> usize {
        self.gene_diffs.iter().map(|g| g.changed_exon_count()).sum()
    }

    /// Human-readable summary string.
    pub fn summary(&self) -> String {
        let mut lines = vec![format!(
            "Gene diff: {}/{} genes changed, {}/{} exons changed",
            self.changed_genes().len(),
            self.total_genes(),
            self.changed_exon_count(),
            self.total_exons(),
        )];
        for gd in self.changed_genes() {
            lines.push(format!(
                "  {} ({} bp): body changed, {}/{} exons differ",
                gd.gene_name, gd.tx_length,
                gd.changed_exon_count(), gd.exon_diffs.len(),
            ));
        }
        lines.join("\n")
    }
}

// ---------------------------------------------------------------------------
// Gene loading
// ---------------------------------------------------------------------------

/// Load gene annotations from a TSV file (chr22_genes.tsv format).
///
/// Expected columns: gene_name, transcript, strand, tx_start, tx_end,
/// tx_len, cds_start, cds_end, exon_count, exon_starts, exon_ends, exon_lengths
pub fn load_gene_regions<P: AsRef<Path>>(path: P) -> Result<Vec<GeneRegion>, String> {
    let file = File::open(path.as_ref())
        .map_err(|e| format!("Cannot open {}: {}", path.as_ref().display(), e))?;
    let reader = BufReader::new(file);
    let mut genes = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| format!("Read error line {}: {}", i, e))?;
        if i == 0 { continue; } // skip header
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 { continue; }

        let name = fields[0].to_string();
        let transcript = fields[1].to_string();
        let strand = fields[2].chars().next().unwrap_or('+');
        let tx_start: u64 = fields[3].parse().map_err(|e| format!("Line {}: {}", i, e))?;
        let tx_end: u64 = fields[4].parse().map_err(|e| format!("Line {}: {}", i, e))?;

        let exon_starts: Vec<u64> = fields[9].split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().unwrap_or(0))
            .collect();
        let exon_ends: Vec<u64> = fields[10].split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().unwrap_or(0))
            .collect();

        let exons: Vec<Exon> = exon_starts.iter().zip(exon_ends.iter())
            .enumerate()
            .map(|(idx, (&s, &e))| Exon { index: idx + 1, start: s, end: e })
            .collect();

        genes.push(GeneRegion {
            name, transcript, chrom: "chr22".into(), strand,
            tx_start, tx_end, exons,
        });
    }

    genes.sort_by_key(|g| g.tx_start);
    Ok(genes)
}

// ---------------------------------------------------------------------------
// Core: gene-level diff
// ---------------------------------------------------------------------------

/// Compare gene regions between two ropes using a single shared arena.
///
/// For each gene, compares the gene body hash. If `include_exons` is true,
/// also compares each exon individually.
///
/// Total cost: O((G + E) · log w) where G = genes, E = exons, w = rope weight.
pub fn diff_genes_single_arena(
    arena: &mut Arena,
    ref_rope: Node,
    sample_rope: Node,
    genes: &[GeneRegion],
    seq_len: Option<u64>,
    include_exons: bool,
) -> GeneDiffReport {
    let mut gene_diffs = Vec::with_capacity(genes.len());

    for gene in genes {
        if let Some(sl) = seq_len {
            if gene.tx_end > sl { continue; }
        }

        let body_ref = arena.substr_hash(ref_rope, gene.tx_start, gene.tx_length());
        let body_sample = arena.substr_hash(sample_rope, gene.tx_start, gene.tx_length());

        let exon_diffs = if include_exons {
            gene.exons.iter().filter_map(|exon| {
                if let Some(sl) = seq_len {
                    if exon.end > sl { return None; }
                }
                if exon.length() == 0 { return None; }
                let ex_ref = arena.substr_hash(ref_rope, exon.start, exon.length());
                let ex_sample = arena.substr_hash(sample_rope, exon.start, exon.length());
                Some(ExonDiff {
                    exon_index: exon.index,
                    start: exon.start,
                    length: exon.length(),
                    ref_hash: ex_ref,
                    sample_hash: ex_sample,
                })
            }).collect()
        } else {
            Vec::new()
        };

        gene_diffs.push(GeneDiff {
            gene_name: gene.name.clone(),
            tx_start: gene.tx_start,
            tx_length: gene.tx_length(),
            body_ref_hash: body_ref,
            body_sample_hash: body_sample,
            exon_diffs,
        });
    }

    GeneDiffReport { gene_diffs }
}

/// Compare gene regions between two ropes in separate arenas.
///
/// Use this when reference and sample are built in different arenas.
pub fn diff_genes(
    arena_ref: &mut Arena,
    ref_rope: Node,
    arena_sample: &mut Arena,
    sample_rope: Node,
    genes: &[GeneRegion],
    seq_len: Option<u64>,
) -> GeneDiffReport {
    let mut gene_diffs = Vec::with_capacity(genes.len());

    for gene in genes {
        if let Some(sl) = seq_len {
            if gene.tx_end > sl { continue; }
        }

        let body_ref = arena_ref.substr_hash(ref_rope, gene.tx_start, gene.tx_length());
        let body_sample = arena_sample.substr_hash(sample_rope, gene.tx_start, gene.tx_length());

        let exon_diffs: Vec<ExonDiff> = gene.exons.iter().filter_map(|exon| {
            if let Some(sl) = seq_len {
                if exon.end > sl { return None; }
            }
            if exon.length() == 0 { return None; }
            let ex_ref = arena_ref.substr_hash(ref_rope, exon.start, exon.length());
            let ex_sample = arena_sample.substr_hash(sample_rope, exon.start, exon.length());
            Some(ExonDiff {
                exon_index: exon.index,
                start: exon.start,
                length: exon.length(),
                ref_hash: ex_ref,
                sample_hash: ex_sample,
            })
        }).collect();

        gene_diffs.push(GeneDiff {
            gene_name: gene.name.clone(),
            tx_start: gene.tx_start,
            tx_length: gene.tx_length(),
            body_ref_hash: body_ref,
            body_sample_hash: body_sample,
            exon_diffs,
        });
    }

    GeneDiffReport { gene_diffs }
}

/// Quick screen: return names of genes whose bodies differ.
///
/// Gene-body-only comparison (no exon detail). Fastest possible
/// screening pass — O(G · log w) total.
pub fn screen_genes_single_arena(
    arena: &mut Arena,
    ref_rope: Node,
    sample_rope: Node,
    genes: &[GeneRegion],
    seq_len: Option<u64>,
) -> Vec<String> {
    genes.iter().filter_map(|gene| {
        if let Some(sl) = seq_len {
            if gene.tx_end > sl { return None; }
        }
        let h_ref = arena.substr_hash(ref_rope, gene.tx_start, gene.tx_length());
        let h_sample = arena.substr_hash(sample_rope, gene.tx_start, gene.tx_length());
        if h_ref != h_sample { Some(gene.name.clone()) } else { None }
    }).collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use hashrope::Arena;

    fn build_rope_chunked(arena: &mut Arena, seq: &[u8], chunk_size: usize) -> Node {
        let mut rope: Node = None;
        for chunk in seq.chunks(chunk_size) {
            let leaf = arena.from_bytes(chunk);
            rope = arena.concat(rope, leaf);
        }
        rope
    }

    fn test_genes() -> Vec<GeneRegion> {
        vec![
            GeneRegion {
                name: "GENE_A".into(), transcript: "NM_001".into(),
                chrom: "chr1".into(), strand: '+',
                tx_start: 0, tx_end: 300,
                exons: vec![Exon { index: 1, start: 0, end: 100 },
                            Exon { index: 2, start: 150, end: 300 }],
            },
            GeneRegion {
                name: "GENE_B".into(), transcript: "NM_002".into(),
                chrom: "chr1".into(), strand: '+',
                tx_start: 400, tx_end: 700,
                exons: vec![Exon { index: 1, start: 400, end: 550 },
                            Exon { index: 2, start: 600, end: 700 }],
            },
            GeneRegion {
                name: "GENE_C".into(), transcript: "NM_003".into(),
                chrom: "chr1".into(), strand: '-',
                tx_start: 750, tx_end: 1000,
                exons: vec![Exon { index: 1, start: 750, end: 900 },
                            Exon { index: 2, start: 920, end: 1000 }],
            },
        ]
    }

    #[test]
    fn test_identical_no_changes() {
        let seq = b"ACGT".repeat(250);
        let genes = test_genes();
        let mut arena = Arena::new();
        let ref_rope = build_rope_chunked(&mut arena, &seq, 64);
        let sample_rope = build_rope_chunked(&mut arena, &seq, 64);

        let report = diff_genes_single_arena(
            &mut arena, ref_rope, sample_rope, &genes, None, true,
        );

        assert_eq!(report.total_genes(), 3);
        assert_eq!(report.changed_genes().len(), 0);
        assert_eq!(report.changed_exon_count(), 0);
    }

    #[test]
    fn test_mutation_detected_in_correct_gene() {
        let ref_seq = b"ACGT".repeat(250);
        let mut sample_seq = ref_seq.clone();
        sample_seq[500] = b'T'; // inside GENE_B (400-700)
        let genes = test_genes();

        let mut arena = Arena::new();
        let ref_rope = build_rope_chunked(&mut arena, &ref_seq, 64);
        let sample_rope = build_rope_chunked(&mut arena, &sample_seq, 64);

        let report = diff_genes_single_arena(
            &mut arena, ref_rope, sample_rope, &genes, None, true,
        );

        let changed: Vec<&str> = report.changed_genes().iter()
            .map(|g| g.gene_name.as_str()).collect();
        assert!(changed.contains(&"GENE_B"));
        assert!(!changed.contains(&"GENE_A"));
        assert!(!changed.contains(&"GENE_C"));
    }

    #[test]
    fn test_mutation_detected_in_correct_exon() {
        let ref_seq = b"ACGT".repeat(250);
        let mut sample_seq = ref_seq.clone();
        sample_seq[500] = b'T'; // inside GENE_B exon 1 (400-550)
        let genes = test_genes();

        let mut arena = Arena::new();
        let ref_rope = build_rope_chunked(&mut arena, &ref_seq, 64);
        let sample_rope = build_rope_chunked(&mut arena, &sample_seq, 64);

        let report = diff_genes_single_arena(
            &mut arena, ref_rope, sample_rope, &genes, None, true,
        );

        let gene_b = report.gene_diffs.iter()
            .find(|g| g.gene_name == "GENE_B").unwrap();
        let changed_idx: Vec<usize> = gene_b.changed_exons().iter()
            .map(|e| e.exon_index).collect();
        assert!(changed_idx.contains(&1));   // exon 1 (400-550) has the mutation
        assert!(!changed_idx.contains(&2));  // exon 2 (600-700) is clean
    }

    #[test]
    fn test_screen_genes() {
        let ref_seq = b"ACGT".repeat(250);
        let mut sample_seq = ref_seq.clone();
        sample_seq[500] = b'T';
        let genes = test_genes();

        let mut arena = Arena::new();
        let ref_rope = build_rope_chunked(&mut arena, &ref_seq, 64);
        let sample_rope = build_rope_chunked(&mut arena, &sample_seq, 64);

        let changed = screen_genes_single_arena(
            &mut arena, ref_rope, sample_rope, &genes, None,
        );
        assert_eq!(changed, vec!["GENE_B"]);
    }

    #[test]
    fn test_summary() {
        let ref_seq = b"ACGT".repeat(250);
        let mut sample_seq = ref_seq.clone();
        sample_seq[500] = b'T';
        let genes = test_genes();

        let mut arena = Arena::new();
        let ref_rope = build_rope_chunked(&mut arena, &ref_seq, 64);
        let sample_rope = build_rope_chunked(&mut arena, &sample_seq, 64);

        let report = diff_genes_single_arena(
            &mut arena, ref_rope, sample_rope, &genes, None, true,
        );
        let summary = report.summary();
        assert!(summary.contains("1/3 genes changed"));
        assert!(summary.contains("GENE_B"));
    }
}
