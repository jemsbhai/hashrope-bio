# hashrope-bio

Computational biology applications of [hashrope](https://github.com/jemsbhai/hashrope) —
O(log w) substring hashing for gene-level change detection, genomic region comparison,
mutation localization, drug resistance screening, and tandem repeat compression.

## Install

```bash
pip install hashrope-bio
```

For optional backends:
```bash
pip install hashrope-bio[genomics]   # pysam + biopython for FASTA I/O
```

## Features

### Gene-level change detection (`genomics.gene_diff`)

Screen hundreds of genes between two genomes in O(G · log w) — orders of magnitude
faster than byte-by-byte comparison for large genes.

```python
from hashrope_bio.genomics import load_fasta_to_rope, load_gene_regions, diff_genes

ref_rope, h, _ = load_fasta_to_rope("reference.fa")
sample_rope, _, _ = load_fasta_to_rope("patient.fa", h=h)
genes = load_gene_regions("chr22_genes.tsv")
report = diff_genes(ref_rope, sample_rope, genes, h)
print(report.summary())
# Gene diff: 12/643 genes changed, 47/5127 exons changed
#   LARGE1 (650,315 bp): body changed, 3/14 exons differ
#   ...
```

Validated on GRCh38 chr22 with 643 real genes (UCSC refGene): **5,770/5,770 queries correct**.
Speedups: **25× for typical genes (10–50 Kbp), 411× for large genes (>500 Kbp)**.

### Drug resistance screening (`cheminformatics.resistance`)

Check N known resistance mutation sites via codon hash comparison.

```python
from hashrope_bio.cheminformatics.resistance import check_resistance_panel, HIV_RT_NNRTI_PANEL

results = check_resistance_panel(ref_rope, sample_rope, HIV_RT_NNRTI_PANEL, h)
for r in results:
    if r.is_mutant:
        print(f"  {r.site.annotation}")
```

Validated on 100 real clinical HIV-1 sequences against Stanford HIVDB:
**103/103 drug resistance mutations detected (100% sensitivity)**.

### Region identity queries (`genomics.region_query`)

O(log w) hash comparison for arbitrary genomic regions.

```python
from hashrope_bio.genomics import region_hash, regions_identical

identical = regions_identical(ref_rope, sample_rope, start=10000, length=50000, h=h)
```

### Mutation localization (`genomics.mutation`)

Binary search via substr_hash localizes a single-nucleotide change in O(log N) comparisons.

Validated on 3,422 ClinVar pathogenic SNVs on chr22: **3,422/3,422 correct (100%)**.

### Tandem repeat compression (`genomics.repeats`)

RepeatNode gives O(log q) hash computation for q-repeat motifs vs O(q·d) naive.
**610× speedup at q=10,000** for clinical trinucleotide repeats (HTT, FMR1).

## Project structure

```
src/hashrope_bio/
├── genomics/          # FASTA loading, region queries, mutation search, gene diff, repeats
├── proteomics/        # MD trajectory frame comparison (stub)
└── cheminformatics/   # Resistance panels, compound lookup (stub)
benchmarks/            # Experiment scripts
scripts/               # Dataset download helpers
data/                  # Downloaded datasets (gitignored)
tests/                 # Unit tests (28 passing)
```

## Depends on

- [hashrope](https://pypi.org/project/hashrope/) — core rope data structure
- pysam / biopython (optional) — FASTA I/O

## Citation

If you use hashrope-bio in published research, please cite:
- Repository: https://github.com/jemsbhai/hashrope-bio
- Core data structure: https://github.com/jemsbhai/hashrope

## License

MIT
