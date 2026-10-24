# hashrope-bio

Computational biology applications of [hashrope](https://github.com/jemsbhai/hashrope) —
Exact byte comparison for gene-level change detection, genomic region comparison,
mutation localization, drug resistance screening, and tandem repeat compression.

## Exact comparison in 0.3.3

This release requires Hashrope 0.3.3 or later in the 0.3 series. Fingerprint
values/defaults are unchanged; equality, mutation and conservation decisions now
verify actual bytes. Report hash fields are metadata, and gene/exon reports carry
separate verified change flags. Manually constructing these reports requires
those flags rather than deriving them from hashes.

Exact verification can scan an entire independent region and uses a per-call
comparison cache where needed. Shared ranges and periodic structure can be skipped after
proof. Historical speedups and datasets below describe the earlier fingerprint
experiments; they do not establish current exact-comparison performance or
universal collision resistance. Historical result files have been preserved.

The [exact-path optimization report](../../../benchmarks/exact_optimization/README.md)
reruns the full stored HIV comparison workload against archived preoptimization
and optimized cores. It includes verification costs and construction exclusions;
it does not revalidate the older chromosome-scale fingerprint speedups below.

## Install

```bash
pip install "hashrope-bio==0.3.3"
```

For optional backends:
```bash
pip install -e ".[genomics]"   # pysam + biopython for FASTA I/O
```

## Features

### Gene-level change detection (`genomics.gene_diff`)

Compare gene bodies and exons exactly, with shortcuts for shared and repeated
structure. Independent equal regions may require reading every byte.

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

Check N known resistance mutation sites by comparing actual codon bytes.

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

Exact byte comparison for arbitrary genomic regions; fingerprint queries remain available separately.

```python
from hashrope_bio.genomics import region_hash, regions_identical

identical = regions_identical(ref_rope, sample_rope, start=10000, length=50000, h=h)
```

### Mutation localization (`genomics.mutation`)

Exact prefix comparison localizes the first divergent byte. Matching fingerprints do not prune mutations.

Validated on 3,422 ClinVar pathogenic SNVs on chr22: **3,422/3,422 correct (100%)**.

### Tandem repeat compression (`genomics.repeats`)

RepeatNode gives O(log q) hash computation for q-repeat motifs vs O(q·d) naive.
**610× speedup at q=10,000** for clinical trinucleotide repeats (HTT, FMR1).

## Project structure

```
src/hashrope_bio/
├── genomics/          # FASTA loading, region queries, mutation search, gene diff, repeats
├── proteomics/        # MD trajectory frame comparison (stub)
└── cheminformatics/   # Resistance panels, compound lookup with retained verification bytes
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

## Compound lookup compatibility

`build_compound_index` retains every candidate and its canonical encoding bytes.
`lookup_compound` checks those bytes inside the matching hash bucket; collisions
and duplicate input rows are retained. Legacy plain dictionaries containing only
hashes/row IDs are rejected because they cannot verify identity; rebuild them
from the source file. Canonicalization remains the caller's responsibility.

Run consumer tests from this directory with local Hashrope installed (or add
`../../../packages/python` and `src` to `PYTHONPATH`): `python -m pytest tests`.
