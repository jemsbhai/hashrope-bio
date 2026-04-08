# hashrope-bio

Computational biology applications of [hashrope](https://github.com/jemsbhai/hashrope) —
O(log w) substring hashing for genomic region comparison, mutation localization,
tandem repeat compression, MD trajectory indexing, and drug resistance screening.

## Install

```bash
pip install -e ".[genomics,bench,dev]"
```

## Project structure

```
src/hashrope_bio/
├── genomics/          # FASTA loading, region queries, mutation search, tandem repeats
├── proteomics/        # MD trajectory frame comparison, periodic detection
└── cheminformatics/   # SMILES/InChI compound lookup, resistance panels
benchmarks/            # Experiment scripts (E-G1 through E-CG1)
scripts/               # Dataset download helpers
data/                  # Downloaded datasets (gitignored)
tests/                 # Unit tests
```

## Experiments

See `../../experiments/EXPERIMENTS.md` for the full protocol.

Run a benchmark:
```bash
python benchmarks/bench_region_query.py
```

## Depends on

- [hashrope](https://pypi.org/project/hashrope/) (PyPI) — core data structure
- pysam / biopython — FASTA I/O
- MDAnalysis — trajectory I/O
- RDKit — cheminformatics

## License

MIT
