# Changelog

## 0.3.3 — 2026-10-24

- Require core Hashrope 0.3.3 and synchronize Python/Rust package versions.
- Verify actual bytes for genomic regions, mutations, resistance codons, gene
  bodies and exons, trajectory frames, and compound candidates. Matching
  fingerprints no longer decide equality or skip differences.
- Use exact range APIs without splitting or rehashing candidate ranges. Shared
  and periodic structure can be skipped only after proof; independent equal
  regions may require a full scan.
- Retain all colliding compound candidates and verify their stored bytes.
- Gene/exon report change flags now record exact comparison outcomes; manual
  report construction must supply those flags. Fingerprints remain metadata.
- Add collision regression tests and distinguish historical fingerprint-only
  benchmark claims from current exact-query performance.
