# Changelog

## 0.3.4 — 2026-10-24

- Require core Hashrope 0.3.4 and synchronize with the Python biology package.
- Keep exact consumer algorithms, fingerprint fields and byte semantics
  unchanged; update benchmark provenance to name the new dependency version.
- Document [release parity and performance caveats](https://github.com/jemsbhai/hashrope/blob/master/releases/0.3.4/README.md), including unresolved core comparison costs.

## 0.3.3 — 2026-10-24

- Require Hashrope 0.3.3 and synchronize with the Python package.
- Replace fingerprint-only gene/exon and resistance decisions with exact byte
  range comparisons, avoiding temporary split nodes and substring rehashing.
- Preserve fingerprint report fields and add verified change flags. Manual
  report construction must provide those flags instead of inferring changes
  from hashes.
- Add deliberate-collision regressions and byte-oracle checks. Exact independent
  comparisons can scan the full region; historical fingerprint speedups are
  not claimed for these operations.
