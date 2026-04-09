#!/usr/bin/env python3
"""E-D2 Real-Data Validation: hashrope resistance panel vs Sierra ground truth.

Validates that hash-based codon comparison correctly identifies drug
resistance mutations in real HIV-1 patient sequences, using Stanford
HIVDB Sierra as gold-standard ground truth.

Pipeline:
  1. Extract HXB2 RT reference (1680 bp) from genome
  2. Load Sierra-aligned patient RT sequences
  3. For each patient: build hashrope, run 23-position panel
  4. Compare hashrope detections vs Sierra mutation calls
  5. Cross-validate: byte comparison must agree with hashrope

Input files (in data/):
  - hiv1_hxb2.fa              — HXB2 full genome (9719 bp)
  - hivdb_rt_aligned.fa        — Sierra-aligned patient RT sequences
  - hivdb_rt_ground_truth.tsv  — Sierra mutation calls

Usage:
  python bench_resistance_realdata.py              # original 100 sequences
  python bench_resistance_realdata.py --expanded   # expanded 2,000+ sequences

Output: results/resistance_realdata[_expanded].json (+ timestamped archive)
"""

import json
import sys
import time
from dataclasses import dataclass, asdict
from pathlib import Path

try:
    from hashrope import PolynomialHash, Leaf, rope_concat, rope_substr_hash
except ImportError:
    print("ERROR: hashrope package required. Install with: pip install hashrope")
    sys.exit(1)

# --- Configuration ---
CODE_DIR = Path(__file__).parent.parent
DATA_DIR = CODE_DIR / "data"
RESULTS_DIR = CODE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

HXB2_RT_START = 2549  # 0-based inclusive
HXB2_RT_END = 4229    # 0-based exclusive
HXB2_RT_LEN = 1680    # 560 codons

# Our 23-position resistance panel
PANEL_POSITIONS = sorted([
    # NNRTI (12)
    100, 101, 103, 106, 108, 138, 181, 188, 190, 225, 227, 230,
    # NRTI (11)
    41, 65, 67, 70, 74, 115, 151, 184, 210, 215, 219,
])

CHUNK_SIZE = 64  # small gene, small chunks

# Standard genetic code (codon -> amino acid)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Known DRMs at our panel positions (position -> set of mutant AAs)
# From Stanford HIVDB
KNOWN_DRMS = {
    41:  {'L'},         # M41L
    65:  {'R', 'E'},    # K65R, K65E
    67:  {'N', 'G', 'E'},  # D67N/G/E
    70:  {'R', 'E', 'G'},  # K70R/E/G
    74:  {'V', 'I'},    # L74V/I
    100: {'I'},         # L100I
    101: {'E', 'P', 'H'},  # K101E/P/H
    103: {'N', 'S', 'H'},  # K103N/S/H
    106: {'A', 'M'},    # V106A/M
    108: {'I'},         # V108I
    115: {'F'},         # Y115F
    138: {'A', 'G', 'K', 'Q', 'R'},  # E138A/G/K/Q/R
    151: {'M'},         # Q151M
    181: {'C', 'I', 'V'},  # Y181C/I/V
    184: {'V', 'I'},    # M184V/I
    188: {'L', 'H', 'C'},  # Y188L/H/C
    190: {'A', 'S', 'E', 'Q'},  # G190A/S/E/Q
    210: {'W'},         # L210W
    215: {'Y', 'F'},    # T215Y/F
    219: {'Q', 'E', 'N', 'R'},  # K219Q/E/N/R
    225: {'H'},         # P225H
    227: {'L', 'C'},    # F227L/C
    230: {'L', 'I'},    # M230L/I
}


def translate_codon(codon: str) -> str:
    """Translate a DNA codon to amino acid. Returns '?' for ambiguous codons."""
    codon = codon.upper()
    return CODON_TABLE.get(codon, '?')


def classify_mismatch(pos: int, ref_codon: str, sample_codon: str) -> dict:
    """Classify a codon mismatch using the annotation pass.

    Returns dict with:
      - ref_aa, sample_aa: translated amino acids
      - is_synonymous: same AA despite different codon
      - is_known_drm: sample AA is in KNOWN_DRMS for this position
      - classification: 'synonymous' | 'known_drm' | 'polymorphism' | 'ambiguous'
      - annotation: human-readable string
    """
    ref_aa = translate_codon(ref_codon)
    sample_aa = translate_codon(sample_codon)

    if sample_aa == '?':
        return {
            "ref_aa": ref_aa, "sample_aa": sample_aa,
            "is_synonymous": False, "is_known_drm": False,
            "classification": "ambiguous",
            "annotation": f"{ref_aa}{pos}? (ambiguous codon: {sample_codon})",
        }

    if ref_aa == sample_aa:
        return {
            "ref_aa": ref_aa, "sample_aa": sample_aa,
            "is_synonymous": True, "is_known_drm": False,
            "classification": "synonymous",
            "annotation": f"{ref_aa}{pos}{sample_aa} (synonymous: {ref_codon}\u2192{sample_codon})",
        }

    drm_aas = KNOWN_DRMS.get(pos, set())
    if sample_aa in drm_aas:
        return {
            "ref_aa": ref_aa, "sample_aa": sample_aa,
            "is_synonymous": False, "is_known_drm": True,
            "classification": "known_drm",
            "annotation": f"{ref_aa}{pos}{sample_aa} (DRUG RESISTANCE MUTATION)",
        }

    return {
        "ref_aa": ref_aa, "sample_aa": sample_aa,
        "is_synonymous": False, "is_known_drm": False,
        "classification": "polymorphism",
        "annotation": f"{ref_aa}{pos}{sample_aa} (natural polymorphism)",
    }

def _build_rope(data: bytes, chunk_size: int, h: PolynomialHash):
    """Build a chunked rope from bytes."""
    rope = None
    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        leaf = Leaf(chunk, h)
        rope = rope_concat(rope, leaf, h)
    return rope


def load_fasta(path: Path) -> list[tuple[str, str]]:
    """Parse FASTA file into (header, sequence) pairs."""
    entries = []
    header = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_lines)))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        entries.append((header, "".join(seq_lines)))
    return entries


def extract_hxb2_rt(genome_path: Path) -> bytes:
    """Extract the RT coding region from HXB2 genome."""
    entries = load_fasta(genome_path)
    if not entries:
        raise RuntimeError(f"No sequences in {genome_path}")
    _, genome_seq = entries[0]
    genome_seq = genome_seq.upper()

    rt_seq = genome_seq[HXB2_RT_START:HXB2_RT_END]
    if len(rt_seq) != HXB2_RT_LEN:
        raise RuntimeError(
            f"Expected RT length {HXB2_RT_LEN}, got {len(rt_seq)}. "
            f"Check HXB2 coordinates."
        )
    print(f"  HXB2 RT first 3 codons: {rt_seq[:9]}")
    print(f"  HXB2 RT length: {len(rt_seq)} bp")
    return rt_seq.encode("ascii")


def load_ground_truth(tsv_path: Path) -> dict[str, dict]:
    """Load Sierra ground truth TSV.

    Returns dict of accession -> {
        'subtype': str,
        'rt_first_aa': int,
        'rt_last_aa': int,
        'panel_mutations': {position: text},
        'all_rt_mutations': [text],
    }
    """
    truth = {}
    with open(tsv_path) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            full_header = parts[0]
            accession = full_header.split()[0]

            panel_pos = parts[4].split(",") if parts[4] else []
            panel_txt = parts[5].split(",") if parts[5] else []
            all_muts = parts[6].split(",") if parts[6] else []

            panel_muts = {}
            for pos_str, txt in zip(panel_pos, panel_txt):
                panel_muts[int(pos_str)] = txt

            truth[accession] = {
                "subtype": parts[1],
                "rt_first_aa": int(parts[2]),
                "rt_last_aa": int(parts[3]),
                "panel_mutations": panel_muts,
                "all_rt_mutations": all_muts,
            }
    return truth


def validate_sequence(
    ref_bytes: bytes,
    sample_seq: str,
    sample_coverage: tuple[int, int],  # (first_aa, last_aa)
    h: PolynomialHash,
) -> dict:
    """Run hashrope panel + byte baseline on one patient sequence.

    Returns dict with per-position results.
    """
    first_aa, last_aa = sample_coverage

    # Convert sample to bytes, uppercase
    sample_upper = sample_seq.upper()
    sample_bytes = sample_upper.encode("ascii")

    # Build ropes
    ref_offset = (first_aa - 1) * 3
    ref_end = min(last_aa * 3, len(ref_bytes))
    ref_region = ref_bytes[ref_offset:ref_end]

    expected_len = (last_aa - first_aa + 1) * 3
    if len(sample_bytes) < expected_len:
        pass

    ref_rope = _build_rope(ref_region, CHUNK_SIZE, h)
    sample_rope = _build_rope(sample_bytes[:len(ref_region)], CHUNK_SIZE, h)

    results = {}
    for pos in PANEL_POSITIONS:
        if pos < first_aa or pos > last_aa:
            results[pos] = {"status": "out_of_range"}
            continue

        codon_offset_in_region = (pos - first_aa) * 3

        if codon_offset_in_region + 3 > len(ref_region) or \
           codon_offset_in_region + 3 > len(sample_bytes):
            results[pos] = {"status": "truncated"}
            continue

        ref_hash = rope_substr_hash(ref_rope, codon_offset_in_region, 3, h)
        sample_hash = rope_substr_hash(sample_rope, codon_offset_in_region, 3, h)
        hashrope_mismatch = (ref_hash != sample_hash)

        ref_codon = ref_region[codon_offset_in_region:codon_offset_in_region + 3]
        sample_codon = sample_bytes[codon_offset_in_region:codon_offset_in_region + 3]
        byte_mismatch = (ref_codon != sample_codon)

        results[pos] = {
            "status": "checked",
            "ref_codon": ref_codon.decode("ascii"),
            "sample_codon": sample_codon.decode("ascii"),
            "hashrope_mismatch": hashrope_mismatch,
            "byte_mismatch": byte_mismatch,
            "hashrope_agrees_with_bytes": hashrope_mismatch == byte_mismatch,
        }

    return results


def main(expanded=False):
    import platform
    import datetime

    # File prefix: original (100 seqs) or expanded (2,000+ seqs)
    prefix = "hivdb_rt_expanded_" if expanded else "hivdb_rt_"
    result_tag = "resistance_realdata_expanded" if expanded else "resistance_realdata"
    label = "Expanded (2,000+)" if expanded else "Original (100)"

    print("=" * 60)
    print(f"E-D2 Real-Data Validation: hashrope vs HIVDB Sierra [{label}]")
    print("=" * 60)

    # Step 1: Load HXB2 RT reference
    print("\nLoading HXB2 RT reference...")
    hxb2_path = DATA_DIR / "hiv1_hxb2.fa"
    if not hxb2_path.exists():
        print(f"ERROR: {hxb2_path} not found. Run download_data.py first.")
        sys.exit(1)
    ref_bytes = extract_hxb2_rt(hxb2_path)

    # Step 2: Load aligned patient sequences
    print("\nLoading Sierra-aligned patient sequences...")
    aligned_path = DATA_DIR / f"{prefix}aligned.fa"
    patient_seqs = load_fasta(aligned_path)
    print(f"  Loaded {len(patient_seqs)} sequences")

    # Step 3: Load Sierra ground truth
    print("\nLoading Sierra ground truth...")
    tsv_path = DATA_DIR / f"{prefix}ground_truth.tsv"
    truth = load_ground_truth(tsv_path)
    print(f"  Loaded truth for {len(truth)} sequences")

    # Step 4: Initialize hashrope
    h = PolynomialHash()

    # Step 5: Validate each sequence
    print("\nRunning validation...")
    all_results = []
    total_checked = 0
    total_hashrope_byte_agree = 0
    total_hashrope_byte_disagree = 0

    sierra_positive_total = 0
    sierra_positive_detected = 0
    sierra_positive_missed = 0
    extra_detections = 0
    true_negatives = 0

    for hdr, seq in patient_seqs:
        accession = hdr.split()[0]
        if accession not in truth:
            continue

        gt = truth[accession]
        coverage = (gt["rt_first_aa"], gt["rt_last_aa"])

        results = validate_sequence(ref_bytes, seq, coverage, h)

        seq_record = {
            "accession": accession,
            "subtype": gt["subtype"],
            "coverage": f"RT:{coverage[0]}-{coverage[1]}",
            "sierra_panel_mutations": gt["panel_mutations"],
            "positions": {},
        }

        for pos in PANEL_POSITIONS:
            r = results[pos]
            if r["status"] != "checked":
                seq_record["positions"][pos] = r
                continue

            total_checked += 1
            sierra_flagged = pos in gt["panel_mutations"]
            hashrope_detected = r["hashrope_mismatch"]

            if r["hashrope_agrees_with_bytes"]:
                total_hashrope_byte_agree += 1
            else:
                total_hashrope_byte_disagree += 1

            if sierra_flagged and hashrope_detected:
                sierra_positive_detected += 1
                sierra_positive_total += 1
                classification = "true_positive"
            elif sierra_flagged and not hashrope_detected:
                sierra_positive_missed += 1
                sierra_positive_total += 1
                classification = "FALSE_NEGATIVE"  # bad!
            elif not sierra_flagged and hashrope_detected:
                extra_detections += 1
                classification = "extra_detection"  # polymorphism
            else:
                true_negatives += 1
                classification = "true_negative"

            r["sierra_flagged"] = sierra_flagged
            r["classification"] = classification

            # Annotation pass: classify mismatches
            if hashrope_detected:
                annot = classify_mismatch(pos, r["ref_codon"], r["sample_codon"])
                r["annotation"] = annot

            seq_record["positions"][pos] = r

        all_results.append(seq_record)

    # Step 6: Compute metrics
    sensitivity = (sierra_positive_detected / sierra_positive_total * 100
                   if sierra_positive_total > 0 else 0.0)
    specificity_denom = true_negatives + extra_detections
    specificity = (true_negatives / specificity_denom * 100
                   if specificity_denom > 0 else 0.0)

    # Annotation summary
    annot_counts = {"synonymous": 0, "known_drm": 0, "polymorphism": 0, "ambiguous": 0}
    for rec in all_results:
        for pos, r in rec["positions"].items():
            if isinstance(r, dict) and "annotation" in r:
                cls = r["annotation"]["classification"]
                annot_counts[cls] = annot_counts.get(cls, 0) + 1
    total_mismatches = sum(annot_counts.values())

    print(f"\n{'=' * 60}")
    print("Results:")
    print(f"  Sequences validated: {len(all_results)}")
    print(f"  Codon positions checked: {total_checked}")
    print(f"  Hashrope-byte agreement: {total_hashrope_byte_agree}/{total_checked} "
          f"({total_hashrope_byte_agree/total_checked*100:.1f}%)")
    print(f"  Hashrope-byte DISAGREEMENT: {total_hashrope_byte_disagree}")
    print()
    print(f"  Sierra DRMs at panel positions: {sierra_positive_total}")
    print(f"    Detected by hashrope (TP): {sierra_positive_detected}")
    print(f"    Missed by hashrope (FN): {sierra_positive_missed}")
    print(f"    Sensitivity: {sensitivity:.1f}%")
    print()
    print(f"  Hashrope extra detections (polymorphisms): {extra_detections}")
    print(f"  True negatives: {true_negatives}")
    print(f"  Specificity (vs Sierra): {specificity:.1f}%")
    print()
    print(f"  --- Annotation Pass (second pass on {total_mismatches} mismatches) ---")
    print(f"    Known DRMs:     {annot_counts['known_drm']}")
    print(f"    Polymorphisms:  {annot_counts['polymorphism']}")
    print(f"    Synonymous:     {annot_counts['synonymous']}")
    print(f"    Ambiguous:      {annot_counts['ambiguous']}")

    if sierra_positive_missed > 0:
        print("\n  *** FALSE NEGATIVES DETECTED — investigate! ***")
        for rec in all_results:
            for pos, r in rec["positions"].items():
                if isinstance(r, dict) and r.get("classification") == "FALSE_NEGATIVE":
                    print(f"    {rec['accession']} pos={pos} "
                          f"ref={r['ref_codon']} sample={r['sample_codon']} "
                          f"sierra={rec['sierra_panel_mutations'].get(pos, '?')}")
    print(f"{'=' * 60}")

    # Step 7: Save results
    output = {
        "benchmark": "resistance_panel_realdata",
        "environment": {
            "timestamp_iso": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "os": platform.system(),
            "arch": platform.machine(),
            "python_version": platform.python_version(),
            "language": "python",
        },
        "config": {
            "reference": "HXB2 RT (K03455, genome positions 2550-4229)",
            "patient_sequences": len(all_results),
            "panel_positions": PANEL_POSITIONS,
            "chunk_size": CHUNK_SIZE,
            "data_source": "NCBI GenBank + Stanford HIVDB Sierra GraphQL",
            "dataset": label,
            "sierra_citation": "Liu TF, Shafer RW. Web Resources for HIV Type 1 Genotypic-Resistance Test Interpretation. Clin Infect Dis. 2006;42(11):1608-1618.",
        },
        "results": {
            "sequences_validated": len(all_results),
            "codon_positions_checked": total_checked,
            "hashrope_byte_agreement": total_hashrope_byte_agree,
            "hashrope_byte_disagreement": total_hashrope_byte_disagree,
            "sierra_drm_total": sierra_positive_total,
            "true_positives": sierra_positive_detected,
            "false_negatives": sierra_positive_missed,
            "sensitivity_pct": round(sensitivity, 2),
            "extra_detections_polymorphisms": extra_detections,
            "true_negatives": true_negatives,
            "specificity_vs_sierra_pct": round(specificity, 2),
            "annotation_pass": {
                "total_mismatches_annotated": total_mismatches,
                "known_drms": annot_counts["known_drm"],
                "polymorphisms": annot_counts["polymorphism"],
                "synonymous": annot_counts["synonymous"],
                "ambiguous": annot_counts["ambiguous"],
            },
        },
        "per_sequence": all_results,
    }

    # Save latest
    out_path = RESULTS_DIR / f"{result_tag}.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {out_path}")

    # Save timestamped archive
    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive_path = RESULTS_DIR / f"{result_tag}_{ts}.json"
    with open(archive_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Saved: {archive_path}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="E-D2 resistance panel validation")
    parser.add_argument("--expanded", action="store_true",
                        help="Use expanded 2,000+ sequence dataset")
    args = parser.parse_args()
    main(expanded=args.expanded)
