#!/usr/bin/env python3
"""E-D2 Drug Resistance Correlation: hashrope detections vs Sierra drug scoring.

Analyzes the correlation between hashrope-detected codon mismatches at
23 resistance panel positions and Sierra's per-drug clinical resistance
levels (Susceptible → High-Level Resistance).

This validates the clinical utility of hashrope's two-pass architecture:
  Pass 1: Hash-based screen identifies codon changes (high sensitivity)
  Pass 2: Annotation classifies clinical significance

If the correlation is strong, it means hashrope's nucleotide-level
detections are predictive of clinically meaningful drug resistance.

Input:
  - data/hivdb_rt_expanded_ground_truth.tsv  — Sierra mutation calls
  - data/hivdb_rt_expanded_drug_scores.tsv   — Sierra per-drug resistance levels

Output: results/drug_resistance_correlation.json (+ timestamped archive)
"""

import json
import sys
import datetime
import platform
from collections import defaultdict
from pathlib import Path

CODE_DIR = Path(__file__).parent.parent
DATA_DIR = CODE_DIR / "data"
RESULTS_DIR = CODE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# Resistance level ordering (Sierra text → numeric rank)
LEVEL_RANK = {
    "Susceptible": 0,
    "Potential Low-Level Resistance": 1,
    "Low-Level Resistance": 2,
    "Intermediate Resistance": 3,
    "High-Level Resistance": 4,
}
LEVEL_LABELS = list(LEVEL_RANK.keys())

# Drug class assignments
NRTI_DRUGS = {"3TC", "ABC", "AZT", "D4T", "DDI", "FTC", "TDF"}
NNRTI_DRUGS = {"DOR", "DPV", "EFV", "ETR", "NVP", "RPV"}


def load_ground_truth(tsv_path: Path) -> dict[str, dict]:
    """Load ground truth TSV → {accession: {panel_mutations, all_rt_mutations}}."""
    truth = {}
    with open(tsv_path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            accession = parts[0].split()[0]
            panel_pos = parts[4].split(",") if parts[4] else []
            panel_txt = parts[5].split(",") if parts[5] else []
            all_muts = parts[6].split(",") if parts[6] else []
            panel_muts = {}
            for p, t in zip(panel_pos, panel_txt):
                panel_muts[int(p)] = t
            truth[accession] = {
                "panel_mutations": panel_muts,
                "panel_count": len(panel_muts),
                "all_rt_mutations": all_muts,
                "total_rt_mutations": len(all_muts),
            }
    return truth


def load_drug_scores(tsv_path: Path) -> dict[str, list[dict]]:
    """Load drug scores TSV → {accession: [{drug_abbr, score, level, text, ...}]}."""
    scores = defaultdict(list)
    with open(tsv_path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            accession = parts[0]
            scores[accession].append({
                "drug_abbr": parts[2],
                "drug_name": parts[3],
                "score": float(parts[4]),
                "level": int(parts[5]),
                "text": parts[6],
                "contributing_mutations": parts[7] if len(parts) > 7 else "",
            })
    return dict(scores)


def main():
    print("=" * 70)
    print("E-D2 Drug Resistance Correlation Analysis")
    print("=" * 70)

    # Load data
    gt_path = DATA_DIR / "hivdb_rt_expanded_ground_truth.tsv"
    ds_path = DATA_DIR / "hivdb_rt_expanded_drug_scores.tsv"
    truth = load_ground_truth(gt_path)
    scores = load_drug_scores(ds_path)
    print(f"Loaded {len(truth)} sequences with ground truth")
    print(f"Loaded {len(scores)} sequences with drug scores")

    # Match sequences
    common = set(truth.keys()) & set(scores.keys())
    print(f"Matched: {len(common)} sequences")

    # === Analysis 1: Panel mutation count vs max resistance level ===
    print(f"\n{'='*70}")
    print("Analysis 1: Panel mutation count → max resistance level")
    print("  (Does more hashrope-detected panel mutations = higher resistance?)")
    print(f"{'='*70}")

    # For each sequence: panel_count, max_resistance_rank, max_resistance_text
    seq_data = []
    for acc in sorted(common):
        gt = truth[acc]
        ds = scores[acc]
        max_rank = max(LEVEL_RANK.get(d["text"], 0) for d in ds)
        max_text = LEVEL_LABELS[max_rank]
        max_nrti_rank = max((LEVEL_RANK.get(d["text"], 0) for d in ds
                             if d["drug_abbr"] in NRTI_DRUGS), default=0)
        max_nnrti_rank = max((LEVEL_RANK.get(d["text"], 0) for d in ds
                              if d["drug_abbr"] in NNRTI_DRUGS), default=0)
        seq_data.append({
            "accession": acc,
            "panel_count": gt["panel_count"],
            "total_rt_mutations": gt["total_rt_mutations"],
            "max_resistance_rank": max_rank,
            "max_resistance_text": max_text,
            "max_nrti_rank": max_nrti_rank,
            "max_nnrti_rank": max_nnrti_rank,
        })

    # Group by panel mutation count
    by_panel_count = defaultdict(list)
    for s in seq_data:
        key = min(s["panel_count"], 4)  # bin 4+ together
        by_panel_count[key].append(s)

    print(f"\n  {'Panel muts':<12} {'N seqs':>8} {'Suscept':>10} {'Pot Low':>10} "
          f"{'Low':>10} {'Interm':>10} {'High':>10} {'% any R':>10}")
    print(f"  {'-'*12:<12} {'-'*8:>8} {'-'*10:>10} {'-'*10:>10} "
          f"{'-'*10:>10} {'-'*10:>10} {'-'*10:>10} {'-'*10:>10}")

    panel_count_summary = []
    for pc in sorted(by_panel_count.keys()):
        seqs = by_panel_count[pc]
        n = len(seqs)
        rank_dist = defaultdict(int)
        for s in seqs:
            rank_dist[s["max_resistance_rank"]] += 1
        pct_resistant = (n - rank_dist[0]) / n * 100 if n > 0 else 0
        label = f"{pc}" if pc < 4 else "4+"
        print(f"  {label:<12} {n:>8} {rank_dist[0]:>10} {rank_dist[1]:>10} "
              f"{rank_dist[2]:>10} {rank_dist[3]:>10} {rank_dist[4]:>10} "
              f"{pct_resistant:>9.1f}%")
        panel_count_summary.append({
            "panel_mutations": label,
            "n_sequences": n,
            "susceptible": rank_dist[0],
            "potential_low": rank_dist[1],
            "low": rank_dist[2],
            "intermediate": rank_dist[3],
            "high": rank_dist[4],
            "pct_any_resistance": round(pct_resistant, 2),
        })

    # === Analysis 2: Resistance level → panel mutation stats ===
    print(f"\n{'='*70}")
    print("Analysis 2: Max resistance level → panel mutation statistics")
    print("  (Do resistant sequences have more hashrope detections?)")
    print(f"{'='*70}")

    by_max_level = defaultdict(list)
    for s in seq_data:
        by_max_level[s["max_resistance_rank"]].append(s)

    print(f"\n  {'Max level':<35} {'N seqs':>8} {'Mean panel':>12} "
          f"{'Median':>8} {'Max':>6} {'Mean total RT':>15}")
    print(f"  {'-'*35:<35} {'-'*8:>8} {'-'*12:>12} "
          f"{'-'*8:>8} {'-'*6:>6} {'-'*15:>15}")

    level_summary = []
    for rank in range(5):
        seqs = by_max_level[rank]
        if not seqs:
            continue
        n = len(seqs)
        panel_counts = sorted(s["panel_count"] for s in seqs)
        total_muts = sorted(s["total_rt_mutations"] for s in seqs)
        mean_panel = sum(panel_counts) / n
        median_panel = panel_counts[n // 2]
        max_panel = panel_counts[-1]
        mean_total = sum(total_muts) / n
        label = LEVEL_LABELS[rank]
        print(f"  {label:<35} {n:>8} {mean_panel:>12.2f} "
              f"{median_panel:>8} {max_panel:>6} {mean_total:>15.1f}")
        level_summary.append({
            "level": label,
            "rank": rank,
            "n_sequences": n,
            "mean_panel_mutations": round(mean_panel, 3),
            "median_panel_mutations": median_panel,
            "max_panel_mutations": max_panel,
            "mean_total_rt_mutations": round(mean_total, 2),
        })

    # === Analysis 3: Per-drug breakdown ===
    print(f"\n{'='*70}")
    print("Analysis 3: Per-drug resistance breakdown")
    print(f"{'='*70}")

    drug_breakdown = {}
    for drug_abbr in sorted(NRTI_DRUGS | NNRTI_DRUGS):
        drug_class = "NRTI" if drug_abbr in NRTI_DRUGS else "NNRTI"
        level_dist = defaultdict(int)
        for acc in common:
            for d in scores[acc]:
                if d["drug_abbr"] == drug_abbr:
                    level_dist[d["text"]] += 1
        total = sum(level_dist.values())
        resistant = total - level_dist.get("Susceptible", 0)
        pct = resistant / total * 100 if total > 0 else 0
        drug_breakdown[drug_abbr] = {
            "class": drug_class,
            "total": total,
            "susceptible": level_dist.get("Susceptible", 0),
            "resistant_any": resistant,
            "pct_resistant": round(pct, 2),
            "high_level": level_dist.get("High-Level Resistance", 0),
            "intermediate": level_dist.get("Intermediate Resistance", 0),
            "low": level_dist.get("Low-Level Resistance", 0),
            "potential_low": level_dist.get("Potential Low-Level Resistance", 0),
        }

    print(f"\n  {'Drug':<6} {'Class':<7} {'Total':>7} {'Suscept':>9} "
          f"{'Pot Low':>9} {'Low':>7} {'Interm':>8} {'High':>7} {'% R':>7}")
    print(f"  {'-'*6:<6} {'-'*7:<7} {'-'*7:>7} {'-'*9:>9} "
          f"{'-'*9:>9} {'-'*7:>7} {'-'*8:>8} {'-'*7:>7} {'-'*7:>7}")
    for drug in sorted(drug_breakdown, key=lambda d: -drug_breakdown[d]["pct_resistant"]):
        db = drug_breakdown[drug]
        print(f"  {drug:<6} {db['class']:<7} {db['total']:>7} {db['susceptible']:>9} "
              f"{db['potential_low']:>9} {db['low']:>7} {db['intermediate']:>8} "
              f"{db['high_level']:>7} {db['pct_resistant']:>6.1f}%")

    # === Analysis 4: Sequences with panel mutations and their drug scores ===
    print(f"\n{'='*70}")
    print("Analysis 4: Hashrope panel detections → per-drug resistance")
    print("  (For sequences WITH vs WITHOUT panel mutations, what's the drug profile?)")
    print(f"{'='*70}")

    with_panel = [s for s in seq_data if s["panel_count"] > 0]
    without_panel = [s for s in seq_data if s["panel_count"] == 0]
    with_accs = {s["accession"] for s in with_panel}
    without_accs = {s["accession"] for s in without_panel}

    print(f"\n  Sequences with panel mutations:    {len(with_panel)}")
    print(f"  Sequences without panel mutations: {len(without_panel)}")

    print(f"\n  {'Drug':<6} {'With panel (% any R)':>22} {'Without panel (% any R)':>25} {'Enrichment':>12}")
    print(f"  {'-'*6:<6} {'-'*22:>22} {'-'*25:>25} {'-'*12:>12}")

    enrichment_data = []
    for drug_abbr in sorted(NRTI_DRUGS | NNRTI_DRUGS):
        # With panel mutations
        with_r = 0
        with_total = 0
        for acc in with_accs:
            for d in scores.get(acc, []):
                if d["drug_abbr"] == drug_abbr:
                    with_total += 1
                    if d["text"] != "Susceptible":
                        with_r += 1
        with_pct = with_r / with_total * 100 if with_total > 0 else 0

        # Without panel mutations
        without_r = 0
        without_total = 0
        for acc in without_accs:
            for d in scores.get(acc, []):
                if d["drug_abbr"] == drug_abbr:
                    without_total += 1
                    if d["text"] != "Susceptible":
                        without_r += 1
        without_pct = without_r / without_total * 100 if without_total > 0 else 0

        enrichment = with_pct / without_pct if without_pct > 0 else float('inf')
        enrichment_str = f"{enrichment:.1f}×" if enrichment != float('inf') else "∞"
        print(f"  {drug_abbr:<6} {with_pct:>20.1f}% {without_pct:>23.1f}% {enrichment_str:>12}")
        enrichment_data.append({
            "drug": drug_abbr,
            "with_panel_pct_resistant": round(with_pct, 2),
            "without_panel_pct_resistant": round(without_pct, 2),
            "enrichment_fold": round(enrichment, 2) if enrichment != float('inf') else None,
        })

    # === Correlation coefficient (Spearman) ===
    print(f"\n{'='*70}")
    print("Analysis 5: Rank correlation (panel mutations vs max resistance)")
    print(f"{'='*70}")

    # Simple Spearman-like: use ranks
    panel_vals = [s["panel_count"] for s in seq_data]
    resist_vals = [s["max_resistance_rank"] for s in seq_data]
    n = len(seq_data)
    # Compute Spearman rank correlation manually (no scipy dependency)
    def rank_data(vals):
        indexed = sorted(enumerate(vals), key=lambda x: x[1])
        ranks = [0.0] * len(vals)
        i = 0
        while i < len(indexed):
            j = i
            while j < len(indexed) and indexed[j][1] == indexed[i][1]:
                j += 1
            avg_rank = (i + j - 1) / 2.0 + 1
            for k in range(i, j):
                ranks[indexed[k][0]] = avg_rank
            i = j
        return ranks

    r_panel = rank_data(panel_vals)
    r_resist = rank_data(resist_vals)
    mean_rp = sum(r_panel) / n
    mean_rr = sum(r_resist) / n
    cov = sum((r_panel[i] - mean_rp) * (r_resist[i] - mean_rr) for i in range(n))
    var_p = sum((r_panel[i] - mean_rp) ** 2 for i in range(n))
    var_r = sum((r_resist[i] - mean_rr) ** 2 for i in range(n))
    spearman_rho = cov / (var_p * var_r) ** 0.5 if var_p > 0 and var_r > 0 else 0
    print(f"\n  Spearman ρ (panel mutation count vs max resistance rank): {spearman_rho:.4f}")
    print(f"  N = {n} sequences")
    if spearman_rho > 0.3:
        print(f"  Interpretation: moderate positive correlation")
    elif spearman_rho > 0.1:
        print(f"  Interpretation: weak positive correlation")
    elif spearman_rho > 0:
        print(f"  Interpretation: very weak positive correlation")
    else:
        print(f"  Interpretation: no positive correlation")
    print(f"  (Note: many sequences are susceptible with 0 panel mutations,")
    print(f"   which compresses the distribution. See Analysis 1 & 4 for detail.)")

    # === Save results ===
    output = {
        "analysis": "drug_resistance_correlation",
        "environment": {
            "timestamp_iso": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "os": platform.system(),
            "python_version": platform.python_version(),
        },
        "config": {
            "dataset": "Expanded 2,000-sequence (NCBI GenBank + Sierra)",
            "sequences_analyzed": len(common),
            "drugs_analyzed": len(NRTI_DRUGS | NNRTI_DRUGS),
            "nrti_drugs": sorted(NRTI_DRUGS),
            "nnrti_drugs": sorted(NNRTI_DRUGS),
        },
        "results": {
            "panel_count_vs_resistance": panel_count_summary,
            "resistance_level_vs_mutations": level_summary,
            "per_drug_breakdown": drug_breakdown,
            "with_vs_without_panel_mutations": {
                "with_panel": len(with_panel),
                "without_panel": len(without_panel),
                "enrichment_per_drug": enrichment_data,
            },
            "spearman_rho": round(spearman_rho, 4),
        },
    }

    out_path = RESULTS_DIR / "drug_resistance_correlation.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved: {out_path}")

    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive = RESULTS_DIR / f"drug_resistance_correlation_{ts}.json"
    with open(archive, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Saved: {archive}")


if __name__ == "__main__":
    main()
