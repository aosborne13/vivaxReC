#!/usr/bin/env python3
"""
summarise_classifications.py

Reads per_infection_classifications TSV files produced by classify_recurrences.R
and produces per-arm summary tables without re-running Pv3Rs.

Reads all available classification files in the output directory:
  per_infection_classifications.tsv             -> table_type: full
  per_infection_classifications_first_rec.tsv   -> table_type: first_rec
  per_infection_classifications_mono_only.tsv   -> table_type: mono_only

For each file found, produces a separate summary TSV:
  summary_by_arm.tsv
  summary_by_arm_first_rec.tsv
  summary_by_arm_mono_only.tsv

Each summary contains two breakdowns per treatment arm:

1. P(L) threshold classification (agreed analytical framework):
     relapse      : P(L) >= 0.75
     reinfection  : P(L) <= 0.25
     uncertain    : 0.25 < P(L) < 0.75 or missing

2. Pv3Rs top-call breakdown (pv3rs_top after recrudescence reassignment):
     L (relapse), I (reinfection), C (recrudescence)

3. Mean and median P(L), P(C), P(I) per arm.

Usage:
  python summarise_classifications.py [flags]

Flags:
  --bl_dir   DIR   Baseline classification output directory
  --con_dir  DIR   Consecutive classification output directory
  --out_dir  DIR   Directory to write summary TSVs (default: same as input dir)
  --trial    STR   Trial name for labelling (optional)
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path

REL_THRESH = 0.75
REI_THRESH = 0.25

CLASSIFICATION_FILES = {
    "full":       "per_infection_classifications.tsv",
    "first_rec":  "per_infection_classifications_first_rec_only.tsv",
    "mono_only":  "per_infection_classifications_mono_only.tsv",
}

SUMMARY_FILES = {
    "full":       "summary_by_arm.tsv",
    "first_rec":  "summary_by_arm_first_rec.tsv",
    "mono_only":  "summary_by_arm_mono_only.tsv",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bl_dir",  default=None, help="Baseline classification output directory")
    p.add_argument("--con_dir", default=None, help="Consecutive classification output directory")
    p.add_argument("--out_dir", default=None,
                   help="Output directory for summary TSVs (default: same as input dir)")
    p.add_argument("--trial",   default=None, help="Trial name for labelling")
    return p.parse_args()


def pl_call(pv3rs_L):
    """Apply P(L) threshold classification."""
    if pd.isna(pv3rs_L):
        return "uncertain"
    if pv3rs_L >= REL_THRESH:
        return "relapse"
    if pv3rs_L <= REI_THRESH:
        return "reinfection"
    return "uncertain"


def summarise_one(df, pair_mode, table_type, trial=None):
    """Produce per-arm summary from a classification dataframe."""

    df = df.copy()

    # Coerce numeric columns
    for col in ["pv3rs_L", "pv3rs_C", "pv3rs_I"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # P(L) threshold classification
    if "pv3rs_L" in df.columns:
        df["pl_call"] = df["pv3rs_L"].apply(pl_call)
    else:
        df["pl_call"] = "uncertain"

    rows = []
    for arm, grp in df.groupby("treatment_arm", dropna=False):
        classified = grp.dropna(subset=["pv3rs_L"])
        n_pairs    = len(grp)
        n_class    = len(classified)

        # ── P(L) threshold breakdown ──────────────────────────────────────────
        n_rel   = (grp["pl_call"] == "relapse").sum()
        n_rei   = (grp["pl_call"] == "reinfection").sum()
        n_unc   = (grp["pl_call"] == "uncertain").sum()
        pct = lambda n: round(100 * n / n_class, 1) if n_class > 0 else None

        # ── Pv3Rs top-call breakdown ──────────────────────────────────────────
        n_top_L = n_top_I = n_top_C = n_top_na = 0
        if "pv3rs_top" in grp.columns:
            n_top_L  = (grp["pv3rs_top"] == "L").sum()
            n_top_I  = (grp["pv3rs_top"] == "I").sum()
            n_top_C  = (grp["pv3rs_top"] == "C").sum()
            n_top_na = grp["pv3rs_top"].isna().sum()

        # ── Mean / median posteriors ──────────────────────────────────────────
        def stat(col):
            if col not in classified.columns:
                return None, None
            v = classified[col].dropna()
            if len(v) == 0:
                return None, None
            return round(float(v.mean()), 4), round(float(v.median()), 4)

        mean_L, med_L = stat("pv3rs_L")
        mean_C, med_C = stat("pv3rs_C")
        mean_I, med_I = stat("pv3rs_I")

        rows.append({
            "trial":                  trial or "",
            "pair_mode":              pair_mode,
            "table_type":             table_type,
            "treatment_arm":          arm,
            # ── pair counts ───────────────────────────────────────────────────
            "n_pairs_total":          n_pairs,
            "n_pairs_classified":     n_class,
            # ── P(L) threshold breakdown ──────────────────────────────────────
            "n_relapse_pl":           int(n_rel),
            "pct_relapse_pl":         pct(n_rel),
            "n_reinfection_pl":       int(n_rei),
            "pct_reinfection_pl":     pct(n_rei),
            "n_uncertain_pl":         int(n_unc),
            "pct_uncertain_pl":       pct(n_unc),
            # ── Pv3Rs top-call breakdown ──────────────────────────────────────
            "n_top_L":                int(n_top_L),
            "pct_top_L":              pct(n_top_L),
            "n_top_I":                int(n_top_I),
            "pct_top_I":              pct(n_top_I),
            "n_top_C":                int(n_top_C),
            "pct_top_C":              pct(n_top_C),
            "n_top_missing":          int(n_top_na),
            # ── Mean / median posteriors ──────────────────────────────────────
            "mean_pv3rs_L":           mean_L,
            "median_pv3rs_L":         med_L,
            "mean_pv3rs_C":           mean_C,
            "median_pv3rs_C":         med_C,
            "mean_pv3rs_I":           mean_I,
            "median_pv3rs_I":         med_I,
        })

    return pd.DataFrame(rows)


def process_dir(input_dir, pair_mode, out_dir, trial):
    input_dir = Path(input_dir)
    out_dir   = Path(out_dir) if out_dir else input_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    found = 0
    for table_type, fname in CLASSIFICATION_FILES.items():
        fp = input_dir / fname
        if not fp.exists():
            print(f"   [{pair_mode}] {fname} not found — skipping {table_type}")
            continue

        print(f"   [{pair_mode}] reading {fname} ...")
        df = pd.read_csv(fp, sep="\t")
        print(f"   [{pair_mode}] {len(df)} rows, {df['treatment_arm'].nunique()} arms")

        summ = summarise_one(df, pair_mode=pair_mode, table_type=table_type, trial=trial)

        out_fp = out_dir / SUMMARY_FILES[table_type]
        summ.to_csv(out_fp, sep="\t", index=False)
        print(f"   [{pair_mode}] wrote: {out_fp}")
        found += 1

    if found == 0:
        print(f"   [{pair_mode}] no classification files found in {input_dir}")


def main():
    args = parse_args()

    if not args.bl_dir and not args.con_dir:
        print("ERROR: at least one of --bl_dir or --con_dir is required.")
        raise SystemExit(1)

    print(f"\n>> summarise_classifications")
    print(f"   trial:   {args.trial or '(not specified)'}")

    if args.bl_dir:
        print(f"\n>> baseline: {args.bl_dir}")
        process_dir(args.bl_dir, "baseline", args.out_dir or args.bl_dir, args.trial)

    if args.con_dir:
        print(f"\n>> consecutive: {args.con_dir}")
        process_dir(args.con_dir, "consecutive", args.out_dir or args.con_dir, args.trial)

    print("\n✓ done")


if __name__ == "__main__":
    main()
