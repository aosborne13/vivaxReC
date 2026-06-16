#!/usr/bin/env python3
"""
summarise_classifications.py

Reads per_infection_classifications TSV files produced by classify_recurrences.R
and produces per-arm summary tables without re-running Pv3Rs.

Reads all available classification files in the output directory:
  per_infection_classifications.tsv              -> table_type: full
  per_infection_classifications_first_rec_only.tsv -> table_type: first_rec
  per_infection_classifications_mono_only.tsv    -> table_type: mono_only

For each file found, produces a separate summary TSV with per-arm rows containing:
  - n_d0_genotyped      : total genotyped D0 patients (constant across all table types)
  - n_no_recurrence     : D0 patients with no classified pair in this subset
  - n_pairs_classified  : pairs with a Pv3Rs result
  - P(L) threshold breakdown (relapse >= 0.75, reinfection <= 0.25, uncertain)
  - Pv3Rs top-call breakdown (L, I, C)
  - Mean and median P(L), P(C), P(I) per arm

Usage:
  python summarise_classifications.py [flags]

Flags:
  --bl_dir      DIR   Baseline classification output directory
  --con_dir     DIR   Consecutive classification output directory
  --master_dir  DIR   Directory containing <trial>/master_samples.tsv (from prep step)
  --trial       STR   Trial name (must match subdirectory in master_dir)
  --out_dir     DIR   Output directory for summary TSVs (default: same as input dir)
"""

import argparse
import pandas as pd
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
    p.add_argument("--bl_dir",     default=None, help="Baseline classification output directory")
    p.add_argument("--con_dir",    default=None, help="Consecutive classification output directory")
    p.add_argument("--master_dir", default=None,
                   help="Directory containing <trial>/master_samples.tsv (from prep step)")
    p.add_argument("--trial",      default=None, help="Trial name for labelling")
    p.add_argument("--out_dir",    default=None,
                   help="Output directory for summary TSVs (default: same as input dir)")
    return p.parse_args()


def load_d0_counts(master_dir, trial):
    """
    Load master_samples.tsv and return a DataFrame with n_d0_genotyped per arm.
    Returns None if master_dir or trial is not provided / file not found.
    """
    if not master_dir or not trial:
        return None
    fp = Path(master_dir) / trial.lower() / "master_samples.tsv"
    if not fp.exists():
        print(f"   WARNING: master_samples.tsv not found at {fp} — n_d0_genotyped will be omitted")
        return None
    ms = pd.read_csv(fp, sep="\t")
    # infection_index == 0 identifies D0 samples
    if "infection_index" not in ms.columns or "treatment_arm" not in ms.columns:
        print(f"   WARNING: master_samples.tsv missing required columns — n_d0_genotyped omitted")
        return None
    d0 = ms[ms["infection_index"] == 0]
    counts = (d0.groupby("treatment_arm", dropna=False)
                .agg(n_d0_genotyped=("patient_id", "nunique"))
                .reset_index())
    return counts


def pl_call(pv3rs_L):
    if pd.isna(pv3rs_L):
        return "uncertain"
    if pv3rs_L >= REL_THRESH:
        return "relapse"
    if pv3rs_L <= REI_THRESH:
        return "reinfection"
    return "uncertain"


def summarise_one(df, pair_mode, table_type, trial=None, d0_counts=None):
    """Produce per-arm summary from a classification dataframe."""
    df = df.copy()

    for col in ["pv3rs_L", "pv3rs_C", "pv3rs_I"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "pv3rs_L" in df.columns:
        df["pl_call"] = df["pv3rs_L"].apply(pl_call)
    else:
        df["pl_call"] = "uncertain"

    # All patient IDs with a classified pair in this file
    classified_pids = set(df.loc[df["pv3rs_L"].notna(), "id"])

    rows = []
    arms = sorted(df["treatment_arm"].dropna().unique())

    for arm in arms:
        grp        = df[df["treatment_arm"] == arm]
        classified = grp.dropna(subset=["pv3rs_L"])
        n_pairs    = len(grp)
        n_class    = len(classified)

        # ── D0 denominator ────────────────────────────────────────────────────
        # n_d0_genotyped is the same regardless of subset — it's the full
        # number of genotyped D0 patients per arm, not the subset count.
        if d0_counts is not None:
            row_d0 = d0_counts[d0_counts["treatment_arm"] == arm]
            n_d0   = int(row_d0["n_d0_genotyped"].iloc[0]) if len(row_d0) else 0
        else:
            n_d0 = None

        # n_no_recurrence: D0 patients in this arm with no classified pair
        if d0_counts is not None and n_d0 is not None:
            arm_d0_pids_approx = n_d0  # we have the count but not the IDs here
            n_with_rec = df[df["treatment_arm"] == arm]["id"].nunique()
            n_no_rec   = max(0, n_d0 - n_with_rec)
        else:
            n_no_rec = None

        # ── P(L) threshold breakdown — applied only to classified pairs ────────
        # "uncertain" here means P(L) genuinely between thresholds (0.25–0.75).
        # Pairs with no pv3rs_L result are reported as n_no_pv3rs, not as uncertain.
        n_no_pv3rs = n_pairs - n_class
        n_rel = (classified["pl_call"] == "relapse").sum()
        n_rei = (classified["pl_call"] == "reinfection").sum()
        n_unc = (classified["pl_call"] == "uncertain").sum()
        pct   = lambda n: round(100 * n / n_class, 1) if n_class > 0 else None

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

        # ── Days to recurrence (baseline first_rec only) ──────────────────────
        mean_days = med_days = None
        if table_type == "first_rec" and pair_mode == "baseline":
            day_col = None
            for c in ["dr_day", "current_day", "days_since_baseline_by_day"]:
                if c in grp.columns:
                    day_col = c
                    break
            if day_col:
                days = pd.to_numeric(grp[day_col], errors="coerce").dropna()
                if len(days) > 0:
                    mean_days = round(float(days.mean()), 1)
                    med_days  = round(float(days.median()), 1)

        rows.append({
            "trial":              trial or "",
            "pair_mode":          pair_mode,
            "table_type":         table_type,
            "treatment_arm":      arm,
            # ── D0 denominator (constant across all table_types) ──────────────
            "n_d0_genotyped":     n_d0,
            "n_no_recurrence":    n_no_rec,
            # ── Pair counts ───────────────────────────────────────────────────
            "n_pairs_total":      n_pairs,
            "n_pairs_classified": n_class,
            "n_no_pv3rs":         n_no_pv3rs,
            # ── Days D0 to first recurrence (first_rec baseline only) ─────────
            "mean_days_to_rec":   mean_days,
            "median_days_to_rec": med_days,
            # ── P(L) threshold breakdown ──────────────────────────────────────
            "n_relapse_pl":       int(n_rel),
            "pct_relapse_pl":     pct(n_rel),
            "n_reinfection_pl":   int(n_rei),
            "pct_reinfection_pl": pct(n_rei),
            "n_uncertain_pl":     int(n_unc),
            "pct_uncertain_pl":   pct(n_unc),
            # ── Pv3Rs top-call breakdown ──────────────────────────────────────
            "n_top_L":            int(n_top_L),
            "pct_top_L":          pct(n_top_L),
            "n_top_I":            int(n_top_I),
            "pct_top_I":          pct(n_top_I),
            "n_top_C":            int(n_top_C),
            "pct_top_C":          pct(n_top_C),
            "n_top_missing":      int(n_top_na),
            # ── Mean / median posteriors ──────────────────────────────────────
            "mean_pv3rs_L":       mean_L,
            "median_pv3rs_L":     med_L,
            "mean_pv3rs_C":       mean_C,
            "median_pv3rs_C":     med_C,
            "mean_pv3rs_I":       mean_I,
            "median_pv3rs_I":     med_I,
        })

    return pd.DataFrame(rows)


def process_dir(input_dir, pair_mode, out_dir, trial, d0_counts):
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

        summ = summarise_one(df, pair_mode=pair_mode, table_type=table_type,
                             trial=trial, d0_counts=d0_counts)

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
    print(f"   trial:      {args.trial or '(not specified)'}")
    print(f"   master_dir: {args.master_dir or '(not provided — n_d0_genotyped will be omitted)'}")

    d0_counts = load_d0_counts(args.master_dir, args.trial)
    if d0_counts is not None:
        print(f"   D0 counts loaded: {len(d0_counts)} arms")

    if args.bl_dir:
        print(f"\n>> baseline: {args.bl_dir}")
        process_dir(args.bl_dir, "baseline", args.out_dir or args.bl_dir,
                    args.trial, d0_counts)

    if args.con_dir:
        print(f"\n>> consecutive: {args.con_dir}")
        process_dir(args.con_dir, "consecutive", args.out_dir or args.con_dir,
                    args.trial, d0_counts)

    print("\n✓ done")


if __name__ == "__main__":
    main()
