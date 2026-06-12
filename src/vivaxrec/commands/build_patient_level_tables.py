#!/usr/bin/env python3
"""
Build patient-level recurrence tables from standardized vivaxReC inputs.

Inputs
------
standardized_inputs/analysis_ready/
    metadata.cleaned.tsv
    coi.cleaned.tsv
    relatedness.cleaned.tsv

Outputs
-------
patient_level_tables/
    baseline_samples.tsv
    patient_infection_sequence.tsv
    consecutive_pairs.tsv
    consecutive_pairs_with_relatedness.tsv
    baseline_pairs_with_relatedness.tsv
    summary_by_trial_site.tsv
"""

import argparse
import pandas as pd
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Build patient-level recurrence tables from standardised vivaxReC inputs."
    )
    p.add_argument("--analysis_dir", default=str(Path("standardized_inputs") / "analysis_ready"),
                   help="Directory containing cleaned analysis-ready TSVs "
                        "(default: standardized_inputs/analysis_ready)")
    p.add_argument("--out_dir", default="patient_level_tables",
                   help="Output directory (default: patient_level_tables)")
    return p.parse_args()


args     = parse_args()
BASE_DIR = Path(args.analysis_dir)
META_FILE = BASE_DIR / "metadata.cleaned.tsv"
COI_FILE  = BASE_DIR / "coi.cleaned.tsv"
REL_FILE  = BASE_DIR / "relatedness.cleaned.tsv"

OUT_DIR = Path(args.out_dir)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# =========================
# LOAD
# =========================
print(f">> reading metadata from:    {META_FILE}")
meta = pd.read_csv(META_FILE, sep="\t")
print(f">> reading COI from:         {COI_FILE}")
coi = pd.read_csv(COI_FILE, sep="\t")
print(f">> reading relatedness from: {REL_FILE}")
rel = pd.read_csv(REL_FILE, sep="\t")

# =========================
# CLEAN TYPES
# =========================
meta["day"] = pd.to_numeric(meta["day"], errors="coerce")
meta["samp_date"] = pd.to_datetime(meta["samp_date"], dayfirst=True, errors="coerce")

for col in ["sample_id", "external_id", "trial", "site", "time_point", "treatment_arm"]:
    if col in meta.columns:
        meta[col] = meta[col].astype(str).str.strip()

if "sample_id" in coi.columns:
    coi["sample_id"] = coi["sample_id"].astype(str).str.strip()

for col in ["sample_id1", "sample_id2"]:
    if col in rel.columns:
        rel[col] = rel[col].astype(str).str.strip()

# =========================
# COI: one row per sample_id
# =========================
coi_keep_cols = [c for c in [
    "sample_id",
    "post_effective_coi_lower",
    "post_effective_coi_med",
    "post_effective_coi_upper",
    "post_effective_coi_mean"
] if c in coi.columns]

coi = coi[coi_keep_cols].drop_duplicates(subset=["sample_id"]).copy()

# =========================
# MERGE COI INTO METADATA
# =========================
meta = meta.merge(coi, on="sample_id", how="left")

# =========================
# BUILD INFECTION SEQUENCE
# =========================
meta = meta.sort_values(
    by=["trial", "site", "external_id", "day", "samp_date", "sample_id"],
    na_position="last"
).copy()

# infection_index = first observed infection = 0, then 1,2,3...
meta["infection_index"] = meta.groupby("external_id").cumcount()

# baseline sample = first observed infection per patient
baseline = meta[meta["infection_index"] == 0].copy()

# add previous infection info for sequential comparisons
meta["prev_sample_id"] = meta.groupby("external_id")["sample_id"].shift(1)
meta["prev_day"] = meta.groupby("external_id")["day"].shift(1)
meta["prev_samp_date"] = meta.groupby("external_id")["samp_date"].shift(1)
meta["prev_time_point"] = meta.groupby("external_id")["time_point"].shift(1)
meta["prev_infection_index"] = meta.groupby("external_id")["infection_index"].shift(1)

meta["days_since_prev_sample_by_date"] = (
    meta["samp_date"] - meta["prev_samp_date"]
).dt.days

meta["days_since_prev_sample_by_day"] = meta["day"] - meta["prev_day"]

# full patient infection sequence
patient_seq = meta.copy()

# consecutive pairs only
pairs = meta.dropna(subset=["prev_sample_id"]).copy()

# =========================
# PREP RELATEDNESS
# =========================
# map sample_id -> patient/trial/site from metadata
sample_map = meta[["sample_id", "external_id", "trial", "site"]].drop_duplicates()

rel = rel.merge(
    sample_map.rename(columns={
        "sample_id": "sample_id1",
        "external_id": "external_id1",
        "trial": "trial1",
        "site": "site1"
    }),
    on="sample_id1",
    how="left"
)

rel = rel.merge(
    sample_map.rename(columns={
        "sample_id": "sample_id2",
        "external_id": "external_id2",
        "trial": "trial2",
        "site": "site2"
    }),
    on="sample_id2",
    how="left"
)

# keep only within-patient comparisons
rel_within = rel[
    rel["external_id1"].notna() &
    rel["external_id2"].notna() &
    (rel["external_id1"] == rel["external_id2"])
].copy()

# symmetrize so pair order does not matter when merging
rel_rev = rel_within.rename(columns={
    "sample_id1": "sample_id2",
    "sample_id2": "sample_id1",
    "coi1": "coi2",
    "coi2": "coi1",
    "external_id1": "external_id2",
    "external_id2": "external_id1",
    "trial1": "trial2",
    "trial2": "trial1",
    "site1": "site2",
    "site2": "site1"
}).copy()

rel_sym = pd.concat([rel_within, rel_rev], ignore_index=True)
rel_sym = rel_sym.drop_duplicates(subset=["sample_id1", "sample_id2"])

rel_keep_cols = [c for c in [
    "sample_id1", "sample_id2",
    "M1", "M", "rtotal",
    "coi1", "coi2",
    "scaled_r", "relatedness"
] if c in rel_sym.columns]

rel_sym = rel_sym[rel_keep_cols].copy()

# =========================
# SEQUENTIAL COMPARISONS
# D0->DRP1, DRP1->DRP2, DRP2->DRP3, ...
# =========================
pairs_rel = pairs.merge(
    rel_sym,
    left_on=["prev_sample_id", "sample_id"],
    right_on=["sample_id1", "sample_id2"],
    how="left"
)

pairs_rel = pairs_rel.rename(columns={
    "sample_id": "current_sample_id",
    "day": "current_day",
    "samp_date": "current_samp_date",
    "time_point": "current_time_point",
    "post_effective_coi_lower": "current_post_effective_coi_lower",
    "post_effective_coi_med": "current_post_effective_coi_med",
    "post_effective_coi_upper": "current_post_effective_coi_upper",
    "post_effective_coi_mean": "current_post_effective_coi_mean",
    "prev_sample_id": "previous_sample_id",
    "prev_day": "previous_day",
    "prev_samp_date": "previous_samp_date",
    "prev_time_point": "previous_time_point"
})

for cutoff in [0.25, 0.5, 0.95]:
    label = str(cutoff).replace(".", "_")
    pairs_rel[f"related_ge_{label}"] = pairs_rel["relatedness"] >= cutoff
    pairs_rel[f"unrelated_lt_{label}"] = pairs_rel["relatedness"] < cutoff

pairs_rel["calendar_year"] = pairs_rel["current_samp_date"].dt.year
pairs_rel["calendar_month"] = pairs_rel["current_samp_date"].dt.month
pairs_rel["calendar_month_name"] = pairs_rel["current_samp_date"].dt.month_name()

# =========================
# BASELINE COMPARISONS
# D0->DRP1, D0->DRP2, D0->DRP3, ...
# more generally: first observed infection -> all later infections
# =========================
baseline_map = (
    meta.sort_values(["external_id", "day", "samp_date", "sample_id"])
    .groupby("external_id")
    .first()
    .reset_index()[["external_id", "sample_id", "day", "samp_date", "time_point"]]
    .rename(columns={
        "sample_id": "baseline_sample_id",
        "day": "baseline_day",
        "samp_date": "baseline_samp_date",
        "time_point": "baseline_time_point"
    })
)

baseline_pairs = meta.merge(baseline_map, on="external_id", how="left")

# remove self-comparisons (baseline vs baseline)
baseline_pairs = baseline_pairs[
    baseline_pairs["sample_id"] != baseline_pairs["baseline_sample_id"]
].copy()

baseline_pairs_rel = baseline_pairs.merge(
    rel_sym,
    left_on=["baseline_sample_id", "sample_id"],
    right_on=["sample_id1", "sample_id2"],
    how="left"
)

baseline_pairs_rel = baseline_pairs_rel.rename(columns={
    "sample_id": "current_sample_id",
    "day": "current_day",
    "samp_date": "current_samp_date",
    "time_point": "current_time_point",
    "post_effective_coi_lower": "current_post_effective_coi_lower",
    "post_effective_coi_med": "current_post_effective_coi_med",
    "post_effective_coi_upper": "current_post_effective_coi_upper",
    "post_effective_coi_mean": "current_post_effective_coi_mean"
})

baseline_pairs_rel["days_since_baseline_by_day"] = (
    baseline_pairs_rel["current_day"] - baseline_pairs_rel["baseline_day"]
)

baseline_pairs_rel["days_since_baseline_by_date"] = (
    baseline_pairs_rel["current_samp_date"] - baseline_pairs_rel["baseline_samp_date"]
).dt.days

for cutoff in [0.25, 0.5, 0.95]:
    label = str(cutoff).replace(".", "_")
    baseline_pairs_rel[f"related_ge_{label}"] = baseline_pairs_rel["relatedness"] >= cutoff
    baseline_pairs_rel[f"unrelated_lt_{label}"] = baseline_pairs_rel["relatedness"] < cutoff

baseline_pairs_rel["calendar_year"] = baseline_pairs_rel["current_samp_date"].dt.year
baseline_pairs_rel["calendar_month"] = baseline_pairs_rel["current_samp_date"].dt.month
baseline_pairs_rel["calendar_month_name"] = baseline_pairs_rel["current_samp_date"].dt.month_name()

# =========================
# SUMMARIES
# =========================
summary = (
    meta.groupby(["trial", "site"], dropna=False)
    .agg(
        n_patients=("external_id", "nunique"),
        n_samples=("sample_id", "count"),
        n_baseline_samples=("infection_index", lambda x: (x == 0).sum())
    )
    .reset_index()
)

seq_summary = (
    pairs_rel.groupby(["trial", "site"], dropna=False)
    .agg(
        n_sequential_pairs=("current_sample_id", "count"),
        n_sequential_pairs_with_relatedness=("relatedness", lambda x: x.notna().sum()),
        mean_sequential_relatedness=("relatedness", "mean"),
        median_sequential_relatedness=("relatedness", "median")
    )
    .reset_index()
)

base_summary = (
    baseline_pairs_rel.groupby(["trial", "site"], dropna=False)
    .agg(
        n_baseline_pairs=("current_sample_id", "count"),
        n_baseline_pairs_with_relatedness=("relatedness", lambda x: x.notna().sum()),
        mean_baseline_relatedness=("relatedness", "mean"),
        median_baseline_relatedness=("relatedness", "median")
    )
    .reset_index()
)

summary = summary.merge(seq_summary, on=["trial", "site"], how="left")
summary = summary.merge(base_summary, on=["trial", "site"], how="left")

# =========================
# WRITE OUTPUTS
# =========================
baseline.to_csv(OUT_DIR / "baseline_samples.tsv", sep="\t", index=False)
patient_seq.to_csv(OUT_DIR / "patient_infection_sequence.tsv", sep="\t", index=False)
pairs.to_csv(OUT_DIR / "consecutive_pairs.tsv", sep="\t", index=False)
pairs_rel.to_csv(OUT_DIR / "consecutive_pairs_with_relatedness.tsv", sep="\t", index=False)
baseline_pairs_rel.to_csv(OUT_DIR / "baseline_pairs_with_relatedness.tsv", sep="\t", index=False)
summary.to_csv(OUT_DIR / "summary_by_trial_site.tsv", sep="\t", index=False)

# =========================
# CONSOLE SUMMARY
# =========================
print("===== PATIENT-LEVEL TABLE SUMMARY =====")
print(f"Baseline samples: {len(baseline)}")
print(f"Patient infection rows: {len(patient_seq)}")
print(f"Sequential pairs: {len(pairs_rel)}")
print(f"Sequential pairs with relatedness found: {pairs_rel['relatedness'].notna().sum()}")
print(f"Baseline pairs: {len(baseline_pairs_rel)}")
print(f"Baseline pairs with relatedness found: {baseline_pairs_rel['relatedness'].notna().sum()}")
print(f"\n✓ Output directory: {OUT_DIR.resolve()}")