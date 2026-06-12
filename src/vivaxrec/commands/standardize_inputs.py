#!/usr/bin/env python3
"""
Standardise vivaxReC input files and split outputs into:

1. analysis-ready files
   - restricted to samples present in metadata

2. dropped files
   - rows/samples not represented in metadata

Inputs
------
- mhap_between_relatedness_estimate.tsv
- mhap_effective_COI_summary.tsv
- metadata.tsv (any tab-separated metadata file)

Outputs
-------
standardized_inputs/
    analysis_ready/
        metadata.cleaned.tsv
        coi.cleaned.tsv
        relatedness.cleaned.tsv

    dropped/
        relatedness_rows_dropped.tsv
        relatedness_sample_ids_not_in_metadata.tsv
        coi_rows_dropped.tsv

    qc/
        metadata_duplicate_sample_ids.tsv
        coi_duplicate_sample_ids.tsv
        metadata_missing_dates.tsv
        day0_nonD0_observed.tsv
        nonzero_day_with_D0.tsv
"""

import argparse
import pandas as pd
import re
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Standardise vivaxReC input files: filter to metadata universe, "
                    "produce analysis-ready, dropped, and QC outputs."
    )
    p.add_argument("--rel",  default="mhap_between_relatedness_estimate.tsv",
                   help="Pairwise relatedness TSV (default: mhap_between_relatedness_estimate.tsv)")
    p.add_argument("--coi",  default="mhap_effective_COI_summary.tsv",
                   help="MOIRE COI summary TSV (default: mhap_effective_COI_summary.tsv)")
    p.add_argument("--meta", default="metadata.tsv",
                   help="Sample metadata TSV (default: metadata.tsv)")
    p.add_argument("--out_dir", default="standardized_inputs",
                   help="Root output directory (default: standardized_inputs)")
    return p.parse_args()


args     = parse_args()
BASE_DIR = Path(args.out_dir)
ANALYSIS_DIR = BASE_DIR / "analysis_ready"
DROPPED_DIR  = BASE_DIR / "dropped"
QC_DIR       = BASE_DIR / "qc"

for d in [ANALYSIS_DIR, DROPPED_DIR, QC_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# =========================
# HELPERS
# =========================
def clean_sample_id(x):
    if pd.isna(x):
        return x
    x = str(x).strip()
    x = re.sub(r"_D0H0$", "_D0", x)
    return x

def get_suffix(sample_id):
    if pd.isna(sample_id):
        return None
    s = str(sample_id)
    if "_" not in s:
        return None
    return s.rsplit("_", 1)[1]

# =========================
# LOAD
# =========================
print(f">> reading relatedness from: {args.rel}")
rel = pd.read_csv(args.rel, sep="\t")
print(f">> reading COI from:         {args.coi}")
coi = pd.read_csv(args.coi, sep="\t")
print(f">> reading metadata from:    {args.meta}")
meta = pd.read_csv(args.meta, sep="\t")

# Preserve originals for audit
rel["sample_id1_original"] = rel["sample_id1"]
rel["sample_id2_original"] = rel["sample_id2"]
coi["sample_id_original"] = coi["sample_id"]
meta["sample_id_original"] = meta["sample_id"]

# =========================
# CLEAN IDS
# =========================
rel["sample_id1"] = rel["sample_id1"].apply(clean_sample_id)
rel["sample_id2"] = rel["sample_id2"].apply(clean_sample_id)
coi["sample_id"] = coi["sample_id"].apply(clean_sample_id)
meta["sample_id"] = meta["sample_id"].apply(clean_sample_id)

# =========================
# METADATA CLEANING
# =========================
meta["samp_date"] = pd.to_datetime(meta["samp_date"], dayfirst=True, errors="coerce")
meta["day"] = pd.to_numeric(meta["day"], errors="coerce")

for col in ["time_point", "trial", "site", "treatment_arm", "external_id"]:
    meta[col] = meta[col].astype(str).str.strip()

meta["sample_suffix"] = meta["sample_id"].apply(get_suffix)

# Informational flags only
meta["is_day0_nonD0_observed"] = (
    (meta["day"] == 0) & (meta["time_point"] != "D0")
)
meta["is_nonzero_day_D0"] = (
    (meta["day"] > 0) & (meta["time_point"] == "D0")
)

# =========================
# QC TABLES
# =========================
meta_dups = meta[meta.duplicated(subset=["sample_id"], keep=False)].sort_values("sample_id")
coi_dups = coi[coi.duplicated(subset=["sample_id"], keep=False)].sort_values("sample_id")
meta_missing_dates = meta[meta["samp_date"].isna()].copy()
day0_nonD0 = meta[meta["is_day0_nonD0_observed"]].copy()
nonzero_day_D0 = meta[meta["is_nonzero_day_D0"]].copy()

# =========================
# DEFINE ANALYSIS UNIVERSE
# =========================
meta_ids = set(meta["sample_id"])
rel_ids = set(rel["sample_id1"]).union(set(rel["sample_id2"]))
coi_ids = set(coi["sample_id"])

rel_sample_ids_not_in_metadata = sorted(rel_ids - meta_ids)

# =========================
# SPLIT INTO ANALYSIS-READY / DROPPED
# =========================
coi_kept = coi[coi["sample_id"].isin(meta_ids)].copy()
coi_dropped = coi[~coi["sample_id"].isin(meta_ids)].copy()

rel_kept = rel[
    rel["sample_id1"].isin(meta_ids) & rel["sample_id2"].isin(meta_ids)
].copy()

rel_dropped = rel[
    ~(
        rel["sample_id1"].isin(meta_ids) & rel["sample_id2"].isin(meta_ids)
    )
].copy()

# =========================
# WRITE ANALYSIS-READY
# =========================
meta.to_csv(ANALYSIS_DIR / "metadata.cleaned.tsv", sep="\t", index=False)
coi_kept.to_csv(ANALYSIS_DIR / "coi.cleaned.tsv", sep="\t", index=False)
rel_kept.to_csv(ANALYSIS_DIR / "relatedness.cleaned.tsv", sep="\t", index=False)

# =========================
# WRITE DROPPED
# =========================
if not coi_dropped.empty:
    coi_dropped.to_csv(DROPPED_DIR / "coi_rows_dropped.tsv", sep="\t", index=False)

if not rel_dropped.empty:
    rel_dropped.to_csv(DROPPED_DIR / "relatedness_rows_dropped.tsv", sep="\t", index=False)

if rel_sample_ids_not_in_metadata:
    pd.DataFrame({"sample_id_not_in_metadata": rel_sample_ids_not_in_metadata}).to_csv(
        DROPPED_DIR / "relatedness_sample_ids_not_in_metadata.tsv",
        sep="\t",
        index=False
    )

# =========================
# WRITE QC
# =========================
if not meta_dups.empty:
    meta_dups.to_csv(QC_DIR / "metadata_duplicate_sample_ids.tsv", sep="\t", index=False)

if not coi_dups.empty:
    coi_dups.to_csv(QC_DIR / "coi_duplicate_sample_ids.tsv", sep="\t", index=False)

if not meta_missing_dates.empty:
    meta_missing_dates.to_csv(QC_DIR / "metadata_missing_dates.tsv", sep="\t", index=False)

if not day0_nonD0.empty:
    day0_nonD0.to_csv(QC_DIR / "day0_nonD0_observed.tsv", sep="\t", index=False)

if not nonzero_day_D0.empty:
    nonzero_day_D0.to_csv(QC_DIR / "nonzero_day_with_D0.tsv", sep="\t", index=False)

# =========================
# SUMMARY
# =========================
print("===== STANDARDISATION SUMMARY =====")
print(f"Metadata rows: {len(meta)}")
print(f"COI rows (raw): {len(coi)}")
print(f"Relatedness rows (raw): {len(rel)}")

print("\nAnalysis-ready outputs")
print(f"COI rows kept: {len(coi_kept)}")
print(f"Relatedness rows kept: {len(rel_kept)}")

print("\nDropped outputs")
print(f"COI rows dropped: {len(coi_dropped)}")
print(f"Relatedness rows dropped: {len(rel_dropped)}")
print(f"Unique relatedness sample IDs not in metadata: {len(rel_sample_ids_not_in_metadata)}")

print("\nQC")
print(f"Metadata duplicate sample_ids: {len(meta_dups)}")
print(f"COI duplicate sample_ids: {len(coi_dups)}")
print(f"Metadata rows with missing samp_date: {len(meta_missing_dates)}")
print(f"day == 0 but time_point != D0 (informational): {len(day0_nonD0)}")
print(f"day > 0 but time_point == D0 (warning): {len(nonzero_day_D0)}")

print(f"\n✓ Outputs written to: {BASE_DIR.resolve()}")