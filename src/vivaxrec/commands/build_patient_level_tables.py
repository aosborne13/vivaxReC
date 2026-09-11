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
    excluded_patients_chain_validity.tsv   [NEW]

── Valid consecutive chain logic ───────────────────────────────────────────────
A patient's samples are only used for pairwise analysis (consecutive pairs AND
baseline-to-recurrence pairs) up to and including the longest unbroken run of
recurrences starting from a TRUE Day 0 (time_point == "D0", i.e. recurrence
number 0). Recurrence number is parsed from time_point ("D0" -> 0, "DRPn" -> n),
NOT inferred from row position or infection_index, since infection_index is
purely positional and does not detect gaps.

Rules:
  - If a patient has no true D0 sample at all, they are excluded from every
    pairwise analysis entirely (no substitution of a later recurrence as a
    stand-in baseline).
  - Walking forward from D0 (0, 1, 2, ...), the first missing recurrence number
    truncates the chain. Every sample at or before the truncation point is
    retained ("in_valid_chain" = True); everything after is discarded, even if
    it was successfully genotyped.
  - A patient with D0 only (no recurrences at all) is retained (in_valid_chain
    = True for that one row) but will not appear in any pair, since pairing
    requires two timepoints.
  - A patient with D0 followed immediately by a gap (e.g. D0, [R1 missing], R2)
    retains only their D0 row; R2 is discarded.

patient_infection_sequence.tsv retains ALL rows (nothing is deleted from the
record) with two new columns added for transparency:
  - recurrence_num           parsed integer recurrence number (D0=0, DRPn=n)
  - in_valid_chain           True if this row is within the patient's longest
                              valid, gap-free, D0-anchored prefix
  - max_valid_recurrence_num the cutoff recurrence number for this patient
                              (-1 if the patient has no true D0 at all)

All pair-building (consecutive_pairs*, baseline_pairs*) is restricted to rows
where in_valid_chain == True, so gapped/orphaned samples can never enter a pair.
excluded_patients_chain_validity.tsv lists every patient affected, for audit.
"""

import argparse
import re
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

# infection_index = first observed infection = 0, then 1,2,3... (POSITIONAL —
# retained for backwards compatibility / other reporting, but NOT used to
# detect gaps; see recurrence_num below for that.)
meta["infection_index"] = meta.groupby("external_id").cumcount()

# ── True recurrence number, parsed from time_point (NOT positional) ──────────
# "D0" -> 0, "DRP1" -> 1, "DRP2" -> 2, etc. Anything unparseable -> NaN, which
# will correctly break a chain (treated as a gap) rather than being silently
# skipped.
_DRP_RE = re.compile(r'^DRP(\d+)$')

def parse_recurrence_num(tp):
    if pd.isna(tp):
        return None
    tp = str(tp).strip().upper()
    if tp == "D0":
        return 0
    m = _DRP_RE.match(tp)
    return int(m.group(1)) if m else None

meta["recurrence_num"] = meta["time_point"].apply(parse_recurrence_num)

# ── Valid consecutive chain cutoff per patient ────────────────────────────────
def valid_chain_cutoff(nums):
    """
    nums: sorted, deduplicated list of recurrence_num values observed for one
    patient (NaNs already dropped).
    Returns the highest recurrence_num that is part of an unbroken 0,1,2,...
    run starting at 0. Returns -1 if no true D0 (recurrence_num==0) is present
    at all.
    """
    if not nums or nums[0] != 0:
        return -1
    expected = 0
    max_valid = -1
    for n in nums:
        if n == expected:
            max_valid = n
            expected += 1
        elif n > expected:
            break  # gap detected — stop here, everything after is discarded
        # n < expected only happens on an already-deduplicated, sorted list
        # if there were duplicate recurrence numbers for this patient, which
        # would itself be a data-quality issue; harmless to skip here.
    return max_valid

chain_cutoff = (
    meta.dropna(subset=["recurrence_num"])
        .groupby("external_id")["recurrence_num"]
        .apply(lambda s: valid_chain_cutoff(sorted(s.unique())))
        .rename("max_valid_recurrence_num")
)

meta = meta.merge(chain_cutoff, on="external_id", how="left")
# Patients with NO parseable recurrence_num at all (e.g. every time_point was
# unparseable) get max_valid_recurrence_num = NaN here; treat as fully invalid.
meta["max_valid_recurrence_num"] = meta["max_valid_recurrence_num"].fillna(-1)

meta["in_valid_chain"] = (
    meta["recurrence_num"].notna() &
    (meta["max_valid_recurrence_num"] >= 0) &
    (meta["recurrence_num"] <= meta["max_valid_recurrence_num"])
)

# ── Audit output: which patients/samples were excluded, and why ──────────────
excluded = meta[~meta["in_valid_chain"]].copy()
if not excluded.empty:
    excluded["exclusion_reason"] = excluded.apply(
        lambda r: "no true D0 for this patient" if r["max_valid_recurrence_num"] == -1
                  else "orphaned: beyond first gap in consecutive chain",
        axis=1
    )
    excluded[[
        "external_id", "trial", "site", "sample_id", "time_point",
        "recurrence_num", "max_valid_recurrence_num", "exclusion_reason"
    ]].sort_values(["trial", "external_id", "recurrence_num"]).to_csv(
        OUT_DIR / "excluded_patients_chain_validity.tsv", sep="\t", index=False
    )

n_patients_total    = meta["external_id"].nunique()
n_patients_no_d0    = meta.loc[meta["max_valid_recurrence_num"] == -1, "external_id"].nunique()
n_samples_orphaned  = ((~meta["in_valid_chain"]) & (meta["max_valid_recurrence_num"] >= 0)).sum()
print(f">> chain validity: {n_patients_no_d0} / {n_patients_total} patients have no true D0 "
      f"(excluded entirely); {n_samples_orphaned} further sample(s) orphaned beyond a gap")

# add previous infection info for sequential comparisons
# (still computed on the FULL meta for informational purposes in
# patient_infection_sequence.tsv; pair-building below restricts to
# in_valid_chain rows so gaps can never leak into a pair)
meta["prev_sample_id"] = meta.groupby("external_id")["sample_id"].shift(1)
meta["prev_day"] = meta.groupby("external_id")["day"].shift(1)
meta["prev_samp_date"] = meta.groupby("external_id")["samp_date"].shift(1)
meta["prev_time_point"] = meta.groupby("external_id")["time_point"].shift(1)
meta["prev_infection_index"] = meta.groupby("external_id")["infection_index"].shift(1)

meta["days_since_prev_sample_by_date"] = (
    meta["samp_date"] - meta["prev_samp_date"]
).dt.days

meta["days_since_prev_sample_by_day"] = meta["day"] - meta["prev_day"]

# full patient infection sequence — ALL rows retained, in_valid_chain flags
# which ones are usable for pairwise analysis
patient_seq = meta.copy()

# baseline sample = first observed TRUE D0 per patient (infection_index==0
# no longer used here — recurrence_num==0 is the authoritative definition)
baseline = meta[meta["recurrence_num"] == 0].copy()

# ── Restrict to valid-chain rows before building ANY pairs ───────────────────
meta_valid = meta[meta["in_valid_chain"]].copy()

# consecutive pairs only, built from the valid-chain-restricted data.
# Because meta_valid contains exactly each patient's unbroken 0,1,2,...
# prefix, a plain positional shift(1) is now safe again — there are no
# gaps left to accidentally bridge.
pairs = meta_valid.dropna(subset=["prev_sample_id"]).copy()
# Guard: prev_sample_id must ALSO be in_valid_chain (it always will be, given
# meta_valid only contains valid-chain rows and shift(1) is computed within
# each patient's already-restricted group) — kept as an explicit check rather
# than an assumption.
pairs = pairs[pairs["prev_sample_id"].isin(meta_valid["sample_id"])].copy()

# =========================
# PREP RELATEDNESS
# =========================
# map sample_id -> patient/trial/site from metadata (full meta, not just
# valid-chain — relatedness lookups themselves are harmless for any sample;
# it's PAIR CONSTRUCTION that must be restricted, which happens below)
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
# (built only from the gap-free, D0-anchored prefix — see meta_valid above)
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
# Restricted to meta_valid: only a patient's TRUE D0 sample is used as
# baseline_sample_id (no substitution), and only samples within that
# patient's valid, gap-free chain are paired against it.
# =========================
baseline_map = (
    meta_valid[meta_valid["recurrence_num"] == 0]
    .sort_values(["external_id", "day", "samp_date", "sample_id"])
    .drop_duplicates(subset="external_id", keep="first")
    [["external_id", "sample_id", "day", "samp_date", "time_point"]]
    .rename(columns={
        "sample_id": "baseline_sample_id",
        "day": "baseline_day",
        "samp_date": "baseline_samp_date",
        "time_point": "baseline_time_point"
    })
)

# inner join: only patients with a valid baseline_map entry (i.e. a true D0
# that is itself in_valid_chain, which it always is when present) contribute
# baseline pairs at all.
baseline_pairs = meta_valid.merge(baseline_map, on="external_id", how="inner")

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
        n_baseline_samples=("recurrence_num", lambda x: (x == 0).sum())
    )
    .reset_index()
)

chain_summary = (
    meta.groupby(["trial", "site"], dropna=False)
    .agg(
        n_samples_in_valid_chain=("in_valid_chain", "sum"),
        n_samples_excluded=("in_valid_chain", lambda x: (~x).sum())
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

summary = summary.merge(chain_summary, on=["trial", "site"], how="left")
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
print(f"Patient infection rows (all, incl. excluded): {len(patient_seq)}")
print(f"  of which in valid chain: {meta['in_valid_chain'].sum()}")
print(f"Sequential pairs (valid-chain only): {len(pairs_rel)}")
print(f"Sequential pairs with relatedness found: {pairs_rel['relatedness'].notna().sum()}")
print(f"Baseline pairs (valid-chain only): {len(baseline_pairs_rel)}")
print(f"Baseline pairs with relatedness found: {baseline_pairs_rel['relatedness'].notna().sum()}")
if not excluded.empty:
    print(f"\n⚠ {excluded['external_id'].nunique()} patient(s) affected by chain-validity exclusion "
          f"— see excluded_patients_chain_validity.tsv")
print(f"\n✓ Output directory: {OUT_DIR.resolve()}")
