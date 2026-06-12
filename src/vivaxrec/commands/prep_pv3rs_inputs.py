#!/usr/bin/env python3
# =============================================================================
# prep_pv3rs_inputs.py
#
# Prepares genotype and metadata tables required by classify_recurrences.R (Pv3Rs pipeline).
#
# Standardised inputs (same files used across the whole pipeline):
#   patient_infection_sequence.tsv   -- all samples, all trials; carries MOIRE
#                                       eMOI (post_effective_coi_mean), metadata,
#                                       and patient linkage via external_id
#   mhap_filtered.tsv                -- single combined long-format genotype file
#                                       (sample_id, locus, allele[, count])
#                                       covering all trials; trial is
#                                       assigned by joining on sample_id from pis
#
# Outputs written to <out_dir>/ and <out_dir>/<trial>/ for each trial found in the data:
#   master_samples.tsv            <- read by classify_recurrences.R
#   per_sample_locus_alleles.tsv
#   d0_dr_pairs.tsv
#   allele_freq_from_D0.tsv
#
# Usage:
#   python prep_pv3rs_inputs.py [flags]
#
# Flags (all optional; defaults shown):
#   --pis       FILE   patient_infection_sequence.tsv       [patient_infection_sequence.tsv]
#   --geno      FILE   combined mhap_filtered.tsv           [mhap_filtered.tsv]
#   --out_dir   DIR    root output directory                 [vivaxrec_inputs]
#
# Full example with explicit paths:
#   python prep_pv3rs_inputs.py \
#     --pis     /data/inputs/patient_level_tables/patient_infection_sequence.tsv \
#     --geno    /data/inputs/standardized_inputs/mhap_filtered.tsv \
#     --out_dir /data/output-analysis/recurrence-classification
# =============================================================================

import argparse
import pandas as pd
import re
from pathlib import Path
from collections import Counter


# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Prepare Pv3Rs-ready input TSVs from standardised pipeline inputs."
    )
    p.add_argument(
        "--pis", default="patient_infection_sequence.tsv",
        help="patient_infection_sequence.tsv — carries sample metadata, MOIRE eMOI "
             "(post_effective_coi_mean), and patient linkage (external_id)"
    )
    p.add_argument(
        "--geno", default="mhap_filtered.tsv",
        help="Combined genotype file covering both trials (sample_id, locus, allele[, count])"
    )
    p.add_argument(
        "--out_dir", default="vivaxrec_inputs",
        help="Root output directory"
    )
    return p.parse_args()


# ── Helpers ───────────────────────────────────────────────────────────────────
def safe_read(path: Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.read_csv(path)


def collapse_to_sets(df: pd.DataFrame) -> pd.DataFrame:
    """Long (sample_id, locus, allele) -> per sample x locus unique allele set ('|' joined)."""
    def uniq_join(series):
        vals = [str(x) for x in series.dropna()]
        return "|".join(sorted({v for v in vals if v and v.lower() != "nan"}))
    return (
        df.groupby(["sample_id", "locus"], as_index=False)
          .agg(alleles_set=("allele", uniq_join))
    )


def alleles_count(alleles_set_str: str) -> int:
    if not isinstance(alleles_set_str, str) or alleles_set_str.lower() in ("", "nan"):
        return 0
    return len([a for a in alleles_set_str.split("|") if a])


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args    = parse_args()
    geno_fp = Path(args.geno)
    out_dir = Path(args.out_dir)

    # ── patient_infection_sequence ────────────────────────────────────────────
    # Single source of truth for sample metadata, MOIRE eMOI, and patient linkage.
    # Replaces the old separate metadata file and both trial-specific COI files.
    # external_id is the authoritative patient identifier -- do not re-derive it
    # from sample_id (the string formats differ).
    print(f">> reading patient_infection_sequence from: {args.pis}")
    pis = safe_read(args.pis)
    need_pis = ["sample_id", "external_id", "trial", "treatment_arm", "site",
                "day", "time_point", "post_effective_coi_mean", "infection_index"]
    missing = [c for c in need_pis if c not in pis.columns]
    if missing:
        raise RuntimeError(
            f"patient_infection_sequence.tsv is missing columns: {missing}\n"
            f"  Present: {list(pis.columns)}"
        )

    pis = pis.copy()
    pis["day"]   = pd.to_numeric(pis["day"], errors="coerce")
    pis["trial"] = pis["trial"].str.upper().str.strip()

    # Build MOIRE eMOI lookup (same source as classify_recurrences.R uses)
    emoi_map = (
        pis[pis["post_effective_coi_mean"].notna()]
           .drop_duplicates(subset="sample_id")
           .set_index("sample_id")["post_effective_coi_mean"]
           .to_dict()
    )
    n_emoi_missing = pis["post_effective_coi_mean"].isna().sum()
    if n_emoi_missing:
        print(f"   ({n_emoi_missing} samples have no MOIRE eMOI; "
              f"naive allele-count LB will be used for those)")

    # ── Raw genotypes ─────────────────────────────────────────────────────────
    # Single combined file covering both trials. Trial is assigned by joining
    # on sample_id from patient_infection_sequence (already loaded above).
    print(f">> reading genotypes from: {geno_fp}")
    geno_raw = safe_read(geno_fp)
    if geno_raw.empty:
        raise RuntimeError(f"Genotype file not found or empty: {geno_fp}")
    missing = [c for c in ["sample_id", "locus", "allele"] if c not in geno_raw.columns]
    if missing:
        raise RuntimeError(
            f"Genotype file missing columns {missing}. Present: {list(geno_raw.columns)}"
        )

    # Build sample_id -> trial map from patient_infection_sequence
    sid_to_trial = (
        pis.drop_duplicates(subset="sample_id")
           .set_index("sample_id")["trial"]
           .to_dict()
    )

    # Collapse to per sample x locus allele sets, then assign trial
    geno_sets = collapse_to_sets(geno_raw)
    geno_sets["trial"] = geno_sets["sample_id"].map(sid_to_trial)

    n_no_trial = geno_sets["trial"].isna().sum()
    if n_no_trial:
        print(f"   ({n_no_trial} genotype rows have no matching sample in "
              f"patient_infection_sequence and will be dropped)")
    geno = geno_sets.dropna(subset=["trial"]).copy()
    print(f"   genotype samples: {geno['sample_id'].nunique()} "
          f"across {geno['trial'].nunique()} trial(s): "
          f"{', '.join(sorted(geno['trial'].unique()))}")

    # ── master_samples ────────────────────────────────────────────────────────
    # All metadata comes from patient_infection_sequence.
    # external_id is used as patient_id throughout (not derived from sample_id).
    meta_cols = ["sample_id", "external_id", "trial", "treatment_arm",
                 "site", "day", "time_point", "infection_index"]
    pis_meta = (pis[meta_cols]
                .drop_duplicates(subset="sample_id")
                .rename(columns={"external_id": "patient_id"})
                .copy())

    # Start from genotyped samples and join metadata
    samples = (
        geno[["trial", "sample_id"]]
        .drop_duplicates()
        .merge(pis_meta, on=["trial", "sample_id"], how="left")
    )

    # Naive COI lower bound: max allele count over loci per sample
    # Used only when MOIRE eMOI is unavailable for a sample
    coi_lb = (
        geno.groupby(["trial", "sample_id"])["alleles_set"]
            .apply(lambda s: max((alleles_count(x) for x in s), default=0))
            .reset_index()
            .rename(columns={"alleles_set": "coi_naive_lb"})
    )
    samples = samples.merge(coi_lb, on=["trial", "sample_id"], how="left")

    # MOIRE eMOI is the primary COI estimate, consistent with classify_recurrences.R
    samples["post_effective_coi_mean"] = samples["sample_id"].map(emoi_map)
    samples["final_coi"] = samples["post_effective_coi_mean"].fillna(
        samples["coi_naive_lb"]
    )
    samples["coi_source"] = samples["post_effective_coi_mean"].apply(
        lambda x: "MOIRE_post_effective_coi_mean" if pd.notna(x) else "naive_lower_bound"
    )

    # ── per_sample_locus_alleles ──────────────────────────────────────────────
    per_sample_locus = geno.merge(
        samples[["trial", "sample_id", "patient_id",
                 "treatment_arm", "site", "day", "time_point",
                 "infection_index", "final_coi", "coi_source"]],
        on=["trial", "sample_id"], how="left"
    )

    # ── D0→DR pairs ───────────────────────────────────────────────────────────
    # infection_index==0 identifies the baseline sample per patient.
    # This correctly handles the small number of cases where the baseline
    # sample is not at day 0 (flagged as is_nonzero_day_D0 in the source data).
    d0_map = (
        pis[pis["infection_index"] == 0]
           .drop_duplicates(subset=["external_id", "trial"])
           [["external_id", "trial", "sample_id"]]
           .rename(columns={"sample_id": "d0_sample_id",
                            "external_id": "patient_id"})
    )
    dr_rows = (
        pis[pis["infection_index"] > 0]
           [["external_id", "trial", "sample_id", "day", "time_point"]]
           .rename(columns={"sample_id": "dr_sample_id",
                            "external_id": "patient_id",
                            "day": "dr_day",
                            "time_point": "dr_time_point"})
    )
    pairs_df = (
        d0_map.merge(dr_rows, on=["patient_id", "trial"], how="inner")
              .sort_values(["trial", "patient_id", "dr_day"])
              .reset_index(drop=True)
    )

    # ── Allele frequencies from baseline (infection_index==0) samples ─────────
    d0_sids = set(d0_map["d0_sample_id"])
    d0_sub  = per_sample_locus[per_sample_locus["sample_id"].isin(d0_sids)]

    pi_rows = []
    for (trial_val, locus), sub in d0_sub.groupby(["trial", "locus"]):
        n = sub["sample_id"].nunique()
        counts = Counter()
        for s in sub["alleles_set"].astype(str):
            alleles = [a for a in s.split("|") if a]
            counts.update(set(alleles))
        for a, cnt in counts.items():
            pi_rows.append({
                "trial":     trial_val,
                "locus":     locus,
                "allele":    a,
                "freq":      cnt / n if n else 0.0,
                "n_samples": n
            })
    pi_df = pd.DataFrame(pi_rows).sort_values(["trial", "locus", "allele"])

    # ── Write outputs ─────────────────────────────────────────────────────────
    out_dir.mkdir(parents=True, exist_ok=True)
    samples.to_csv(          out_dir / "master_samples.tsv",           sep="\t", index=False)
    per_sample_locus.to_csv( out_dir / "per_sample_locus_alleles.tsv", sep="\t", index=False)
    pairs_df.to_csv(         out_dir / "d0_dr_pairs.tsv",              sep="\t", index=False)
    pi_df.to_csv(            out_dir / "allele_freq_from_D0.tsv",      sep="\t", index=False)

    for trial_val in sorted(samples["trial"].unique()):
        subdir = out_dir / trial_val.lower()
        subdir.mkdir(parents=True, exist_ok=True)
        t = trial_val
        samples.query("trial == @t").to_csv(
            subdir / "master_samples.tsv",           sep="\t", index=False)
        per_sample_locus.query("trial == @t").to_csv(
            subdir / "per_sample_locus_alleles.tsv", sep="\t", index=False)
        pairs_df.query("trial == @t").to_csv(
            subdir / "d0_dr_pairs.tsv",              sep="\t", index=False)
        pi_df.query("trial == @t").to_csv(
            subdir / "allele_freq_from_D0.tsv",      sep="\t", index=False)

    print(f"\n✓ Wrote Pv3Rs-ready TSVs in {out_dir.resolve()}")
    print(f"  master_samples rows    : {len(samples)}")
    print(f"  per_sample_locus rows  : {len(per_sample_locus)}")
    print(f"  D0→DR pair rows        : {len(pairs_df)}")
    print(f"  allele-freq rows       : {len(pi_df)}")
    print(f"  MOIRE eMOI coverage    : "
          f"{samples['post_effective_coi_mean'].notna().sum()} / {len(samples)} samples")


if __name__ == "__main__":
    main()