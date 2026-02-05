from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd


def safe_read(p: str | Path) -> pd.DataFrame:
    p = Path(p)
    if not p.exists() or p.stat().st_size == 0:
        raise FileNotFoundError(f"File not found or empty: {p}")
    try:
        return pd.read_csv(p, sep="\t")
    except Exception:
        return pd.read_csv(p)


def parse_external_from_sample_id(sample_id: str, split_char: str = "_") -> str:
    if not isinstance(sample_id, str):
        return ""
    return sample_id.split(split_char, 1)[0]


def parse_tp_index(tp: str) -> int:
    """Map time_point strings: D0->0, DRP1->1, DRP2->2, ... else parse trailing number if present."""
    if not isinstance(tp, str):
        return 10**9
    s = tp.strip().upper()
    if s == "D0":
        return 0
    m = re.search(r"DRP(\d+)", s)
    if m:
        return int(m.group(1))
    m = re.search(r"(\d+)$", s)
    return int(m.group(1)) if m else 10**9


def collapse_to_sets(df: pd.DataFrame) -> pd.DataFrame:
    """From long rows (sample_id, locus, allele) -> per sample×locus unique allele set ('|' joined)."""

    def uniq_join(series: pd.Series) -> str:
        vals = [str(x) for x in series.dropna().tolist()]
        vals = [v for v in vals if v and v.lower() != "nan"]
        return "|".join(sorted(set(vals)))

    return (
        df.groupby(["sample_id", "locus"], as_index=False)
        .agg(alleles_set=("allele", uniq_join))
    )


def alleles_count(alleles_set_str: str) -> int:
    if (
        not isinstance(alleles_set_str, str)
        or alleles_set_str == ""
        or alleles_set_str.lower() == "nan"
    ):
        return 0
    return len([a for a in alleles_set_str.split("|") if a])


def _normalize_two_cols(cols: list[str]) -> Tuple[str, str]:
    """
    Accept either:
      --coi-cols sample_id coi
    or (less recommended):
      --coi-cols sample_id,coi
    """
    if len(cols) == 1 and "," in cols[0]:
        parts = [x.strip() for x in cols[0].split(",") if x.strip()]
        if len(parts) != 2:
            raise ValueError(f"--coi-cols must contain exactly two columns, got: {cols}")
        return parts[0], parts[1]
    if len(cols) != 2:
        raise ValueError(f"--coi-cols must contain exactly two columns, got: {cols}")
    return cols[0], cols[1]


def run_prepare_input_data_single_dataset(
    *,
    dataset_name: str,
    genotypes_file: str | Path,
    metadata_file: str | Path,
    outdir: str | Path,
    # Optional COI source
    coi_source_file: Optional[str | Path] = None,
    coi_cols: Optional[list[str]] = None,  # [sample_col, coi_col]
    # Metadata controls
    metadata_dataset_col: Optional[str] = None,  # e.g. "trial" or "dataset" if present
    sample_id_col: str = "sample_id",
    locus_col: str = "locus",
    allele_col: str = "allele",
    day_col: str = "day",
    time_point_col: str = "time_point",
    patient_id_col: str = "patient_id",
    split_char: str = "_",
) -> None:
    """
    Prepare Aster-ready inputs for ONE dataset.

    Required genotype columns (defaults):
      sample_id, locus, allele

    Required metadata columns (defaults):
      sample_id, day, time_point
    Optional metadata:
      patient_id (if absent, derived from sample_id split on split_char)

    Optional COI:
      provide coi_source_file + coi_cols: [sample_col, coi_col]
      If absent, final_coi falls back to a naive lower bound from allele counts per locus.
    """
    dataset_name = str(dataset_name)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- load genotypes ----
    g_raw = safe_read(genotypes_file)
    for c in [sample_id_col, locus_col, allele_col]:
        if c not in g_raw.columns:
            raise RuntimeError(f"Genotypes missing column '{c}'. Present: {list(g_raw.columns)}")

    g_raw = g_raw.rename(
        columns={sample_id_col: "sample_id", locus_col: "locus", allele_col: "allele"}
    )
    geno = collapse_to_sets(g_raw).assign(dataset=dataset_name)  # sample_id,locus,alleles_set,dataset

    # ---- load metadata ----
    meta = safe_read(metadata_file)

    # If the metadata has a dataset column, optionally filter
    if metadata_dataset_col:
        if metadata_dataset_col not in meta.columns:
            raise RuntimeError(
                f"--metadata-dataset-col '{metadata_dataset_col}' not found in metadata. "
                f"Columns: {list(meta.columns)}"
            )
        meta = meta.loc[meta[metadata_dataset_col].astype(str).str.lower() == dataset_name.lower()].copy()

    # Ensure required metadata columns exist
    for c in [sample_id_col, day_col, time_point_col]:
        if c not in meta.columns:
            raise RuntimeError(f"Metadata missing column '{c}'. Present: {list(meta.columns)}")

    # Standardize names used downstream
    meta = meta.rename(
        columns={
            sample_id_col: "sample_id",
            day_col: "day",
            time_point_col: "time_point",
        }
    )
    meta["day"] = pd.to_numeric(meta["day"], errors="coerce")

    # patient_id: use provided column if present, else derive from sample_id
    if patient_id_col in meta.columns:
        meta = meta.rename(columns={patient_id_col: "patient_id"})
    else:
        meta["patient_id"] = meta["sample_id"].astype(str).apply(
            lambda s: parse_external_from_sample_id(s, split_char=split_char)
        )

    # ---- COI table (optional) ----
    coi_tbl = pd.DataFrame(columns=["sample_id", "coi_estimate", "coi_source"])
    if coi_source_file is not None:
        if not coi_cols:
            raise RuntimeError("If --coi-source is provided, you must also provide --coi-cols.")
        coi_sample_col, coi_value_col = _normalize_two_cols(coi_cols)

        coi_in = safe_read(coi_source_file)
        for c in [coi_sample_col, coi_value_col]:
            if c not in coi_in.columns:
                raise RuntimeError(
                    f"COI file missing column '{c}'. Present: {list(coi_in.columns)}"
                )

        coi_tbl = (
            coi_in[[coi_sample_col, coi_value_col]]
            .rename(columns={coi_sample_col: "sample_id", coi_value_col: "coi_estimate"})
            .copy()
        )
        coi_tbl["coi_estimate"] = pd.to_numeric(coi_tbl["coi_estimate"], errors="coerce")
        coi_tbl["coi_source"] = Path(coi_source_file).name

    # ---- master samples ----
    samples = (
        geno[["dataset", "sample_id"]]
        .drop_duplicates()
        .merge(meta, on="sample_id", how="left")
    )

    # naive COI lower bound from genotype data: max allele count per locus
    coi_lb = (
        geno.groupby(["dataset", "sample_id"])["alleles_set"]
        .apply(lambda s: max((alleles_count(x) for x in s), default=0))
        .reset_index()
        .rename(columns={"alleles_set": "coi_naive_lb"})
    )

    samples = samples.merge(coi_lb, on=["dataset", "sample_id"], how="left").merge(
        coi_tbl, on="sample_id", how="left"
    )

    samples["final_coi"] = samples["coi_estimate"].fillna(samples["coi_naive_lb"])
    samples["coi_source"] = samples["coi_source"].fillna("naive_lower_bound")

    # ---- per-sample × locus with metadata ----
    per_sample_locus = geno.merge(
        samples[["dataset", "sample_id", "patient_id", "day", "time_point", "final_coi", "coi_source"]],
        on=["dataset", "sample_id"],
        how="left",
    )

    # ---- D0 → DR pairs per patient ----
    pairs = []
    for patient, sub in per_sample_locus.groupby(["patient_id"], dropna=False):
        # order by day then time_point index then sample_id
        sids = (
            sub[["sample_id", "day", "time_point"]]
            .drop_duplicates(subset=["sample_id"])
            .copy()
        )
        sids["tp_idx"] = sids["time_point"].apply(parse_tp_index)
        sids = sids.sort_values(["day", "tp_idx", "sample_id"], na_position="last")
        if len(sids) == 0:
            continue
        d0_sid = sids.iloc[0]["sample_id"]
        later = sids.loc[sids["sample_id"] != d0_sid].copy()
        for _, r in later.iterrows():
            pairs.append(
                {
                    "dataset": dataset_name,
                    "patient_id": patient,
                    "d0_sample_id": d0_sid,
                    "dr_sample_id": r["sample_id"],
                    "dr_day": r["day"],
                    "dr_time_point": r["time_point"],
                }
            )
    pairs_df = pd.DataFrame(pairs)

    # ---- allele frequencies (π) from D0 only ----
    # Identify D0 sample per patient
    d0_rows = []
    for patient, sub in per_sample_locus.groupby(["patient_id"], dropna=False):
        sids = (
            sub[["sample_id", "day", "time_point"]]
            .drop_duplicates(subset=["sample_id"])
            .copy()
        )
        sids["tp_idx"] = sids["time_point"].apply(parse_tp_index)
        sids = sids.sort_values(["day", "tp_idx", "sample_id"], na_position="last")
        if len(sids) == 0:
            continue
        d0_rows.append({"sample_id": sids.iloc[0]["sample_id"]})

    d0_df = pd.DataFrame(d0_rows)
    d0_sub = per_sample_locus.merge(d0_df, on="sample_id", how="inner")

    pi_rows = []
    for locus, sub in d0_sub.groupby(["locus"]):
        n = sub["sample_id"].nunique()
        counts: Counter[str] = Counter()
        for s in sub["alleles_set"].astype(str):
            alleles = [a for a in s.split("|") if a]
            counts.update(set(alleles))  # presence per sample
        for a, cnt in counts.items():
            pi_rows.append(
                {
                    "dataset": dataset_name,
                    "locus": locus,
                    "allele": a,
                    "freq": cnt / n if n else 0.0,
                    "n_samples": n,
                }
            )

    pi_df = pd.DataFrame(pi_rows).sort_values(["dataset", "locus", "allele"])

    # ---- write outputs ----
    samples.to_csv(outdir / "master_samples.tsv", sep="\t", index=False)
    per_sample_locus.to_csv(outdir / "per_sample_locus_alleles.tsv", sep="\t", index=False)
    pairs_df.to_csv(outdir / "d0_dr_pairs.tsv", sep="\t", index=False)
    pi_df.to_csv(outdir / "allele_freq_from_D0.tsv", sep="\t", index=False)

    print(f"✓ Wrote Aster-ready TSVs in {outdir.resolve()}")
