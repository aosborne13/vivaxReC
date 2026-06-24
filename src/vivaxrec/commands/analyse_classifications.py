#!/usr/bin/env python3
"""
analyse_classifications.py

Reads per_infection_classifications TSV files from consecutive pairs and
produces patient-level analytical summaries. Does not re-run Pv3Rs.

All analyses use consecutive pair classifications (step-to-step Pv3Rs
posteriors). Days are taken from D0 using the to_day column, which records
the absolute day from enrolment for each recurrence sample.

Three analyses:

1. same_clone_relapses
   Consecutive steps classified as relapse (p(L) >= rel_thresh) with
   DCiFER r >= r_thresh — suggesting the same clone persisting across
   multiple recurrences rather than a new hypnozoite activation.

2. longest_relapses
   Relapse-classified recurrences sorted by to_day (days from D0).
   Reports the longest-persisting relapses and counts those exceeding
   day_thresh.

3. patient_trajectories
   Full per-patient recurrence sequence with classifications ordered
   by day from D0. Per-patient summary of relapse/reinfection/uncertain
   counts and trajectory patterns.

Usage:
  python analyse_classifications.py [flags]

Flags:
  --con_dir     DIR     Consecutive classification output directory (required)
  --out_dir     DIR     Output directory (default: same as con_dir)
  --trial       STR     Trial name for labelling
  --rel_thresh  FLOAT   p(L) threshold for relapse (default: 0.75)
  --rei_thresh  FLOAT   p(L) threshold for reinfection (default: 0.25)
  --r_thresh    FLOAT   DCiFER r threshold for same-clone definition (default: 0.95)
  --day_thresh  INT     Days from D0 threshold for long relapse flag (default: 300)
"""

import argparse
import pandas as pd
from pathlib import Path

CLASSIFICATION_FILES = {
    "full":       "per_infection_classifications.tsv",
    "first_rec":  "per_infection_classifications_first_rec_only.tsv",
    "mono_only":  "per_infection_classifications_mono_only.tsv",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--con_dir",    required=True,
                   help="Consecutive classification output directory")
    p.add_argument("--out_dir",    default=None)
    p.add_argument("--master_dir", default=None,
                   help="Directory containing <trial>/master_samples.tsv (from prep step)")
    p.add_argument("--trial",      default=None)
    p.add_argument("--rel_thresh", default=0.75, type=float)
    p.add_argument("--rei_thresh", default=0.25, type=float)
    p.add_argument("--r_thresh",   default=0.95, type=float)
    p.add_argument("--day_thresh", default=300,  type=int)
    p.add_argument("--pq_arms",    default=None, nargs="+",
                   help="PQ treatment arm name(s) e.g. --pq_arms '7-day' '14-day'")
    return p.parse_args()


def pl_call(pv3rs_L, rel_thresh, rei_thresh):
    if pd.isna(pv3rs_L):
        return "uncertain"
    if pv3rs_L >= rel_thresh:
        return "relapse"
    if pv3rs_L <= rei_thresh:
        return "reinfection"
    return "uncertain"


def load_d0_counts(master_dir, trial):
    """
    Load per-arm D0 genotyped patient counts from master_samples.tsv.
    Uses time_point == 'D0' as the denominator.
    Returns a dict of {treatment_arm: n_d0} or None if unavailable.
    """
    if not master_dir or not trial:
        return None
    fp = Path(master_dir) / trial.lower() / "master_samples.tsv"
    if not fp.exists():
        print(f"   WARNING: master_samples.tsv not found at {fp} — no-recurrence counts omitted")
        return None
    ms = pd.read_csv(fp, sep="\t")
    if not {"time_point", "treatment_arm", "patient_id"}.issubset(ms.columns):
        print("   WARNING: master_samples.tsv missing required columns — no-recurrence counts omitted")
        return None
    d0 = ms[ms["time_point"] == "D0"]
    counts = (d0.groupby("treatment_arm", dropna=False)["patient_id"]
                .nunique()
                .to_dict())
    # Normalise keys to lowercase for case-insensitive matching downstream
    return {str(k): v for k, v in counts.items()}



    if pd.isna(pv3rs_L):
        return "uncertain"
    if pv3rs_L >= rel_thresh:
        return "relapse"
    if pv3rs_L <= rei_thresh:
        return "reinfection"
    return "uncertain"


def load_file(fp):
    df = pd.read_csv(fp, sep="\t")
    for col in ["pv3rs_L", "pv3rs_C", "pv3rs_I", "dcifer_r",
                "emoi_s1", "emoi_s2", "to_day", "from_day",
                "step_no", "shared_loci", "loci_tested", "share_frac"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


# =============================================================================
# Analysis 1: Same-clone relapses
# =============================================================================
def analyse_same_clone_relapses(df, args, table_type):
    """
    Identify consecutive steps classified as relapse where DCiFER r >= r_thresh.
    Days are from D0 via to_day.
    """
    df = df.copy()
    df["pl_call"]       = df["pv3rs_L"].apply(lambda x: pl_call(x, args.rel_thresh, args.rei_thresh))
    df["is_relapse"]    = df["pl_call"] == "relapse"
    df["is_same_clone"] = df["is_relapse"] & (df["dcifer_r"] >= args.r_thresh)
    df["day"]           = df["to_day"]

    # Patient-level
    pat = df.groupby(["id", "treatment_arm", "site"], dropna=False).agg(
        n_steps            = ("step_no",       "count"),
        n_relapse_steps    = ("is_relapse",    "sum"),
        n_same_clone_steps = ("is_same_clone", "sum"),
        max_dcifer_r       = ("dcifer_r",      "max"),
        last_day           = ("day",           "max"),
    ).reset_index()

    pat["has_same_clone_relapse"] = pat["n_same_clone_steps"] >= 1
    pat["has_multi_same_clone"]   = pat["n_same_clone_steps"] >= 2
    pat["trial"]        = args.trial or ""
    pat["table_type"]   = table_type
    pat["r_threshold"]  = args.r_thresh
    pat["rel_threshold"] = args.rel_thresh

    # Arm-level
    arm = pat.groupby("treatment_arm", dropna=False).agg(
        n_patients           = ("id",                     "count"),
        n_with_relapse       = ("n_relapse_steps",        lambda x: (x >= 1).sum()),
        n_same_clone         = ("has_same_clone_relapse", "sum"),
        n_multi_same_clone   = ("has_multi_same_clone",   "sum"),
        mean_same_clone_steps = ("n_same_clone_steps",    "mean"),
        max_same_clone_steps  = ("n_same_clone_steps",    "max"),
    ).reset_index()
    arm["pct_same_clone"]  = round(100 * arm["n_same_clone"] / arm["n_patients"].clip(lower=1), 1)
    arm["trial"]           = args.trial or ""
    arm["table_type"]      = table_type
    arm["r_threshold"]     = args.r_thresh

    return pat, arm


# =============================================================================
# Analysis 2: Longest relapses
# =============================================================================
def analyse_longest_relapses(df, args, table_type):
    """
    Find relapse-classified recurrences with the highest to_day value.
    to_day is the absolute day from D0 for that recurrence.
    """
    df = df.copy()
    df["pl_call"] = df["pv3rs_L"].apply(lambda x: pl_call(x, args.rel_thresh, args.rei_thresh))
    df["day"]     = df["to_day"]

    relapses = df[df["pl_call"] == "relapse"].copy()
    if not len(relapses):
        print("   longest_relapses: no relapses found")
        return None, None

    relapses = relapses.sort_values("day", ascending=False)
    relapses["exceeds_day_thresh"] = relapses["day"] > args.day_thresh
    relapses["trial"]       = args.trial or ""
    relapses["table_type"]  = table_type
    relapses["day_threshold"] = args.day_thresh

    arm = relapses.groupby("treatment_arm", dropna=False).agg(
        n_relapses      = ("id",                 "count"),
        n_patients      = ("id",                 "nunique"),
        n_exceed_thresh = ("exceeds_day_thresh", "sum"),
        max_day         = ("day",                "max"),
        median_day      = ("day",                "median"),
        mean_day        = ("day",                "mean"),
    ).reset_index()
    arm["pct_exceed_thresh"] = round(
        100 * arm["n_exceed_thresh"] / arm["n_relapses"].clip(lower=1), 1)
    arm["trial"]        = args.trial or ""
    arm["table_type"]   = table_type
    arm["day_threshold"] = args.day_thresh

    return relapses, arm


# =============================================================================
# Analysis 3: Patient trajectories
# =============================================================================
def build_patient_trajectories(df, args, table_type):
    """
    Full per-patient recurrence sequence ordered by day from D0 (to_day).
    """
    df = df.copy()
    df["pl_call"] = df["pv3rs_L"].apply(lambda x: pl_call(x, args.rel_thresh, args.rei_thresh))
    df["day"]     = df["to_day"]
    df["trial"]   = args.trial or ""
    df["table_type"] = table_type

    keep = ["trial", "table_type", "id", "treatment_arm", "site",
            "step_no", "day", "from_day", "to_day",
            "pv3rs_L", "pv3rs_I", "pv3rs_C", "pv3rs_top",
            "dcifer_r", "emoi_s1", "emoi_s2", "share_frac", "pl_call"]
    keep = [c for c in keep if c in df.columns]
    traj = df[keep].sort_values(["id", "day"], na_position="last")

    # Per-patient summary
    pat = traj.groupby(["id", "treatment_arm", "site", "trial"], dropna=False).agg(
        n_recurrences  = ("day",      "count"),
        n_relapse      = ("pl_call",  lambda x: (x == "relapse").sum()),
        n_reinfection  = ("pl_call",  lambda x: (x == "reinfection").sum()),
        n_uncertain    = ("pl_call",  lambda x: (x == "uncertain").sum()),
        first_rec_day  = ("day",      "min"),
        last_rec_day   = ("day",      "max"),
        max_pv3rs_L    = ("pv3rs_L",  "max"),
        max_dcifer_r   = ("dcifer_r", "max"),
    ).reset_index()

    pat["all_relapse"]      = pat["n_relapse"] == pat["n_recurrences"]
    pat["mixed_trajectory"] = (pat["n_relapse"] > 0) & (pat["n_reinfection"] > 0)
    pat["table_type"]       = table_type

    return traj, pat


# =============================================================================
# Narrative summary
# =============================================================================
def write_narrative_summary(pat_traj, pat_sc, pat_lr, args, table_type, out_dir,
                            d0_counts=None, pq_arms=None):
    """
    Write a plain-English narrative summary of key findings per treatment arm,
    suitable for use as draft manuscript text.
    """
    trial = args.trial or "Unknown trial"
    lines = []

    def h(text):
        lines.append(f"\n{'='*70}")
        lines.append(text)
        lines.append('='*70)

    def p(text=""):
        lines.append(text)

    h(f"NARRATIVE SUMMARY — {trial}  |  table: {table_type}")
    p(f"Relapse threshold: p(L) >= {args.rel_thresh}  |  "
      f"Reinfection threshold: p(L) <= {args.rei_thresh}")
    p(f"Same-clone DCiFER r threshold: >= {args.r_thresh}  |  "
      f"Long relapse day threshold: > {args.day_thresh} days")

    arms = sorted(pat_traj["treatment_arm"].dropna().unique())

    for arm in arms:
        h(f"ARM: {arm}")

        # ── Trajectory stats ──────────────────────────────────────────────────
        tp = pat_traj[pat_traj["treatment_arm"] == arm]
        n_total       = len(tp)
        n_any_rel     = (tp["n_relapse"] >= 1).sum()
        n_multi_rel   = (tp["n_relapse"] >= 2).sum()
        n_all_rel     = tp["all_relapse"].sum()
        n_mixed       = tp["mixed_trajectory"].sum()
        n_any_rei     = (tp["n_reinfection"] >= 1).sum()
        max_rel       = int(tp["n_relapse"].max())
        max_rel_id    = tp.loc[tp["n_relapse"].idxmax(), "id"] if n_total > 0 else "N/A"
        max_day       = tp["last_rec_day"].max()

        # Pair-level counts (match table denominators)
        n_pairs_total = int(tp["n_recurrences"].sum())
        n_pairs_rel   = int(tp["n_relapse"].sum())
        n_pairs_rei   = int(tp["n_reinfection"].sum())
        n_pairs_unc   = int(tp["n_uncertain"].sum()) if "n_uncertain" in tp.columns else None

        p(f"\nPATIENT-LEVEL (one count per individual):")
        n_d0 = 0
        n_no_rec = 0
        if d0_counts:
            # Case-insensitive lookup
            arm_key = next((k for k in d0_counts if k.lower() == str(arm).lower()), None)
            n_d0    = d0_counts.get(arm_key, 0) if arm_key else 0
            n_no_rec = max(0, n_d0 - n_total)
            p(f"  D0 genotyped patients:                 {n_d0}")
            p(f"  No recurrence (or no classified pair): {n_no_rec} "
              f"({100*n_no_rec/max(n_d0,1):.0f}%)")
        p(f"  Patients with >= 1 classified consecutive pair: {n_total}"
          + (f" ({100*n_total/max(n_d0,1):.0f}% of D0 genotyped)"
             if d0_counts and n_d0 > 0 else ""))
        p(f"  - Any relapse step (p(L) >= {args.rel_thresh}):          {n_any_rel} "
          f"({100*n_any_rel/max(n_total,1):.0f}%)")
        p(f"  - >= 2 relapse steps:                  {n_multi_rel} "
          f"({100*n_multi_rel/max(n_total,1):.0f}%)")
        p(f"  - All steps relapse (pure relapsers):  {n_all_rel} "
          f"({100*n_all_rel/max(n_total,1):.0f}%)")
        p(f"  - Any reinfection step (p(L) <= {args.rei_thresh}):        {n_any_rei} "
          f"({100*n_any_rei/max(n_total,1):.0f}%)")
        p(f"  - Mixed (relapse AND reinfection):     {n_mixed} "
          f"({100*n_mixed/max(n_total,1):.0f}%)")

        p(f"\nPAIR-LEVEL (one count per classified pair, matches summary table):")
        p(f"  Total classified pairs:   {n_pairs_total}")
        p(f"  - Relapse pairs:          {n_pairs_rel} "
          f"({100*n_pairs_rel/max(n_pairs_total,1):.1f}%)")
        p(f"  - Reinfection pairs:      {n_pairs_rei} "
          f"({100*n_pairs_rei/max(n_pairs_total,1):.1f}%)")
        if n_pairs_unc is not None:
            p(f"  - Indeterminate pairs:    {n_pairs_unc} "
              f"({100*n_pairs_unc/max(n_pairs_total,1):.1f}%)")

        if max_rel > 0:
            p(f"\n  Maximum relapse steps in one patient: {max_rel} (patient {max_rel_id})")
        if pd.notna(max_day):
            p(f"  Latest observed recurrence: day {int(max_day)} from D0")

        # Draft sentence
        p(f"\n  DRAFT: In the {arm} arm, {n_any_rel} of {n_total} patients "
          f"({100*n_any_rel/max(n_total,1):.0f}%) experienced at least one relapse. "
          f"{n_multi_rel} patients had two or more consecutive relapse steps")
        if d0_counts and n_d0 > 0:
            p(f"\nD0-DENOMINATOR PERSPECTIVE (all genotyped patients at risk):")
            p(f"  This uses all {n_d0} D0-genotyped patients as the denominator,")
            p(f"  including {n_no_rec} who had no classified recurrence.")
            p(f"  - Patients with >= 1 relapse:     {n_any_rel} / {n_d0} "
              f"({100*n_any_rel/n_d0:.1f}%)")
            p(f"  - Patients with >= 1 reinfection: {n_any_rei} / {n_d0} "
              f"({100*n_any_rei/n_d0:.1f}%)")
            p(f"  - Patients with no recurrence:    {n_no_rec} / {n_d0} "
              f"({100*n_no_rec/n_d0:.1f}%)")
            p(f"  - Total relapse pairs (all steps): {n_pairs_rel} across {n_any_rel} patients")
            p()
            p(f"  DRAFT (relapses avoided framing): Of {n_d0} patients genotyped")
            p(f"  at D0 in the {arm} arm, {n_any_rel} ({100*n_any_rel/n_d0:.1f}%) experienced")
            p(f"  at least one relapse, contributing {n_pairs_rel} relapse-classified")
            p(f"  recurrences in total. {n_no_rec} patients ({100*n_no_rec/n_d0:.1f}%)")
            p(f"  had no detected recurrence during follow-up.")

        # ── Same-clone relapses ───────────────────────────────────────────────
        if pat_sc is not None:
            sc = pat_sc[pat_sc["treatment_arm"] == arm]
            n_sc       = sc["has_same_clone_relapse"].sum()
            n_sc_multi = sc["has_multi_same_clone"].sum()
            max_sc     = int(sc["n_same_clone_steps"].max()) if len(sc) else 0

            p(f"\nSame-clone relapses (consecutive relapse steps, DCiFER r >= {args.r_thresh}):")
            p(f"  - Patients with >= 1 same-clone step: {n_sc} "
              f"({100*n_sc/max(n_total,1):.0f}%)")
            p(f"  - Patients with >= 2 same-clone steps: {n_sc_multi} "
              f"({100*n_sc_multi/max(n_total,1):.0f}%)")
            p(f"  - Maximum same-clone steps: {max_sc}")

            if n_sc > 0:
                p(f"\n  DRAFT: {n_sc} patient(s) in the {arm} arm had consecutive relapses "
                  f"with DCiFER r >= {args.r_thresh}, consistent with persistence of the same "
                  f"parasite clone across multiple recurrences.")

        # ── Longest relapses ──────────────────────────────────────────────────
        if pat_lr is not None:
            lr = pat_lr[pat_lr["treatment_arm"] == arm]
            n_exceed  = int((lr["day"] > args.day_thresh).sum())
            max_day_r = lr["day"].max()

            p(f"\nLongest relapses (days from D0):")
            p(f"  - Total relapse-classified recurrences: {len(lr)}")
            p(f"  - Exceeding day {args.day_thresh}: {n_exceed}")
            if pd.notna(max_day_r):
                p(f"  - Latest relapse: day {int(max_day_r)}")

            if n_exceed > 0:
                p(f"\n  DRAFT: {n_exceed} relapse(s) in the {arm} arm occurred beyond "
                  f"day {args.day_thresh} from enrolment, the latest at day "
                  f"{int(max_day_r) if pd.notna(max_day_r) else 'unknown'}.")

    # ── Overall summary ───────────────────────────────────────────────────────
    h("OVERALL SUMMARY")

    pq_arms_lower = [a.lower() for a in (pq_arms or [])]

    def arm_group_summary(label, subset_arms):
        tp_g = pat_traj[pat_traj["treatment_arm"].isin(subset_arms)]
        if not len(tp_g):
            return
        n_tot_g       = len(tp_g)
        n_rel_g       = (tp_g["n_relapse"] >= 1).sum()
        n_rei_g       = (tp_g["n_reinfection"] >= 1).sum()
        n_mix_g       = tp_g["mixed_trajectory"].sum()
        n_all_rel_g   = tp_g["all_relapse"].sum()
        n_pairs_rel_g = int(tp_g["n_relapse"].sum())
        n_pairs_rei_g = int(tp_g["n_reinfection"].sum())
        n_pairs_unc_g = int(tp_g["n_uncertain"].sum()) if "n_uncertain" in tp_g.columns else 0
        n_pairs_g     = int(tp_g["n_recurrences"].sum())

        p(f"\n{label} ({', '.join(sorted(subset_arms))}):")
        n_d0_g = 0
        n_no_g = 0
        if d0_counts:
            n_d0_g = sum(d0_counts.get(
                next((k for k in d0_counts if k.lower() == a.lower()), a), 0)
                for a in subset_arms)
            n_no_g = max(0, n_d0_g - n_tot_g)

        denom = n_d0_g if n_d0_g > 0 else n_tot_g
        denom_label = f"of {n_d0_g} D0 genotyped" if (d0_counts and n_d0_g > 0) \
                      else f"of {n_tot_g} with classified pair"

        if d0_counts and n_d0_g > 0:
            p(f"  D0 genotyped:                   {n_d0_g}")
            p(f"  No recurrence / unclassified:   {n_no_g} "
              f"({100*n_no_g/max(denom,1):.1f}%)")
            p(f"  With classified pair:           {n_tot_g} "
              f"({100*n_tot_g/max(denom,1):.1f}%)")
        else:
            p(f"  Patients with classified pair:  {n_tot_g}")
        p(f"  Any relapse [{denom_label}]:  {n_rel_g} "
          f"({100*n_rel_g/max(denom,1):.1f}%)")
        p(f"  Any reinfection [{denom_label}]: {n_rei_g} "
          f"({100*n_rei_g/max(denom,1):.1f}%)")
        p(f"  Mixed trajectory [{denom_label}]: {n_mix_g} "
          f"({100*n_mix_g/max(denom,1):.1f}%)")
        p(f"  Pure relapse, all steps [{denom_label}]: {n_all_rel_g} "
          f"({100*n_all_rel_g/max(denom,1):.1f}%)")
        p(f"  PAIR-LEVEL — total pairs:       {n_pairs_g}")
        p(f"    Relapse:                      {n_pairs_rel_g} "
          f"({100*n_pairs_rel_g/max(n_pairs_g,1):.1f}%)")
        p(f"    Reinfection:                  {n_pairs_rei_g} "
          f"({100*n_pairs_rei_g/max(n_pairs_g,1):.1f}%)")
        p(f"    Indeterminate:                {n_pairs_unc_g} "
          f"({100*n_pairs_unc_g/max(n_pairs_g,1):.1f}%)")
        if pat_sc is not None:
            sc_g   = pat_sc[pat_sc["treatment_arm"].isin(subset_arms)]
            n_sc_g = sc_g["has_same_clone_relapse"].sum() if len(sc_g) else 0
            n_sc_pairs_g = int(sc_g["n_same_clone_steps"].sum()) if len(sc_g) else 0
            p(f"  Same-clone relapses (DCiFER r >= {args.r_thresh}):")
            p(f"    Patients with >= 1 same-clone step: {n_sc_g} / {n_rel_g} relapsers "
              f"({100*n_sc_g/max(n_rel_g,1):.1f}% of those with any relapse)")
            p(f"    Same-clone relapse pairs:           {n_sc_pairs_g} / {n_pairs_rel_g} "
              f"({100*n_sc_pairs_g/max(n_pairs_rel_g,1):.1f}% of all relapse pairs)")
        if pat_lr is not None:
            lr_g  = pat_lr[pat_lr["treatment_arm"].isin(subset_arms)]
            n_ex_g = int((lr_g["day"] > args.day_thresh).sum()) if len(lr_g) else 0
            p(f"  Relapses > day {args.day_thresh}:          {n_ex_g}")

    all_arms     = sorted(pat_traj["treatment_arm"].dropna().unique())
    pq_arm_set   = {a for a in all_arms if a.lower() in pq_arms_lower}
    nonpq_arm_set = {a for a in all_arms if a.lower() not in pq_arms_lower}

    if pq_arms and pq_arm_set:
        arm_group_summary("NON-PQ ARMS", nonpq_arm_set)
        arm_group_summary("PQ ARMS", pq_arm_set)
    else:
        arm_group_summary("ALL ARMS", set(all_arms))

    if d0_counts:
        p(f"\nTOTAL D0 genotyped across all arms: {sum(d0_counts.values())}")

    p()
    return "\n".join(lines)

    if pat_sc is not None:
        n_sc_tot = pat_sc["has_same_clone_relapse"].sum()
        p(f"  Same-clone relapses:  {n_sc_tot} ({100*n_sc_tot/max(n_tot,1):.0f}%)")

    if pat_lr is not None:
        n_exceed_tot = int((pat_lr["day"] > args.day_thresh).sum())
        p(f"  Relapses > day {args.day_thresh}:    {n_exceed_tot}")

    p()
    return "\n".join(lines)


# =============================================================================
# Main
# =============================================================================
def main():
    args = parse_args()
    con_dir = Path(args.con_dir)
    out_dir = Path(args.out_dir) if args.out_dir else con_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n>> analyse_classifications")
    print(f"   trial:       {args.trial or '(not specified)'}")
    print(f"   con_dir:     {con_dir}")
    print(f"   rel_thresh:  p(L) >= {args.rel_thresh}")
    print(f"   rei_thresh:  p(L) <= {args.rei_thresh}")
    print(f"   r_thresh:    DCiFER r >= {args.r_thresh} (same-clone definition)")
    print(f"   day_thresh:  > {args.day_thresh} days from D0")
    print(f"   pq_arms:     {args.pq_arms or '(none specified — no PQ/non-PQ split)'}")

    d0_counts = load_d0_counts(args.master_dir, args.trial)
    if d0_counts:
        print(f"   D0 counts loaded: {sum(d0_counts.values())} patients across "
              f"{len(d0_counts)} arms")
        print(f"   D0 arms: {list(d0_counts.keys())}")

    for table_type, fname in CLASSIFICATION_FILES.items():
        fp = con_dir / fname
        if not fp.exists():
            print(f"\n   {fname} not found — skipping {table_type}")
            continue

        print(f"\n>> reading {fname} ({table_type}) ...")
        df = load_file(fp)
        print(f"   {len(df)} rows, {df['id'].nunique()} patients")

        if "pv3rs_L" not in df.columns:
            print("   pv3rs_L not found — skipping (Pv3Rs may not have run)")
            continue

        suffix = f"_{table_type}"

        # ── 1. Same-clone relapses ────────────────────────────────────────────
        pat_sc, arm_sc = analyse_same_clone_relapses(df, args, table_type)
        pat_sc.to_csv(out_dir / f"same_clone_relapses_patient{suffix}.tsv",  sep="\t", index=False)
        arm_sc.to_csv(out_dir / f"same_clone_relapses_summary{suffix}.tsv", sep="\t", index=False)
        n_sc = pat_sc["has_same_clone_relapse"].sum()
        print(f"   same-clone relapses: {n_sc} patients with >= 1 same-clone step "
              f"(r >= {args.r_thresh})")

        # ── 2. Longest relapses ───────────────────────────────────────────────
        pat_lr, arm_lr = analyse_longest_relapses(df, args, table_type)
        if pat_lr is not None:
            pat_lr.to_csv(out_dir / f"longest_relapses_patient{suffix}.tsv",  sep="\t", index=False)
            arm_lr.to_csv(out_dir / f"longest_relapses_summary{suffix}.tsv", sep="\t", index=False)
            print(f"   longest relapse: {int(pat_lr['day'].max())} days from D0; "
                  f"{int((pat_lr['day'] > args.day_thresh).sum())} exceed {args.day_thresh} days")

        # ── 3. Patient trajectories ───────────────────────────────────────────
        traj, pat_traj = build_patient_trajectories(df, args, table_type)
        traj.to_csv(    out_dir / f"patient_trajectories{suffix}.tsv",        sep="\t", index=False)
        pat_traj.to_csv(out_dir / f"patient_trajectory_summary{suffix}.tsv", sep="\t", index=False)
        print(f"   trajectories: {len(pat_traj)} patients; "
              f"{pat_traj['all_relapse'].sum()} all-relapse, "
              f"{pat_traj['mixed_trajectory'].sum()} mixed")

        # ── 4. Narrative summary ──────────────────────────────────────────────
        narrative = write_narrative_summary(
            pat_traj, pat_sc, pat_lr, args, table_type, out_dir,
            d0_counts=d0_counts, pq_arms=args.pq_arms)
        fp_narr = out_dir / f"narrative_summary{suffix}.txt"
        fp_narr.write_text(narrative)
        print(f"   wrote narrative: {fp_narr.name}")

    print("\n✓ done")

if __name__ == '__main__':
    main()
