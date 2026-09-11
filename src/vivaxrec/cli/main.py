#!/usr/bin/env python3
"""
vivaxrec — Plasmodium vivax recurrence classification pipeline.

Steps:
  1. standardize   Clean raw input files and split into analysis-ready, dropped, and QC outputs
  2. build         Build patient-level infection sequences and recurrence pair tables
  3. prep          Format data for Pv3Rs (produces master_samples.tsv etc.)
  4. classify      Run recurrence classification (DCiFER + Pv3Rs)
  5. summarise     Summarise classification outputs without re-running Pv3Rs
  6. plot          Generate classification figures
  7. run-all       Run the full pipeline end to end (steps 1–6)

Usage:
  vivaxrec standardize [options]
  vivaxrec build       [options]
  vivaxrec prep        [options]
  vivaxrec classify    [options]
  vivaxrec plot        [options]
  vivaxrec run-all     --trials TRIAL [TRIAL ...] [options]
"""

import argparse
import sys
import subprocess
from pathlib import Path

COMMANDS_DIR = Path(__file__).parent.parent / "commands"

R_CMD  = "Rscript"
PY_CMD = sys.executable


def _run(cmd):
    print(f"[vivaxrec] Running: {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(result.returncode)


def cmd_standardize(extra):
    _run([PY_CMD, str(COMMANDS_DIR / "standardize_inputs.py")] + extra)


def cmd_build(extra):
    _run([PY_CMD, str(COMMANDS_DIR / "build_patient_level_tables.py")] + extra)


def cmd_prep(extra):
    _run([PY_CMD, str(COMMANDS_DIR / "prep_pv3rs_inputs.py")] + extra)


def cmd_analyse(extra):
    _run([PY_CMD, str(COMMANDS_DIR / "analyse_classifications.py")] + extra)


def cmd_summarise(extra):
    _run([PY_CMD, str(COMMANDS_DIR / "summarise_classifications.py")] + extra)


def cmd_classify(extra):
    _run([R_CMD, str(COMMANDS_DIR / "classify_recurrences.R")] + extra)


def cmd_plot(extra):
    _run([R_CMD, str(COMMANDS_DIR / "plot_recurrence_classifications.R")] + extra)


def cmd_run_all(extra):
    # Extract --trials; remaining flags forwarded to each step
    trials = []
    remaining = []
    i = 0
    while i < len(extra):
        if extra[i] == "--trials":
            i += 1
            while i < len(extra) and not extra[i].startswith("--"):
                trials.append(extra[i])
                i += 1
        else:
            remaining.append(extra[i])
            i += 1

    if not trials:
        print("[vivaxrec] ERROR: --trials is required for run-all.")
        print("  Example: vivaxrec run-all --trials mytrial --meta /data/metadata.tsv")
        sys.exit(1)

    print(f"\n[vivaxrec] === Step 1/5: standardize ===")
    _run([PY_CMD, str(COMMANDS_DIR / "standardize_inputs.py")] + remaining)

    print(f"\n[vivaxrec] === Step 2/5: build ===")
    _run([PY_CMD, str(COMMANDS_DIR / "build_patient_level_tables.py")] + remaining)

    print(f"\n[vivaxrec] === Step 3/5: prep ===")
    _run([PY_CMD, str(COMMANDS_DIR / "prep_pv3rs_inputs.py")] + remaining)

    for trial in trials:
        for pair_mode in ["baseline", "consecutive"]:
            print(f"\n[vivaxrec] === Step 4/5: classify — trial={trial} pair_mode={pair_mode} ===")
            _run([R_CMD, str(COMMANDS_DIR / "classify_recurrences.R"),
                  "trial", trial, "pair_mode", pair_mode] + remaining)

    for trial in trials:
        print(f"\n[vivaxrec] === Step 5/5: plot — trial={trial} ===")
        _run([R_CMD, str(COMMANDS_DIR / "plot_recurrence_classifications.R"),
              "trial", trial] + remaining)

    print(f"\n[vivaxrec] ✓ Pipeline complete for trial(s): {', '.join(trials)}")


# Subcommand registry: name -> (function, one-line help, full description)
SUBCOMMANDS = {
    "standardize": (
        cmd_standardize,
        "Clean raw input files; produce analysis-ready, dropped, and QC outputs",
        """\
Clean raw input files and filter to the metadata universe.

Options:
  --rel      PATH   Pairwise relatedness TSV (mhap_between_relatedness_estimate.tsv)
  --coi      PATH   MOIRE COI summary TSV (mhap_effective_COI_summary.tsv)
  --meta     PATH   Sample metadata TSV
  --out_dir  DIR    Root output directory (default: standardized_inputs)

Outputs:
  standardized_inputs/analysis_ready/   cleaned metadata, COI, relatedness
  standardized_inputs/dropped/          samples not in metadata
  standardized_inputs/qc/               duplicate IDs, missing dates, flag mismatches"""
    ),
    "build": (
        cmd_build,
        "Build patient-level infection sequences and recurrence pair tables",
        """\
Build patient-level infection sequences and recurrence pair tables.

Options:
  --analysis_dir  DIR   Path to standardized_inputs/analysis_ready/ (default: standardized_inputs/analysis_ready)
  --out_dir       DIR   Output directory (default: patient_level_tables)

Outputs:
  patient_infection_sequence.tsv
  baseline_pairs_with_relatedness.tsv
  consecutive_pairs_with_relatedness.tsv
  baseline_samples.tsv
  summary_by_trial_site.tsv

NOTE: for exclude_substituted_first to work in pair_mode=baseline in the
classify step, baseline_pairs_with_relatedness.tsv must carry a
baseline_time_point column (the true time_point label of whichever sample was
used as each patient's baseline, e.g. "D0" or "DRP1" if substituted). If this
column is absent, exclude_substituted_first is a no-op in baseline mode
(consecutive mode has no such dependency)."""
    ),
    "prep": (
        cmd_prep,
        "Prepare Pv3Rs-ready input TSVs from patient-level tables",
        """\
Prepare Pv3Rs-ready input TSVs from patient-level tables.

Options:
  --pis      PATH   patient_infection_sequence.tsv
  --geno     PATH   mhap_filtered.tsv (combined long-format genotype file)
  --out_dir  DIR    Root output directory (default: vivaxrec_inputs)

Outputs per trial found in the data:
  master_samples.tsv
  per_sample_locus_alleles.tsv
  d0_dr_pairs.tsv
  allele_freq_from_D0.tsv"""
    ),
    "analyse": (
        cmd_analyse,
        "Patient-level analyses: same-clone relapses, longest relapses, trajectories",
        """\
Patient-level analyses from classification outputs. Does not re-run Pv3Rs.

Three analyses are produced:

  same_clone_relapses  (consecutive pairs)
    Consecutive relapse pairs where DCiFER r >= r_thresh — suggesting the same
    clone persisting across multiple recurrences. Per-patient and per-arm summaries.

  longest_relapses  (baseline pairs)
    Relapses sorted by days since D0. Flags those exceeding day_thresh.
    Per-patient and per-arm summaries.

  patient_trajectories  (baseline and/or consecutive)
    Full per-patient recurrence sequence with classifications. Per-patient
    summary of relapse/reinfection/uncertain counts and trajectory patterns.

Options:
  --bl_dir      DIR    Baseline classification output directory
  --con_dir     DIR    Consecutive classification output directory (required)
  --out_dir     DIR    Output directory (default: same as con_dir)
  --master_dir  DIR    Directory containing <trial>/master_samples.tsv (for no-recurrence counts)
  --trial       STR    Trial name for labelling
  --rel_thresh  FLOAT  P(L) threshold for relapse (default: 0.75)
  --rei_thresh  FLOAT  P(L) threshold for reinfection (default: 0.25)
  --r_thresh    FLOAT  DCiFER r threshold for same-clone definition (default: 0.95)
  --day_thresh  INT    Days threshold for long relapse flag (default: 300)

Example:
  vivaxrec analyse \\
    --bl_dir  /data/recurrence-classification/mytrial/baseline-to-recurrence \\
    --con_dir /data/recurrence-classification/mytrial/consecutive-pairs \\
    --trial   mytrial \\
    --r_thresh 0.95 --day_thresh 300"""
    ),
    "summarise": (
        cmd_summarise,
        "Summarise classification outputs without re-running Pv3Rs",
        """\
Summarise per_infection_classifications TSV files into per-arm tables.
Reads all available classification files (full, first_rec, mono_only) and
writes separate summary TSVs for each. Does not re-run Pv3Rs.

Each summary contains:
  - P(L) threshold breakdown (relapse >= 0.75, reinfection <= 0.25, uncertain)
  - Pv3Rs top-call breakdown (L, I, C counts per arm)
  - Mean and median P(L), P(C), P(I) per arm

Options:
  --bl_dir      DIR   Baseline classification output directory
  --con_dir     DIR   Consecutive classification output directory
  --master_dir  DIR   Directory containing <trial>/master_samples.tsv (from prep step)
  --trial       STR   Trial name for labelling
  --out_dir     DIR   Output directory for summary TSVs (default: same as input)
  --mono_emoi_th FLOAT  eMOI threshold for monoclonal D0 denominator (default: 1.1)
  --d0_by_sample      Count D0 denominator by unique sample_id rather than patient_id
                      (non-default; use if multiple D0 samples per patient should each count separately)

Example:
  vivaxrec summarise \\
    --bl_dir     /data/recurrence-classification/mytrial/baseline-to-recurrence \\
    --con_dir    /data/recurrence-classification/mytrial/consecutive-pairs \\
    --master_dir /data/recurrence-classification \\
    --trial      mytrial"""
    ),
    "classify": (
        cmd_classify,
        "Run recurrence classification (DCiFER + Pv3Rs)",
        """\
Run recurrence classification (DCiFER + Pv3Rs).

Flags are passed as space-separated key value pairs (not --flag style):

Required:
  trial           TRIAL_NAME          Trial name (must match trial column in patient_infection_sequence.tsv)

Common:
  pair_mode       baseline|consecutive   Which pair type to classify (default: baseline)
  input_dir       DIR                    patient_level_tables directory
  geno            PATH                   mhap_filtered.tsv
  master_dir      DIR                    Root directory containing <trial>/master_samples.tsv
  out_dir         DIR                    Output directory for this run

Prior / model:
  max_recrud_days INT          Days threshold for recrudescence reassignment (default: 45)
  window_days     INT INT      Window boundaries in days, e.g. 35 120 (default: 35 120)
  pC              F F F        Per-window recrudescence priors (default: 0.08 0.01 0.005)
  pL              F F F        Per-window relapse priors (default: 0.25 0.70 0.40)
  lambda_foi_w    FLOAT        Daily reinfection rate (default: 0.002)
  pq_discount     FLOAT        Relapse prior multiplier for PQ arms in windows 2-3 (default: 0.19)
  site            SITE [SITE]  Restrict to one or more sites (space-separated; default: all sites)
  pq_arms         ARM [ARM]    PQ treatment arm names (space-separated)

Drug-specific window/prior overrides (CHANGES Pv3Rs posteriors — use when arms
use drugs with meaningfully different half-lives, e.g. AL vs CQ). Any arm not
listed under al_arms/cq_arms uses the global window_days/pC/pL/max_recrud_days.
Composes independently with pq_arms — e.g. a patient in ALPQ gets both the
AL window shape AND the PQ relapse-prior discount at once:
  al_arms              ARM [ARM]   Arm names to receive the AL-specific overrides below
  cq_arms              ARM [ARM]   Arm names to receive the CQ-specific overrides below
  window_days_al       INT INT     w1 w2 override for al_arms (default: falls back to window_days)
  window_days_cq       INT INT     w1 w2 override for cq_arms (default: falls back to window_days)
  pC_al / pC_cq        F F F       Per-window recrudescence prior override for al_arms/cq_arms
  pL_al / pL_cq        F F F       Per-window relapse prior override for al_arms/cq_arms
  max_recrud_days_al   INT         Recrudescence reassignment threshold override for al_arms
  max_recrud_days_cq   INT         Recrudescence reassignment threshold override for cq_arms

Reporting-only arm grouping (does NOT affect classification — only adds an
additional summary_by_arm_group*.tsv alongside the existing per-arm summary).
Arms not listed keep their own name as their group label:
  arm_group   "ARM:GROUP" ["ARM:GROUP" ...]
              e.g. OPRA:   arm_group "ALPQ:PQ" "CQPQ:PQ"      (AL, CQ stay separate)
              e.g. IMPROV: arm_group "7-day:PQ" "14-day:PQ"   (Control stays separate)

Substituted-baseline reporting exclusion. Some patients are missing a true D0
but have a genotyped first recurrence (R1) and a directly consecutive second
recurrence (R2); upstream processing substitutes R1 as a stand-in baseline so
the R1->R2 interval can still be analysed. This is treatment-equivalent where
every recurrence gets the same regimen (e.g. IMPROV — leave off), but NOT
treatment-equivalent where only the true enrolment treatment is directly
observed (e.g. OPRA, where later doses are unsupervised — turn on):
  exclude_substituted_first   on|off   Suppress the first interval from reporting
                                       for patients with a substituted (non-true-
                                       D0) baseline (default: off). Does NOT
                                       remove the substituted episode from Pv3Rs's
                                       joint model — later intervals for that
                                       patient are unaffected.

Subset analyses (each triggers an independent Pv3Rs re-run on the restricted data):
  first_rec_only  on|off       Re-run on D0 + first recurrence only (default: off)
  mono_only       on|off       Re-run on patients with monoclonal D0 only (default: off)
  mono_emoi_th    FLOAT        eMOI threshold for monoclonal definition (default: 1.1; accounts for MOIRE rounding)

Example — OPRA, consecutive: PQ arms grouped (AL/CQ stay separate), AL/CQ
half-life-specific windows, and substituted-baseline exclusion on:
  vivaxrec classify trial opra pair_mode consecutive \\
    input_dir /data/patient_level_tables \\
    geno /data/mhap_filtered.tsv \\
    master_dir /data/recurrence-classification \\
    out_dir /data/recurrence-classification/opra/consecutive-pairs \\
    pq_arms "ALPQ" "CQPQ" \\
    arm_group "ALPQ:PQ" "CQPQ:PQ" \\
    al_arms "AL" "ALPQ" window_days_al 14 60 max_recrud_days_al 14 \\
    cq_arms "CQ" "CQPQ" window_days_cq 35 120 max_recrud_days_cq 45 \\
    exclude_substituted_first on

Example — IMPROV, consecutive: PQ arms grouped (Control stays separate);
exclude_substituted_first left off (PQ given at every recurrence, so the
substitution is treatment-equivalent there):
  vivaxrec classify trial improv pair_mode consecutive \\
    input_dir /data/patient_level_tables \\
    geno /data/mhap_filtered.tsv \\
    master_dir /data/recurrence-classification \\
    out_dir /data/recurrence-classification/improv/consecutive-pairs \\
    pq_arms "7-day" "14-day" \\
    arm_group "7-day:PQ" "14-day:PQ\""""
    ),
    "plot": (
        cmd_plot,
        "Generate classification figures",
        """\
Generate classification figures.

Flags are passed as space-separated key value pairs (not --flag style):

  trial       TRIAL_NAME                Trial name (required)
  pair_mode   baseline|consecutive|both Which outputs to plot (default: both)
  bl_dir      DIR                       Baseline classification output directory
  con_dir     DIR                       Consecutive classification output directory
  figs_dir    DIR                       Output figures directory
  pis         PATH                      patient_infection_sequence.tsv

Example:
  vivaxrec plot trial mytrial pair_mode both \\
    bl_dir  /data/recurrence-classification/mytrial/baseline-to-recurrence \\
    con_dir /data/recurrence-classification/mytrial/consecutive-pairs \\
    figs_dir /data/recurrence-classification/mytrial/figures \\
    pis /data/patient_level_tables/patient_infection_sequence.tsv"""
    ),
    "run-all": (
        cmd_run_all,
        "Run the full pipeline end to end (standardize → build → prep → classify → plot)",
        """\
Run all five pipeline steps in sequence.

Required:
  --trials  TRIAL [TRIAL ...]   Trial name(s) to run through classify and plot

All other flags are forwarded to the relevant step. Use --flag style for Python
steps (standardize, build, prep) and key value style for R steps (classify, plot).

Example:
  vivaxrec run-all \\
    --trials mytrial \\
    --rel  /data/mhap_between_relatedness_estimate.tsv \\
    --coi  /data/mhap_effective_COI_summary.tsv \\
    --meta /data/metadata.tsv \\
    --out_dir /data/standardized_inputs"""
    ),
}


def main():
    parser = argparse.ArgumentParser(
        prog="vivaxrec",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", metavar="{" + ",".join(SUBCOMMANDS) + "}")

    for name, (fn, help_text, description) in SUBCOMMANDS.items():
        sub.add_parser(
            name,
            help=help_text,
            description=description,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )

    # parse_known_args lets flags meant for the underlying script (--rel, --coi etc.)
    # pass through without argparse rejecting them as unrecognised
    args, extra = parser.parse_known_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    fn, _, _ = SUBCOMMANDS[args.command]
    fn(extra)


if __name__ == "__main__":
    main()
