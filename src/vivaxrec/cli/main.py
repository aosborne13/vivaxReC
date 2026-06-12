#!/usr/bin/env python3
"""
vivaxrec — Plasmodium vivax recurrence classification pipeline.

Steps:
  1. standardize   Clean raw input files and split into analysis-ready, dropped, and QC outputs
  2. build         Build patient-level infection sequences and recurrence pair tables
  3. prep          Format data for Pv3Rs (produces master_samples.tsv etc.)
  4. classify      Run recurrence classification (DCiFER + Pv3Rs)
  5. plot          Generate classification figures
  6. run-all       Run the full pipeline end to end (steps 1–5)

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

SCRIPTS_DIR = Path(__file__).parent / "commands"

R_CMD  = "Rscript"
PY_CMD = sys.executable


def _run(cmd, **kwargs):
    print(f"[vivaxrec] Running: {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, **kwargs)
    if result.returncode != 0:
        sys.exit(result.returncode)


def _get_extra(args):
    """Strip a leading '--' separator if present — users shouldn't need to type it."""
    extra = args.extra or []
    if extra and extra[0] == "--":
        extra = extra[1:]
    return extra


def cmd_standardize(args):
    script = SCRIPTS_DIR / "standardize_inputs.py"
    _run([PY_CMD, str(script)] + _get_extra(args))


def cmd_build(args):
    script = SCRIPTS_DIR / "build_patient_level_tables.py"
    _run([PY_CMD, str(script)] + _get_extra(args))


def cmd_prep(args):
    script = SCRIPTS_DIR / "prep_pv3rs_inputs.py"
    _run([PY_CMD, str(script)] + _get_extra(args))


def cmd_classify(args):
    script = SCRIPTS_DIR / "classify_recurrences.R"
    _run([R_CMD, str(script)] + _get_extra(args))


def cmd_plot(args):
    script = SCRIPTS_DIR / "plot_recurrence_classifications.R"
    _run([R_CMD, str(script)] + _get_extra(args))


def cmd_run_all(args):
    """Run the full pipeline end to end (steps 1–5).

    --trials is required: pass one or more trial names matching the trial
    column in patient_infection_sequence.tsv.

    Example:
      vivaxrec run-all -- --trials mytrial --meta /data/metadata.tsv \\
        --rel   /data/mhap_between_relatedness_estimate.tsv \\
        --coi   /data/mhap_effective_COI_summary.tsv \\
        --out_dir /data/output
    """
    extra = _get_extra(args)

    # Extract --trials values; remaining args are passed to each step
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
        print("  Example: vivaxrec run-all -- --trials mytrial --meta /data/metadata.tsv")
        sys.exit(1)

    print(f"\n[vivaxrec] === Step 1/5: standardize ===")
    _run([PY_CMD, str(SCRIPTS_DIR / "standardize_inputs.py")] + remaining)

    print(f"\n[vivaxrec] === Step 2/5: build ===")
    _run([PY_CMD, str(SCRIPTS_DIR / "build_patient_level_tables.py")] + remaining)

    print(f"\n[vivaxrec] === Step 3/5: prep ===")
    _run([PY_CMD, str(SCRIPTS_DIR / "prep_pv3rs_inputs.py")] + remaining)

    for trial in trials:
        for pair_mode in ["baseline", "consecutive"]:
            print(f"\n[vivaxrec] === Step 4/5: classify — trial={trial} pair_mode={pair_mode} ===")
            _run([R_CMD, str(SCRIPTS_DIR / "classify_recurrences.R"),
                  "trial", trial, "pair_mode", pair_mode] + remaining)

    for trial in trials:
        print(f"\n[vivaxrec] === Step 5/5: plot — trial={trial} ===")
        _run([R_CMD, str(SCRIPTS_DIR / "plot_recurrence_classifications.R"),
              "trial", trial] + remaining)

    print(f"\n[vivaxrec] ✓ Pipeline complete for trial(s): {', '.join(trials)}")


def main():
    parser = argparse.ArgumentParser(
        prog="vivaxrec",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    commands = [
        ("standardize", cmd_standardize,
         "Clean raw input files; produce analysis-ready, dropped, and QC outputs"),
        ("build",       cmd_build,
         "Build patient-level infection sequences and recurrence pair tables"),
        ("prep",        cmd_prep,
         "Prepare Pv3Rs-ready input TSVs from patient-level tables"),
        ("classify",    cmd_classify,
         "Run recurrence classification (DCiFER + Pv3Rs)"),
        ("plot",        cmd_plot,
         "Generate classification figures"),
        ("run-all",     cmd_run_all,
         "Run the full pipeline end to end (standardize → build → prep → classify → plot)"),
    ]

    for name, fn, help_text in commands:
        p = sub.add_parser(name, help=help_text)
        p.add_argument("extra", nargs=argparse.REMAINDER,
                       help="Additional flags passed directly to the underlying script")
        p.set_defaults(func=fn)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
