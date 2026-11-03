#!/usr/bin/env Rscript
# =============================================================================
# classify_recurrences.R
#
# Recurrence classification combining DCiFER (pairwise IBD) and Pv3Rs
# (probabilistic joint inference). Allele-sharing statistics are computed
# for all pairs and retained as supplementary descriptive output.
#
# ── Input files ───────────────────────────────────────────────────────────────
# All under analysis directory:
#
#   input-files/patient_level_tables/
#     patient_infection_sequence.tsv          (eMOI source + patient linkage)
#     baseline_pairs_with_relatedness.tsv     (pair_mode baseline)
#     consecutive_pairs_with_relatedness.tsv  (pair_mode consecutive)
#
#   input-files/standardized_inputs/
#     mhap_filtered.tsv                       (raw genotypes; can cover multiple trials)
#
#   output-analysis/                          (written by prep_pv3rs_inputs.py)
#     <trial>/master_samples.tsv
#
# ── Output ────────────────────────────────────────────────────────────────────
#   output-analysis/baseline-to-recurrence/   (pair_mode baseline)
#   output-analysis/consecutive-pairs/        (pair_mode consecutive)
#     per_infection_classifications.tsv
#     summary_by_arm.tsv
#     summary_by_arm_group.tsv                (only if arm_group supplied)
#     sanity_check.txt
#     checks_emoi_missing.tsv                 (only if samples lack eMOI)
#
# ── CLI flags ─────────────────────────────────────────────────────────────────
#   trial           trial name matching the trial column in patient_infection_sequence.tsv (case-insensitive)
#   pair_mode       baseline | consecutive
#   input_dir       path to patient_level_tables/
#   geno            path to combined mhap_filtered.tsv (both trials)
#   master_dir      path containing <trial>/master_samples.tsv
#   out_dir         exact output directory for results
#   pprior          numeric  (default 0.5)
#   pdet            numeric  (default 0.9)
#   strict          strict | loose
#   pv3rs           on | off
#   af_model        d0 | pair | mix
#   alpha           numeric  (mix weight, default 0.15)
#   eps             numeric  (smoothing multiplier, default 0.05)
#   maxrec          integer  (Pv3Rs max recurrences per patient, default 6)
#   max_recrud_days numeric  (post-hoc C->L reassignment threshold, default 45)
#
# ── Independent pC / pL window schedules (primary prior system) ───────────────
#   Based on P. vivax biology, per Supplementary Table 1.
#   Defaults derived from published P. vivax radical cure trial data and should
#   be reviewed and adjusted to reflect the epidemiology of the setting being analysed.
#
#   *** pC and pL use INDEPENDENT window-day schedules ***
#   Recrudescence plausibility (pC) uses the blood-stage partner drug's
#   terminal half-life, so pC's window boundaries are drug-specific
#   (see al_arms/cq_arms below). 
#   Relapse timing (pL) is governed by hypnozoite
#   reactivation biology - no link to the blood-stage
#   partner drug - pL ALWAYS follows a single schedule
#   (window_days / pL below), identically for every arm, PQ-discounted or not.
#
#     CQ efficacy example (Getachew et al. 2015): 0.7% early, 8% late treatment failure
#     PQ efficacy example (Taylor et al. 2019 Lancet): ~5.3x reduction in relapse
#     AL efficacy example (Byakika-Kibwika et al. 2010, PMC2817784): lumefantrine
#     terminal half-life ~3-6 days (range 4-10 days in patients), vs CQ's 20-60 days
#
#   window_days     w1 w2  (default: 35 120) - pC AND pL default schedule;
#                   for pC specifically, this is OVERRIDDEN per-drug via
#                   window_days_al / window_days_cq below.
#     Window 1 ≤w1:   recrudescence possible, early relapse possible
#     Window 2 w1-w2: recrudescence very unlikely, primary relapse window
#     Window 3 >w2:   recrudescence implausible
#
#   pC              pC_w1 pC_w2 pC_w3  (default: 0.08 0.01 0.005)
#   pL              pL_w1 pL_w2 pL_w3  (default: 0.25 0.70 0.40)
#     Note: pL values shown are for non-PQ arms. PQ arms receive pL*pq_discount
#     in windows 2 and 3 (window 1 undiscounted — PQ course still being taken).
#     This discount is applied on top of the global pL window (see above)
#
#   lambda_foi_w    daily reinfection rate, uncapped  (default: 0.002/day)
#     pI grows as 1-exp(-lambda_foi_w * t); no ceiling — reinfection possible
#     at any time regardless of treatment arm
#
#   pq_arms         names of PQ treatment arms (space-separated; used for pq_discount)
#                   e.g.: pq_arms "7-day" "14-day"
#                   If omitted, no arms receive the PQ discount.
#                   Combines with al_arms/cq_arms independently — e.g. a patient
#                   in ALPQ gets BOTH the AL-specific pC window (from al_arms)
#                   AND the PQ relapse-prior discount (from pq_arms) simultaneously.
#
#   armtime         legacy 5-value log-normal prior (back-compat only)
#                   lambda_foi bump_center bump_sdlog pC_base pI_base
#
# ── Drug-specific pC window/prior overrides (pC ONLY — see note above) ────────
#   These OVERRIDE window_days / pC / max_recrud_days, for the pC (recrudescence)
#   calculation ONLY, for patients in the named arms, to account for drugs with
#   materially different half-lives (e.g. AL's few-day terminal half-life vs
#   CQ's 20-60 day half-life). Any arm not listed under al_arms or cq_arms uses
#   the global window_days/pC/max_recrud_days values unchanged for pC.
#   pL is NEVER affected by al_arms/cq_arms — every arm's pL always follows the
#   single global window_days/pL schedule. This changes Pv3Rs posteriors!
#   Composes independently with pq_arms (see above).
#
#   al_arms              names of arms to receive the AL-specific pC overrides below
#                        e.g.: al_arms "AL" "ALPQ"
#   cq_arms               names of arms to receive the CQ-specific pC overrides below
#                        e.g.: cq_arms "CQ" "CQPQ"
#   window_days_al    w1 w2 override for al_arms' pC windows (default: falls back to global window_days)
#   window_days_cq    w1 w2 override for cq_arms' pC windows (default: falls back to global window_days)
#   pC_al / pC_cq        pC_w1 pC_w2 pC_w3 override for al_arms / cq_arms
#   max_recrud_days_al   recrudescence reassignment threshold override for al_arms
#   max_recrud_days_cq   recrudescence reassignment threshold override for cq_arms
#
# ── Reporting-only arm grouping ─────────────────────────────────────────────────
#   Purely additive — does NOT affect Pv3Rs computation, allele frequencies, r_bg,
#   or any classification. Only relabels treatment_arm for an additional
#   summary_by_arm_group*.tsv output, alongside the existing per-arm summary.
#
#   arm_group       "raw_arm:group_label" pairs, space-separated
#                   e.g. OPRA: arm_group "ALPQ:PQ" "CQPQ:PQ"   (AL, CQ stay separate)
#                   e.g. IMPROV: arm_group "7-day:PQ" "14-day:PQ"   (Control stays separate)
#                   Arms not listed keep their own name as their group label.
#
# ── Substituted-baseline handling (reporting-level exclusion) ──────────────────
#   OPTIONAL SETTING - DO NOT USE BY DEFAULT - Some patients are missing a true Day 0 sample but have a genotyped first
#   recurrence (R1) and a consecutive second recurrence (R2); upstream
#   processing substitutes R1 as a stand-in "baseline" so the R1->R2 interval can
#   still be analysed DEPENING ON THE SITUATION. This is treatment-equivalent in trials where every
#   recurrence receives the same regimen (e.g. where PQ is given at every
#   recurrence), but NOT treatment-equivalent where only the true enrolment
#   treatment is directly observed and later doses are unsupervised -
#   in that setting we cannot know whether the substituted "baseline" visit itself
#   received supervised treatment.
#
#   exclude_substituted_first   on | off  (default: off)
#     When on: the first interval for any patient with a substituted (non-true-D0)
#     baseline is dropped from reporting output (per_infection_classifications.tsv,
#     summary_by_arm.tsv). The substituted episode's genotype is NOT removed from
#     Pv3Rs's joint model — later intervals for that same patient (R2->R3 etc.)
#     are retained and still benefit from the full episode context.
#     Detection: baseline mode requires a d0_time_point column populated from an
#     upstream baseline_time_point field (see build_patient_level_tables.py); a
#     patient is flagged substituted if that value is not "D0". Consecutive mode
#     requires no new upstream column — it uses the existing from_time_point on
#     step_no==1 rows, on the same "not D0" logic.
#
# ── Pv3Rs stepwise fallback ───────────────────────────────────────────────────
#   When sum(MOIs) > 8 for a patient, Pv3Rs runs each consecutive pair
#   independently rather than returning all-NA. Results tagged
#   pv3rs_note = "stepwise_fallback" in output.
#
# ── Example run commands ──────────────────────────────────────────────────────
#   BASE=/path/to/your/analysis/directory
#   SCRIPTS=/path/to/vivaxReC/src/vivaxrec/commands
#   TRIAL=mytrial   # must match the trial column in patient_infection_sequence.tsv
#
#   # Baseline pairs
#   Rscript $SCRIPTS/classify_recurrences.R \
#     trial $TRIAL pair_mode baseline max_recrud_days 45 \
#     input_dir  $BASE/input-files/patient_level_tables \
#     geno       $BASE/input-files/standardized_inputs/mhap_filtered.tsv \
#     master_dir $BASE/output-analysis/recurrence-classification \
#     out_dir    $BASE/output-analysis/recurrence-classification/$TRIAL/baseline-to-recurrence
#
#   # Consecutive pairs
#   Rscript $SCRIPTS/classify_recurrences.R \
#     trial $TRIAL pair_mode consecutive max_recrud_days 45 \
#     input_dir  $BASE/input-files/patient_level_tables \
#     geno       $BASE/input-files/standardized_inputs/mhap_filtered.tsv \
#     master_dir $BASE/output-analysis/recurrence-classification \
#     out_dir    $BASE/output-analysis/recurrence-classification/$TRIAL/consecutive-pairs
#
#   # OPRA, consecutive: group PQ arms (AL/CQ stay separate), AL/CQ pC half-life
#   # windows ONLY (pL stays identical across all arms — see note above), and
#   # exclude substituted-baseline first intervals (OPRA-specific, since only
#   # the true D0 treatment is directly observed there)
#   Rscript $SCRIPTS/classify_recurrences.R \
#     trial opra pair_mode consecutive \
#     input_dir  $BASE/input-files/patient_level_tables \
#     geno       $BASE/input-files/standardized_inputs/mhap_filtered.tsv \
#     master_dir $BASE/output-analysis/recurrence-classification \
#     out_dir    $BASE/output-analysis/recurrence-classification/opra/consecutive-pairs \
#     pq_arms "ALPQ" "CQPQ" \
#     arm_group "ALPQ:PQ" "CQPQ:PQ" \
#     al_arms "AL" "ALPQ" window_days_al 14 60 pC_al 0.08 0.01 0.005 max_recrud_days_al 14 \
#     cq_arms "CQ" "CQPQ" window_days_cq 35 120 pC_cq 0.08 0.01 0.005 max_recrud_days_cq 45 \
#     exclude_substituted_first on
#
#   # IMPROV, consecutive: group PQ arms (Control stays separate); substitution
#   # exclusion left off (default) since PQ is given at every recurrence there,
#   # making the substitution treatment-equivalent
#   Rscript $SCRIPTS/classify_recurrences.R \
#     trial improv pair_mode consecutive \
#     input_dir  $BASE/input-files/patient_level_tables \
#     geno       $BASE/input-files/standardized_inputs/mhap_filtered.tsv \
#     master_dir $BASE/output-analysis/recurrence-classification \
#     out_dir    $BASE/output-analysis/recurrence-classification/improv/consecutive-pairs \
#     pq_arms "7-day" "14-day" \
#     arm_group "7-day:PQ" "14-day:PQ"
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(purrr); library(rlang)
  library(dcifer)
})

# ── Helpers ───────────────────────────────────────────────────────────────────
`%||%` <- function(x, y) if (is.null(x)) y else x

parse_cli <- function(args) {
  spec <- c(
    trial      = 1,   # any trial name matching a value in the trial column of patient_infection_sequence.tsv
    pair_mode  = 1,   # baseline | consecutive
    input_dir  = 1,   # dir containing patient_infection_sequence.tsv + pairs TSVs
    geno       = 1,   # path to combined mhap_filtered.tsv (both trials)
    master_dir = 1,   # dir containing <trial>/master_samples.tsv (prep output root)
    out_dir    = 1,   # exact directory for output files
    pprior     = 1,
    pdet       = 1,
    strict     = 1,   # strict | loose
    pv3rs      = 1,   # on | off
    af_model   = 1,   # d0 | pair | mix
    alpha      = 1,
    eps        = 1,
    cap            = 1,   # kept for back-compat, unused
    max_recrud_days = 1,  # post-hoc C->L/I reassignment threshold in days (default: 45)
    maxrec         = 1,   # maximum number of recurrences per patient passed to Pv3Rs
    armtime        = 5,   # legacy: lambda_foi, bump_center, bump_sdlog, pC_base, pI_base
    # ── Global window/prior schedule — used by pL ALWAYS, and by pC unless
    # overridden per-drug via window_days_al / window_days_cq below ─────
    # Window 1: 0–w1 days     — recrudescence possible, short-latency relapse possible
    # Window 2: w1+1–w2 days  — recrudescence very unlikely, primary relapse window
    # Window 3: >w2 days      — recrudescence implausible, reinfection accumulates
    window_days    = 2,   # w1 w2  (default: 35 120)
    pC             = 3,   # pC_w1 pC_w2 pC_w3  (default: 0.08 0.01 0.005)
    pL             = 3,   # pL_w1 pL_w2 pL_w3  (default: 0.25 0.70 0.40) — ALWAYS applies, every arm
    # pI_cap removed: reinfection risk is uncapped (grows freely with lambda_foi)
    lambda_foi_w   = 1,   # daily force of infection for pI growth (default: 0.002)
    pq_discount    = 1,   # multiplier on pL for PQ arms in windows 2+3 (default: 0.19)
    pq_arms        = -1,  # names of PQ treatment arms (variable length; space-separated)
    first_rec_only = 1,   # on | off — re-run Pv3Rs on first-recurrence patients only
    mono_only      = 1,   # on | off — re-run Pv3Rs on monoclonal-D0 patients only (D0 eMOI <= threshold)
    mono_emoi_th   = 1,   # eMOI threshold for monoclonal definition (default: 1.1; accounts for MOIRE rounding)
    site           = -1,  # site name(s) to restrict analysis to (variable length; space-separated)
    # ── Reporting-only arm grouping (does not affect Pv3Rs computation) ────────
    arm_group      = -1,  # "raw_arm:group_label" pairs, e.g. arm_group "ALPQ:PQ" "CQPQ:PQ"
    # ── Drug-specific pC window/prior overrides (pC ONLY — never pL) ───────────
    al_arms             = -1,  # arm names treated with AL (short half-life) pC overrides
    cq_arms             = -1,  # arm names treated with CQ (long half-life) pC overrides
    window_days_al   = 2,   # w1 w2 override for al_arms' pC windows (falls back to global window_days)
    window_days_cq   = 2,   # w1 w2 override for cq_arms' pC windows (falls back to global window_days)
    pC_al                = 3,  # pC_w1 pC_w2 pC_w3 override for al_arms
    pC_cq                = 3,  # pC_w1 pC_w2 pC_w3 override for cq_arms
    max_recrud_days_al  = 1,   # recrudescence reassignment threshold override for al_arms
    max_recrud_days_cq  = 1,   # recrudescence reassignment threshold override for cq_arms
    # ── Substituted-baseline reporting exclusion ────────────────────────────────
    exclude_substituted_first = 1   # on | off — suppress reporting of the first
                                     # interval for patients with a substituted
                                     # (non-true-D0) baseline (default: off)
  )
  out <- list(); i <- 1L
  while (i <= length(args)) {
    key <- tolower(args[[i]])
    if (!key %in% names(spec)) stop("Unknown flag: '", key, "'")
    n <- spec[[key]]
    if (n == -1L) {
      # Variable-length: consume all remaining values until the next known flag
      vals <- character(0)
      j <- i + 1L
      while (j <= length(args) && !tolower(args[[j]]) %in% names(spec)) {
        vals <- c(vals, args[[j]]); j <- j + 1L
      }
      out[[key]] <- vals; i <- j
    } else {
      if (i + n > length(args)) stop("Flag '", key, "' expects ", n, " value(s)")
      out[[key]] <- args[(i + 1L):(i + n)]
      i <- i + n + 1L
    }
  }
  out
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))

trial       <- tolower(cli$trial     %||% stop("'trial' flag is required"))
pair_mode   <- tolower(cli$pair_mode %||% "baseline")
input_dir   <- cli$input_dir  %||% "patient_level_tables"
geno_fp_arg <- cli$geno       %||% file.path("standardized_inputs", "mhap_filtered.tsv")
master_dir  <- cli$master_dir %||% "output-analysis"
out_dir     <- cli$out_dir    %||% file.path(
  "output-analysis",
  if (pair_mode == "baseline") "baseline-to-recurrence" else "consecutive-pairs"
)
pprior      <- as.numeric(cli$pprior %||% 0.5)
pdet  <- as.numeric(cli$pdet   %||% 0.9)
strict_flag <- tolower(cli$strict    %||% "loose")
pv3rs_flag  <- tolower(cli$pv3rs    %||% "on")
af_model    <- tolower(cli$af_model  %||% "mix")
alpha       <- as.numeric(cli$alpha  %||% 0.15)
eps         <- as.numeric(cli$eps    %||% 0.05)
maxrec          <- as.integer(cli$maxrec           %||% 6L)
max_recrud_days <- as.numeric(cli$max_recrud_days  %||% 45)

use_armtime <- !is.null(cli$armtime)
if (use_armtime) {
  lambda_foi  <- as.numeric(cli$armtime[1])
  bump_center <- as.numeric(cli$armtime[2])
  bump_sdlog  <- as.numeric(cli$armtime[3])
  pC_base     <- as.numeric(cli$armtime[4])
  pI_base     <- as.numeric(cli$armtime[5])
} else {
  lambda_foi <- bump_center <- bump_sdlog <- pC_base <- pI_base <- NA_real_
}

# ── Global window/prior parameters ────────────────────────────────────────────
# Defaults using P. vivax biology and published Ethiopian trial data.
# Users should review and adjust these values to reflect their own study setting.
#
# w1/w2 and pL_w below are the ONE AND ONLY schedule pL ever uses, for every
# arm, PQ or not — there is no per-drug pL override anywhere in this script.
# w1/w2 and pC_w are pC's DEFAULT schedule, used for any arm not listed under
# al_arms/cq_arms; arms that ARE listed use their window_days_*/pC_*
# override instead (see resolve_drug_params() below).
#
# pq_discount: multiplier on pL for PQ arms in windows 2 & 3 of pL's (global)
# schedule — never applied against a drug-specific window, since pL has none.
#   default 0.19 = 1/5.3 from published ~5-fold reduction in relapse with PQ
use_window_prior <- !is.null(cli$window_days)
if (use_window_prior) {
  w1           <- as.numeric(cli$window_days[1])  # default 35
  w2           <- as.numeric(cli$window_days[2])  # default 120
  pC_w         <- as.numeric(cli$pC)              # 3 values: w1 w2 w3
  pL_w         <- as.numeric(cli$pL)              # 3 values: w1 w2 w3 — ALWAYS applies
  lambda_foi_w <- as.numeric(cli$lambda_foi_w %||% 0.002)
  pq_discount  <- as.numeric(cli$pq_discount  %||% 0.19)
} else {
  w1 <- 35; w2 <- 120
  pC_w     <- c(0.08, 0.01, 0.005)
  pL_w     <- c(0.25, 0.70, 0.40)
  lambda_foi_w <- 0.002
  pq_discount  <- 0.19
}

# pq_arms_flag: treatment arm names to receive PQ discount in windows 2-3.
# Supplied via CLI as: pq_arms "7-day" "14-day"
# If not supplied, no discount is applied.
pq_arms_flag     <- cli$pq_arms %||% character(0)

# ── Reporting-only arm grouping ────────────────────────────────────────────────
# Parsed as "raw_arm:group_label" tokens. Purely additive — never used in any
# Pv3Rs computation, only to label output rows for an additional summary table.
arm_group_flag <- cli$arm_group %||% character(0)
arm_group_map  <- character(0)
if (length(arm_group_flag) > 0) {
  parts <- strsplit(arm_group_flag, ":", fixed = TRUE)
  bad   <- vapply(parts, length, integer(1)) != 2L
  if (any(bad))
    stop("arm_group entries must be 'raw_arm:group_label', got: ",
         paste(arm_group_flag[bad], collapse = ", "))
  arm_group_map <- setNames(vapply(parts, `[`, character(1), 2),
                            vapply(parts, `[`, character(1), 1))
  message(sprintf(">> arm_group mapping: %s",
                  paste(names(arm_group_map), arm_group_map, sep="->", collapse=", ")))
}
map_arm_group <- function(raw_arm) {
  if (raw_arm %in% names(arm_group_map)) arm_group_map[[raw_arm]] else raw_arm
}

# ── Drug-specific pC window/prior overrides (pC ONLY) ──────────────────────────
# Must be resolved AFTER w1/w2/pC_w/max_recrud_days are set above, since any
# per-drug value not explicitly overridden falls back to those globals.
# Composes independently with pq_arms: a patient's arm can appear in BOTH
# al_arms/cq_arms (controls pC window shape) and pq_arms (controls pL's
# relapse-prior discount) at once — e.g. ALPQ, CQPQ.
#
# NOTE: there is deliberately no pL_al/pL_cq or window_days_pL_al/
# window_days_pL_cq anywhere in this script. pL always uses the single
# global w1/w2/pL_w schedule, for every arm — see header note above.
al_arms_flag <- cli$al_arms %||% character(0)
cq_arms_flag <- cli$cq_arms %||% character(0)

resolve_drug_params <- function(suffix, arms_flag) {
  if (!length(arms_flag)) return(NULL)
  wdC <- cli[[paste0("window_days_", suffix)]]
  pc  <- cli[[paste0("pC_", suffix)]]
  mrd <- cli[[paste0("max_recrud_days_", suffix)]]
  list(
    arms            = arms_flag,
    w1C             = if (!is.null(wdC)) as.numeric(wdC[1]) else w1,
    w2C             = if (!is.null(wdC)) as.numeric(wdC[2]) else w2,
    pC_w            = if (!is.null(pc))  as.numeric(pc)     else pC_w,
    max_recrud_days = if (!is.null(mrd)) as.numeric(mrd)    else max_recrud_days
  )
}

drug_params <- list(
  al = resolve_drug_params("al", al_arms_flag),
  cq = resolve_drug_params("cq", cq_arms_flag)
)
if (!is.null(drug_params$al))
  message(sprintf(">> AL-arm pC overrides: w1C=%g w2C=%g pC=%s max_recrud_days=%g (arms: %s) | pL: global, unaffected",
                  drug_params$al$w1C, drug_params$al$w2C,
                  paste(drug_params$al$pC_w, collapse=" "),
                  drug_params$al$max_recrud_days,
                  paste(al_arms_flag, collapse=", ")))
if (!is.null(drug_params$cq))
  message(sprintf(">> CQ-arm pC overrides: w1C=%g w2C=%g pC=%s max_recrud_days=%g (arms: %s) | pL: global, unaffected",
                  drug_params$cq$w1C, drug_params$cq$w2C,
                  paste(drug_params$cq$pC_w, collapse=" "),
                  drug_params$cq$max_recrud_days,
                  paste(cq_arms_flag, collapse=", ")))

# Subset analysis flags — each triggers a separate Pv3Rs run on the restricted dataset.
# The full-dataset run always happens first; subsets are additional outputs.
first_rec_only   <- tolower(cli$first_rec_only %||% "off")  # on | off
mono_only        <- tolower(cli$mono_only      %||% "off")  # on | off
mono_emoi_th     <- as.numeric(cli$mono_emoi_th %||% 1.1)  # raw eMOI threshold; 1.1 used to account for MOIRE rounding
site_filter      <- cli$site %||% character(0)              # site name(s); empty = no filter

# ── Substituted-baseline reporting exclusion ────────────────────────────────────
exclude_substituted_first <- tolower(cli$exclude_substituted_first %||% "off")

stopifnot(first_rec_only %in% c("on", "off"))
stopifnot(mono_only      %in% c("on", "off"))
stopifnot(exclude_substituted_first %in% c("on", "off"))

stopifnot(pair_mode %in% c("baseline", "consecutive"))
stopifnot(strict_flag %in% c("strict", "loose"))
stopifnot(pv3rs_flag  %in% c("on", "off"))
stopifnot(af_model    %in% c("d0", "pair", "mix"))

set.seed(123)
TRIAL_UP <- toupper(trial)

# ── Paths ─────────────────────────────────────────────────────────────────────
fp_pis       <- file.path(input_dir, "patient_infection_sequence.tsv")
fp_pairs_bl  <- file.path(input_dir, "baseline_pairs_with_relatedness.tsv")
fp_pairs_con <- file.path(input_dir, "consecutive_pairs_with_relatedness.tsv")
geno_fp      <- geno_fp_arg
master_fp    <- file.path(master_dir, trial, "master_samples.tsv")

# Validate all required inputs exist before doing any work
required <- c(
  fp_pis,
  if (pair_mode == "baseline") fp_pairs_bl else fp_pairs_con,
  geno_fp,
  master_fp
)
missing_fps <- required[!file.exists(required)]
if (length(missing_fps))
  stop("Required input file(s) not found:\n",
       paste(" ", missing_fps, collapse = "\n"),
       "\nCheck --input_dir / --geno_dir / --master_dir flags.")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message(">> trial: ", TRIAL_UP, " | pair_mode: ", pair_mode)
message(">> input_dir:  ", input_dir)
message(">> geno:       ", geno_fp)
message(">> master_dir: ", master_dir)
message(">> out_dir:    ", out_dir)

# ── Load patient_infection_sequence & eMOI lookup ─────────────────────────────
message(">> loading patient_infection_sequence")
pis_all <- read_tsv(fp_pis, show_col_types = FALSE)

# eMOI lookup: sample_id -> MOIRE post_effective_coi_mean
# Single source of truth for COI across the whole pipeline
emoi_lookup <- pis_all %>%
  filter(!is.na(post_effective_coi_mean)) %>%
  transmute(sample_id, emoi = post_effective_coi_mean) %>%
  distinct(sample_id, .keep_all = TRUE)

get_emoi <- function(sid, fallback = NA_real_) {
  val <- emoi_lookup$emoi[match(sid, emoi_lookup$sample_id)]
  dplyr::coalesce(val, fallback)
}

# ── Load pairs file ───────────────────────────────────────────────────────────
message(">> loading pairs (pair_mode=", pair_mode, ")")
pairs_raw <- read_tsv(
  if (pair_mode == "baseline") fp_pairs_bl else fp_pairs_con,
  show_col_types = FALSE
)
pairs <- pairs_raw %>%
  filter(toupper(trial) == TRIAL_UP, !is.na(relatedness))

if (length(site_filter) > 0) {
  pairs <- pairs %>% filter(site %in% site_filter)
  message(sprintf(">> site filter applied: %s", paste(site_filter, collapse = ", ")))
}

message(sprintf(">> pairs loaded: %d rows (trial=%s, after dropping NA relatedness)",
                nrow(pairs), TRIAL_UP))
if (nrow(pairs) == 0)
  stop("No pairs found for trial ", TRIAL_UP, " in ", pair_mode, " mode.")

# For baseline mode, check whether the upstream file carries the baseline-side
# time_point/day columns needed for substituted-baseline detection. If absent,
# exclude_substituted_first=on will be a no-op for this pair_mode with a warning
# (rather than a hard failure), since consecutive mode does not depend on this.
if (pair_mode == "baseline") {
  has_bl_tp  <- "baseline_time_point" %in% names(pairs)
  if (exclude_substituted_first == "on" && !has_bl_tp) {
    message(">> WARNING: exclude_substituted_first=on but 'baseline_time_point' column ",
            "not found in baseline_pairs_with_relatedness.tsv — substituted-baseline ",
            "rows cannot be detected in baseline mode. Add this column upstream ",
            "(build_patient_level_tables.py) or use pair_mode=consecutive instead.")
  }
}

# ── Load raw genotypes ────────────────────────────────────────────────────────
message(">> loading raw genotypes for ", TRIAL_UP)
geno_raw <- read_tsv(geno_fp, show_col_types = FALSE)

bad_alleles <- c("", "NA", ".", "na", "Na", "n/a")
n_before <- nrow(geno_raw)
geno <- geno_raw %>%
  mutate(
    allele    = trimws(as.character(allele)),
    locus     = as.character(locus),
    sample_id = as.character(sample_id)
  ) %>%
  filter(!is.na(allele), !(allele %in% bad_alleles)) %>%
  distinct(sample_id, locus, allele, .keep_all = TRUE)
n_after <- nrow(geno)
message(sprintf(">> allele rows: kept %d of %d (dropped %d missing/blank/dup)",
                n_after, n_before, n_before - n_after))

# ── Load master_samples (produced by prep_aster_inputs.py) ───────────────────
# patient_id here is external_id from patient_infection_sequence — the correct
# patient grouping key (not a string split of sample_id, which gives a different format).
samples <- read_tsv(master_fp, show_col_types = FALSE)

meta <- samples %>%
  transmute(
    sample_id,
    id            = as.character(patient_id),
    day           = suppressWarnings(as.integer(day)),
    time_point    = as.character(time_point),
    treatment_arm = as.character(treatment_arm),
    site          = as.character(site)
  )

# ── Site filter ───────────────────────────────────────────────────────────────
# When site_filter is supplied, restrict meta (and therefore dfsmp, dsmp, coi,
# and all allele frequency estimates) to samples from the specified site(s) only.
# This ensures allele frequencies are estimated from site-specific D0 samples
# and Pv3Rs only sees episodes from that site — consistent with the mono_only logic.
if (length(site_filter) > 0) {
  n_before_site <- nrow(meta)
  meta <- meta %>% filter(site %in% site_filter)
  message(sprintf(">> site filter: retained %d of %d samples (site(s): %s)",
                  nrow(meta), n_before_site, paste(site_filter, collapse = ", ")))
  if (nrow(meta) == 0)
    stop("No samples found for site(s): ", paste(site_filter, collapse = ", "))
}

# ── Format genotypes via dcifer ───────────────────────────────────────────────
dfgen    <- geno %>% transmute(smp_id = sample_id, locus, allele)
dsmp_all <- dcifer::formatDat(dfgen, svar = "smp_id", lvar = "locus", avar = "allele")

has_calls_sid <- function(lst, sid) {
  x <- lst[[sid]]
  is.list(x) && length(x) > 0 && any(vapply(x, function(v) length(v) > 0, logical(1)))
}

present <- intersect(meta$sample_id, names(dsmp_all))
present <- present[vapply(present, \(s) has_calls_sid(dsmp_all, s), logical(1))]

dfsmp <- meta %>%
  filter(.data$sample_id %in% present, !is.na(.data$day)) %>%
  arrange(.data$id, .data$day, .data$sample_id) %>%
  transmute(id = .data$id, smp_id = .data$sample_id, day = .data$day)

dsmp <- dsmp_all[dfsmp$smp_id]
stopifnot(length(dsmp) == nrow(dfsmp))
stopifnot(!any(vapply(dsmp, is.null, logical(1))))

# ── COI vector: MOIRE eMOI rounded to integer, floored at 1 ──────────────────
coi_raw <- get_emoi(dfsmp$smp_id, fallback = NA_real_)
coi     <- pmax(1L, round(coi_raw))

n_coi_na <- sum(is.na(coi))
if (n_coi_na > 0) {
  warn_path <- file.path(out_dir, "checks_emoi_missing.tsv")
  tibble(
    smp_id = dfsmp$smp_id[is.na(coi)],
    id     = dfsmp$id[is.na(coi)],
    day    = dfsmp$day[is.na(coi)]
  ) %>% write_tsv(warn_path)
  msg <- sprintf(
    ">> eMOI CHECK: %d samples have no MOIRE estimate -> falling back to COI=1. See: %s",
    n_coi_na, warn_path)
  if (strict_flag == "strict") stop(msg) else message(msg)
  coi[is.na(coi)] <- 1L
}

dfsmp$pdet <- 0.9

# ── Allele frequency computation (reusable for subset runs) ──────────────────
# Encapsulated so the mono_only run can recompute from monoclonal D0 only.
# panel_all is built from ALL samples (ensures every observed allele has an
# entry with at least eps_floor, keeping frequency tables pairwise-stable).
panel_all <- list()
for (s in dsmp) for (L in names(s))
  panel_all[[L]] <- union(panel_all[[L]], unique(names(s[[L]])))

calc_afreq <- function(sub_dsmp, sub_coi, sub_is_d0) {
  if (!any(sub_is_d0)) stop("calc_afreq: no D0 samples in subset.")
  af <- dcifer::calcAfreq(sub_dsmp[sub_is_d0], sub_coi[sub_is_d0], tol = 1e-5)
  minfr_sub  <- 1 / sum(sub_coi[sub_is_d0])
  eps_fl_sub <- max(minfr_sub * eps, 1e-12)
  af <- lapply(af, function(v) { v[v == 0] <- minfr_sub; v / sum(v) })
  # Augment with full panel so every allele in the data has a frequency entry
  af_aug <- lapply(names(panel_all), function(L) {
    alleles <- panel_all[[L]]
    v <- af[[L]]
    if (is.null(v)) {
      v <- rep(eps_fl_sub, length(alleles)); names(v) <- alleles
    } else {
      miss <- setdiff(alleles, names(v))
      if (length(miss)) v <- c(v, setNames(rep(eps_fl_sub, length(miss)), miss))
    }
    v / sum(v)
  })
  names(af_aug) <- names(panel_all)
  list(afreq = af_aug, eps_floor = eps_fl_sub)
}

# ── Allele frequencies (π) and r_bg ──────────────────────────────────────────
# When mono_only is on: subset to monoclonal D0 samples first, then compute
# allele frequencies and r_bg from that subset only.
# When mono_only is off: use all D0 samples.
is_d0 <- dfsmp$day == 0
if (!any(is_d0)) stop("No D0 samples found for ", TRIAL_UP, ".")

if (mono_only == "on") {
  d0_emoi     <- get_emoi(dfsmp$smp_id[is_d0], fallback = NA_real_)
  is_d0_use   <- is_d0
  is_d0_use[is_d0] <- !is.na(d0_emoi) & d0_emoi <= mono_emoi_th
  n_d0        <- sum(is_d0_use)
  if (n_d0 == 0) stop("No monoclonal D0 samples found (eMOI <= ", mono_emoi_th, ")")
  message(sprintf(">> subsetting to %d monoclonal D0 samples (eMOI <= %.2f) for afreq and r_bg",
                  n_d0, mono_emoi_th))
  af_mono  <- calc_afreq(dsmp, coi, is_d0_use)
  af_full  <- NULL   # not needed
} else {
  is_d0_use <- is_d0
  n_d0      <- sum(is_d0_use)
  af_full   <- calc_afreq(dsmp, coi, is_d0_use)
  af_mono   <- NULL  # not needed
  message(sprintf(">> allele frequencies estimated from %d D0 samples (all)", n_d0))
}

# afreq_base and eps_floor are the frequencies used by the current run's
# Pv3Rs and dcifer calls (make_afreq_pair, r_bg, build_pv3rs_inputs)
af_run     <- if (mono_only == "on") af_mono else af_full
afreq_base <- af_run$afreq
eps_floor  <- af_run$eps_floor

# r_bg: computed from the same D0 sample set used for afreq
message(sprintf(">> computing r_bg from %d D0 samples", n_d0))
rbg <- mean(
  dcifer::ibdDat(dsmp[is_d0_use], coi[is_d0_use], afreq_base, pval = FALSE),
  na.rm = TRUE
)
message(sprintf(">> r_bg = %.4f (n_D0 = %d)", rbg, n_d0))

# ── Pair-specific allele-frequency helper ─────────────────────────────────────
make_afreq_pair <- function(i, j) {
  li <- dsmp[[i]]; lj <- dsmp[[j]]
  loci <- intersect(names(li), names(lj))
  out  <- vector("list", length(loci)); names(out) <- loci
  for (L in loci) {
    alle_obs <- union(names(li[[L]]), names(lj[[L]]))
    v_base <- afreq_base[[L]]
    if (is.null(v_base)) {
      v_base <- rep(eps_floor, length(alle_obs)); names(v_base) <- alle_obs
    } else {
      miss <- setdiff(alle_obs, names(v_base))
      if (length(miss)) v_base <- c(v_base, setNames(rep(eps_floor, length(miss)), miss))
    }
    v_base <- v_base / sum(v_base)
    if (af_model == "d0") {
      out[[L]] <- v_base
    } else {
      v_pair <- rep(1 / length(alle_obs), length(alle_obs)); names(v_pair) <- alle_obs
      v_pair <- (v_pair + eps_floor) / sum(v_pair + eps_floor)
      out[[L]] <- if (af_model == "pair") v_pair else {
        v_mix <- (1 - alpha) * v_base + alpha * v_pair; v_mix / sum(v_mix)
      }
    }
  }
  out
}

allele_share <- function(si, sj) {
  li <- dsmp[[si]]; lj <- dsmp[[sj]]
  loci <- intersect(names(li), names(lj))
  if (!length(loci)) return(c(shared = 0L, tested = 0L))
  shared <- sum(vapply(loci, function(L) {
    any((as.numeric(li[[L]]) > 0) & (as.numeric(lj[[L]]) > 0))
  }, logical(1)))
  c(shared = as.integer(shared), tested = as.integer(length(loci)))
}

# ── Build pair table ──────────────────────────────────────────────────────────
if (pair_mode == "baseline") {
  has_bl_tp <- "baseline_time_point" %in% names(pairs)
  pair_tbl <- pairs %>%
    transmute(
      trial          = TRIAL_UP,
      pair_type      = "D0_to_DR",
      id             = external_id,                  # correct patient identifier
      d0_sample_id   = dplyr::coalesce(baseline_sample_id, sample_id1),
      d0_time_point  = if (has_bl_tp) as.character(baseline_time_point) else NA_character_,
      dr_sample_id   = current_sample_id,
      dr_day         = current_day,
      dr_time_point  = current_time_point,
      treatment_arm,
      arm_group      = vapply(treatment_arm, map_arm_group, character(1)),
      site,
      emoi_s1        = get_emoi(dplyr::coalesce(baseline_sample_id, sample_id1)),
      emoi_s2        = current_post_effective_coi_mean,
      pdet_s1        = 0.9,
      pdet_s2        = 0.9,
      dcifer_r       = relatedness,
      .s1            = dplyr::coalesce(baseline_sample_id, sample_id1),
      .s2            = current_sample_id
    ) %>%
    arrange(id, dr_day, dr_sample_id) %>%
    group_by(id) %>%
    mutate(recurrence_no = dplyr::row_number()) %>%
    ungroup() %>%
    mutate(
      # Only meaningful (and only ever TRUE) on the first interval per patient;
      # later recurrence_no values are never substituted-baseline candidates.
      is_substituted_baseline = recurrence_no == 1L &
                                !is.na(d0_time_point) &
                                d0_time_point != "D0"
    )
} else {
  pair_tbl <- pairs %>%
    transmute(
      trial               = TRIAL_UP,
      pair_type           = "DR_to_nextDR",
      id                  = external_id,             # correct patient identifier
      from_sample_id      = dplyr::coalesce(previous_sample_id, sample_id1),
      to_sample_id        = current_sample_id,
      from_day            = previous_day,
      to_day              = current_day,
      from_time_point     = previous_time_point,
      to_time_point       = current_time_point,
      treatment_arm,
      arm_group           = vapply(treatment_arm, map_arm_group, character(1)),
      site,
      emoi_s1             = get_emoi(dplyr::coalesce(previous_sample_id, sample_id1)),
      emoi_s2             = current_post_effective_coi_mean,
      pdet_s1             = 0.9,
      pdet_s2             = 0.9,
      dcifer_r            = relatedness,
      .s1                 = dplyr::coalesce(previous_sample_id, sample_id1),
      .s2                 = current_sample_id
    ) %>%
    arrange(id, to_day, to_sample_id) %>%
    group_by(id) %>%
    mutate(step_no = dplyr::row_number()) %>%
    ungroup() %>%
    mutate(
      # A "substituted baseline" consecutive pair is the first step for a patient
      # whose FROM sample is not itself a true D0 (i.e. their real D0 failed
      # genotyping and an upstream process substituted their first recurrence
      # as the effective baseline). Only ever TRUE for step_no == 1.
      is_substituted_baseline = step_no == 1L &
                                !is.na(from_time_point) &
                                from_time_point != "D0"
    )
}

pair_tbl <- pair_tbl %>%
  mutate(
    .idx_s1 = match(.s1, dfsmp$smp_id),
    .idx_s2 = match(.s2, dfsmp$smp_id)
  )

n_substituted <- sum(pair_tbl$is_substituted_baseline, na.rm = TRUE)
message(sprintf(">> substituted-baseline first intervals detected: %d (exclude_substituted_first=%s)",
                n_substituted, exclude_substituted_first))

# Build set of patient IDs in PQ arms for window prior discounting.
# pq_arms is supplied via the CLI flag (space-separated arm names).
# If not supplied, no discount is applied (all arms treated as non-PQ).
pq_arm_pids <- if (length(pq_arms_flag) > 0) {
  pair_tbl %>%
    filter(treatment_arm %in% pq_arms_flag) %>%
    pull(id) %>% unique()
} else {
  character(0)
}
message(sprintf(">> PQ arm patients: %d (pq_discount=%.2f applied to pL in windows 2-3 of its global schedule)",
                length(pq_arm_pids), pq_discount))
message(sprintf(">> PQ arms specified: %s",
                if (length(pq_arms_flag)) paste(pq_arms_flag, collapse=", ") else "none"))

# ── Per-patient drug-group lookup (for al_arms/cq_arms pC window/prior overrides) ─
# Looks up each patient's treatment_arm from the persistent, full pair_tbl
# (never reassigned across subset runs, so this remains correct regardless of
# which subset — full / first_rec_only / mono_only — is currently active).
# Returns pC-specific parameters ONLY — pL is never looked up per-patient here,
# since it never varies by drug; callers use the global w1/w2/pL_w directly.
get_patient_drug_params <- function(pid) {
  arm <- pair_tbl$treatment_arm[match(pid, pair_tbl$id)]
  if (!is.null(drug_params$al) && arm %in% drug_params$al$arms) return(drug_params$al)
  if (!is.null(drug_params$cq) && arm %in% drug_params$cq$arms) return(drug_params$cq)
  list(w1C = w1, w2C = w2, pC_w = pC_w, max_recrud_days = max_recrud_days)
}

has_geno <- !is.na(pair_tbl$.idx_s1) & !is.na(pair_tbl$.idx_s2)
message(sprintf(">> pairs with genotypes for allele-sharing: %d / %d",
                sum(has_geno), nrow(pair_tbl)))

# ── Window-based prior (primary) ──────────────────────────────────────────────
# pC and pL now follow INDEPENDENT window-day schedules:
#   - pC uses w1C_use/w2C_use (drug-specific if the patient's arm is in
#     al_arms/cq_arms, otherwise the global w1/w2).
#   - pL ALWAYS uses the global w1/w2 — there is no w1L_use/w2L_use argument,
#     deliberately, because pL never varies by drug (see header note).
# Falls back to legacy armtime log-normal if window_days not specified but
# armtime is, for backwards compatibility.
tt_prior_row <- function(dt, is_pq_arm = FALSE,
                          w1C_use = w1, w2C_use = w2,
                          pC_w_use = pC_w) {
  if (!is.finite(dt) || dt < 0) return(c(C = 1/3, L = 1/3, I = 1/3))

  # pC window — may be drug-specific
  winC <- if (dt <= w1C_use) 1L else if (dt <= w2C_use) 2L else 3L
  pC   <- pC_w_use[winC]

  # pL window — ALWAYS the global schedule, identical for every arm
  winL <- if (dt <= w1) 1L else if (dt <= w2) 2L else 3L
  pL   <- pL_w[winL]
  # PQ arms: discount pL in windows 2 and 3 of its (global) schedule
  # (hypnozoite reservoir cleared)
  if (is_pq_arm && winL >= 2L) pL <- pL * pq_discount

  # lambda_foi-derived reinfection probability — uncapped, grows freely with time
  # Reinfection can occur at any time regardless of arm; influenced only by transmission
  pI  <- 1 - exp(-lambda_foi_w * max(dt, 0))

  vec <- c(C = pC, L = pL, I = pI)
  vec <- pmax(vec, 1e-9)
  vec / sum(vec)
}

# ── Legacy log-normal prior (kept for --armtime back-compat) ──────────────────
tt_prior_row_lognormal <- function(dt) {
  if (!use_armtime || !is.finite(dt) || dt < 0) return(c(C = 1/3, L = 1/3, I = 1/3))
  pI   <- 1 - exp(-lambda_foi * max(dt, 0))
  mu   <- log(bump_center) - 0.5 * bump_sdlog^2
  dens <- suppressWarnings(dlnorm(max(dt, 1e-6), meanlog = mu, sdlog = bump_sdlog))
  dens <- dens / max(dens, 1e-12)
  pL   <- pmax(0, pC_base * dens)
  if (is.finite(pI_base)) pI <- max(min(pI, pI_base), 0)
  vec  <- c(C = pC_base, L = pL, I = pI); vec <- pmax(vec, 1e-9); vec / sum(vec)
}

build_pv3rs_inputs <- function(pid, afreq, eps_fl) {
  idx <- which(dfsmp$id == pid)
  if (length(idx) < 2) return(NULL)
  idx  <- idx[order(dfsmp$day[idx], dfsmp$smp_id[idx])]
  sids <- dfsmp$smp_id[idx]
  y_raw <- lapply(sids, function(sid) {
    idx_s <- match(sid, dfsmp$smp_id)
    if (is.na(idx_s)) return(NULL)
    s <- dsmp[[idx_s]]
    if (is.null(s) || !length(s)) return(NULL)
    lapply(s, function(v) {
      nm <- names(v)[v > 0]; if (length(nm)) unique(nm) else NA
    })
  })
  keep_ep <- vapply(y_raw, Negate(is.null), logical(1))
  if (!any(keep_ep)) return(NULL)
  y_raw <- y_raw[keep_ep]
  if (length(y_raw) - 1 > maxrec)
    return(list(y = NULL, fs = NULL, sids = sids,
                note = sprintf("skip: n_rec>%d", maxrec)))
  all_markers <- Reduce(union, lapply(y_raw, names))
  if (!length(all_markers))
    return(list(y = NULL, fs = NULL, sids = sids, note = "no markers"))
  y <- lapply(y_raw, function(ep) {
    o <- setNames(vector("list", length(all_markers)), all_markers)
    for (m in all_markers)
      o[[m]] <- if (!is.null(ep[[m]]) && length(ep[[m]]) > 0) ep[[m]] else NA
    o
  })
  fs <- lapply(all_markers, function(m) {
    alle_epi <- unique(unlist(lapply(y, \(ep) ep[[m]]))); alle_epi <- alle_epi[!is.na(alle_epi)]
    base <- afreq[[m]]
    if (is.null(base)) {
      if (!length(alle_epi)) return(setNames(numeric(0), character(0)))
      vec <- rep(eps_fl, length(alle_epi)); names(vec) <- alle_epi; return(vec / sum(vec))
    }
    miss <- setdiff(alle_epi, names(base))
    if (length(miss)) base <- c(base, setNames(rep(eps_fl, length(miss)), miss))
    base / sum(base)
  })
  names(fs) <- all_markers
  keep_m2 <- names(fs)[vapply(fs, length, integer(1)) > 0]
  if (!length(keep_m2))
    return(list(y = NULL, fs = NULL, sids = sids, note = "empty fs"))
  list(y = lapply(y, \(ep) ep[keep_m2]), fs = fs[keep_m2], sids = sids, note = NULL)
}

run_pv3rs_one <- function(pid, afreq, eps_fl) {
  if (!requireNamespace("Pv3Rs", quietly = TRUE)) return(NULL)
  inp <- build_pv3rs_inputs(pid, afreq, eps_fl); if (is.null(inp)) return(NULL)
  n_rec <- max(0L, length(inp$y) - 1L)
  if (!is.null(inp$note)) {
    if (n_rec < 1L) return(NULL)
    return(tibble(pid = pid, rec_idx = seq_len(n_rec),
                  C = NA_real_, L = NA_real_, I = NA_real_, pv3rs_note = inp$note))
  }
  if (n_rec < 1L) return(NULL)
  idx_ord   <- which(dfsmp$id == pid)
  idx_ord   <- idx_ord[order(dfsmp$day[idx_ord], dfsmp$smp_id[idx_ord])]
  day_diffs <- diff(dfsmp$day[idx_ord])
  is_pq     <- pid %in% pq_arm_pids
  dparams   <- get_patient_drug_params(pid)   # pC-only (w1C, w2C, pC_w, max_recrud_days)
  MOIs_min <- tryCatch(
    Pv3Rs::determine_MOIs(inp$y),
    error = function(e) rep(1L, length(inp$y))
  )
  if (sum(MOIs_min) > 8L) {
    message(sprintf("  >> %s: sum(MOI)=%d>8; running stepwise consecutive pairs",
                    pid, sum(MOIs_min)))
    rows <- vector("list", n_rec)
    for (k in seq_len(n_rec)) {
      ep_pair <- inp$y[c(k, k + 1L)]
      moi_pair <- tryCatch(Pv3Rs::determine_MOIs(ep_pair),
                           error = function(e) c(1L, 1L))
      if (sum(moi_pair) > 8L) {
        rows[[k]] <- tibble(pid = pid, rec_idx = k,
                            C = NA_real_, L = NA_real_, I = NA_real_,
                            pv3rs_note = sprintf("stepwise_sum(MOI)=%d>8", sum(moi_pair)))
        next
      }
      dt_k  <- if (length(day_diffs) >= k) day_diffs[k] else NA_real_
      pr_k  <- matrix(tt_prior_row(dt_k, is_pq_arm = is_pq,
                                    w1C_use = dparams$w1C, w2C_use = dparams$w2C,
                                    pC_w_use = dparams$pC_w),
                      nrow = 1L, ncol = 3L,
                      dimnames = list(NULL, c("C","L","I")))
      res_k <- tryCatch(
        Pv3Rs::compute_posterior(y = ep_pair, fs = inp$fs, prior = pr_k,
                                 MOIs = moi_pair, progress.bar = FALSE),
        error = function(e) NULL
      )
      if (is.null(res_k) || is.null(res_k$marg)) {
        rows[[k]] <- tibble(pid = pid, rec_idx = k,
                            C = NA_real_, L = NA_real_, I = NA_real_,
                            pv3rs_note = "stepwise_compute_failed")
      } else {
        m <- as.data.frame(res_k$marg); colnames(m) <- c("C","L","I")
        rows[[k]] <- tibble(pid = pid, rec_idx = k,
                            C = m[nrow(m), "C"],
                            L = m[nrow(m), "L"],
                            I = m[nrow(m), "I"],
                            pv3rs_note = "stepwise_fallback")
      }
    }
    return(dplyr::bind_rows(rows))
  }
  prior <- matrix(1/3, nrow = n_rec, ncol = 3, dimnames = list(NULL, c("C", "L", "I")))
  if (length(day_diffs) == n_rec) {
    for (k in seq_len(n_rec)) {
      prior[k, ] <- if (use_armtime && !use_window_prior)
        tt_prior_row_lognormal(day_diffs[k])
      else
        tt_prior_row(day_diffs[k], is_pq_arm = is_pq,
                     w1C_use = dparams$w1C, w2C_use = dparams$w2C,
                     pC_w_use = dparams$pC_w)
    }
  } else if (is.finite(pC_base) && is.finite(pI_base)) {
    v <- c(C = pC_base, L = 1 - (pC_base + pI_base), I = pI_base)
    v <- pmax(v, 1e-9); v <- v / sum(v)
    prior[, ] <- rep(v, each = n_rec)
  }
  res <- tryCatch(
    Pv3Rs::compute_posterior(y = inp$y, fs = inp$fs, prior = prior,
                             MOIs = MOIs_min, progress.bar = FALSE),
    error = function(e) NULL
  )
  if (is.null(res) || is.null(res$marg))
    return(tibble(pid = pid, rec_idx = seq_len(n_rec),
                  C = NA_real_, L = NA_real_, I = NA_real_,
                  pv3rs_note = "compute_posterior_failed"))
  df <- as.data.frame(res$marg); colnames(df) <- c("C", "L", "I")
  df$pid <- pid; df$rec_idx <- seq_len(n_rec)
  as_tibble(df)
}

# ── Recrudescence time-based reassignment ─────────────────────────────────────
# P. vivax recrudescence (treatment failure) is biologically implausible beyond
# a drug-specific interval. When Pv3Rs calls C (recrudescence) but the interval
# between the paired samples exceeds the applicable threshold, the call is
# reassigned to the next most probable state (L or I) based on the pv3rs_L and
# pv3rs_I posteriors. The raw pv3rs_C/L/I posteriors are retained unchanged;
# only pv3rs_top and the derived agreement columns are updated.
#
# The threshold is resolved per row from treatment_arm: rows whose arm is in
# al_arms use max_days_al, rows in cq_arms use max_days_cq, everything else
# uses max_days_default (the global max_recrud_days).
reassign_recrud <- function(df, max_days_default,
                            al_arms = character(0), cq_arms = character(0),
                            max_days_al = NULL, max_days_cq = NULL) {
  if (!"pv3rs_top" %in% names(df)) return(df)
  interval <- if ("dr_day" %in% names(df)) {
    df$dr_day
  } else if (all(c("from_day", "to_day") %in% names(df))) {
    df$to_day - df$from_day
  } else {
    rep(NA_real_, nrow(df))
  }

  max_days_vec <- rep(max_days_default, nrow(df))
  if (length(al_arms) && !is.null(max_days_al))
    max_days_vec[df$treatment_arm %in% al_arms] <- max_days_al
  if (length(cq_arms) && !is.null(max_days_cq))
    max_days_vec[df$treatment_arm %in% cq_arms] <- max_days_cq

  reassign_mask <- !is.na(df$pv3rs_top) &
                   df$pv3rs_top == "C" &
                   !is.na(interval) &
                   interval > max_days_vec
  n_reassigned <- sum(reassign_mask)
  if (n_reassigned > 0) {
    message(sprintf(
      ">> Pv3Rs recrudescence reassignment: %d pair(s) reassigned (per-row threshold applied)",
      n_reassigned))
    new_top <- ifelse(
      df$pv3rs_L[reassign_mask] >= df$pv3rs_I[reassign_mask], "L", "I"
    )
    df$pv3rs_top[reassign_mask] <- new_top
    df$pv3rs_reassigned <- FALSE
    df$pv3rs_reassigned[reassign_mask] <- TRUE
  } else {
    df$pv3rs_reassigned <- FALSE
  }
  df
}

# ── make_arm_summary: produce per-arm classification summary TSV ───────────────
#
# Classification is based directly on pv3rs_L using the agreed thresholds:
#   relapse      : pv3rs_L >= 0.75
#   reinfection  : pv3rs_L <= 0.25
#   indeterminate: 0.25 < pv3rs_L < 0.75
#
# This is applied to pv3rs_L after recrudescence reassignment (pv3rs_L will have
# been updated where a C call was reassigned to L or I).
#
# Arguments:
#   df            the classified output tibble (sub_out after Pv3Rs + reassignment
#                 + substituted-baseline filtering)
#   meta_d0       ALL genotyped D0 patients (full population, never subset-filtered).
#                 Used as the denominator for n_d0_genotyped regardless of which
#                 subset is being summarised.
#   pair_mode     "baseline" or "consecutive"
#   table_type    "full" | "first_rec" | "mono"
#   mono_emoi_th  threshold used for mono filter (only used when table_type == "mono")
#
# NOTE: this function groups strictly by whatever is in df$treatment_arm /
# meta_d0$treatment_arm — it has no knowledge of arm_group. The grouped
# summary (summary_by_arm_group*.tsv) is produced by remapping treatment_arm
# to the group label before calling this function again (see run_and_write).
#
# "No recurrences" = patients with a genotyped D0 who have no classified pair.
# NOTE: when exclude_substituted_first=on, a patient whose only interval was the
# suppressed substituted-baseline one will now correctly show as having "no
# recurrence" in this summary, since their (excluded) pair no longer counts.
make_arm_summary <- function(df, meta_d0, pair_mode, table_type, mono_emoi_th = 1.1) {

  rel_thresh  <- 0.75
  rei_thresh  <- 0.25

  # D0 genotyped patients per arm (denominator)
  d0_by_arm <- meta_d0 %>%
    count(treatment_arm, name = "n_d0_genotyped")

  # Patients with at least one classified pair
  classified_pids <- unique(df$id[!is.na(df$pv3rs_L)])

  # For mono table: pairs are already filtered to both-monoclonal before Pv3Rs ran.
  # No additional pair filtering needed here.

  # Classified pairs: rows with a pv3rs_L value
  df_class <- df %>% filter(!is.na(pv3rs_L))

  # Apply P(L) threshold classification
  df_class <- df_class %>%
    mutate(pl_call = case_when(
      pv3rs_L >= rel_thresh ~ "relapse",
      pv3rs_L <= rei_thresh ~ "reinfection",
      TRUE                  ~ "indeterminate"
    ))

  # Per-arm: classified pair counts and call breakdown
  class_by_arm <- df_class %>%
    group_by(treatment_arm) %>%
    summarise(
      classified_pairs  = n(),
      n_relapse         = sum(pl_call == "relapse"),
      n_reinfection     = sum(pl_call == "reinfection"),
      n_indeterminate   = sum(pl_call == "indeterminate"),
      mean_days         = if ("dr_day" %in% names(.)) mean(dr_day, na.rm = TRUE) else NA_real_,
      median_days       = if ("dr_day" %in% names(.)) median(dr_day, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )

  # Patients with ≥1 classified pair (for "genotyped recurrence ≥1" column)
  rec_pids_by_arm <- df %>%
    filter(id %in% classified_pids) %>%
    distinct(id, treatment_arm) %>%
    count(treatment_arm, name = "n_with_recurrence")

  # For mono table: no-pairs column = mono D0 patients with no mono recurrence pair
  if (table_type == "mono") {
    mono_class_pids <- unique(df_class$id)
    no_pairs_by_arm <- meta_d0 %>%
      group_by(treatment_arm) %>%
      summarise(
        n_no_mono_pairs = sum(!id %in% mono_class_pids),
        .groups = "drop"
      )
  }

  # No recurrences = D0 genotyped patients with no classified pair at all
  no_rec_pids_by_arm <- meta_d0 %>%
    group_by(treatment_arm) %>%
    summarise(
      n_no_recurrence = sum(!id %in% classified_pids),
      .groups = "drop"
    )

  # Assemble
  summ <- d0_by_arm %>%
    left_join(rec_pids_by_arm,   by = "treatment_arm") %>%
    left_join(no_rec_pids_by_arm, by = "treatment_arm") %>%
    left_join(class_by_arm,       by = "treatment_arm") %>%
    mutate(
      trial                    = TRIAL_UP,
      pair_mode                = pair_mode,
      table_type               = table_type,
      n_with_recurrence        = tidyr::replace_na(n_with_recurrence, 0L),
      pct_with_recurrence      = round(100 * n_with_recurrence / n_d0_genotyped, 1),
      n_no_recurrence          = tidyr::replace_na(n_no_recurrence, 0L),
      pct_no_recurrence        = round(100 * n_no_recurrence / n_d0_genotyped, 1),
      classified_pairs         = tidyr::replace_na(classified_pairs, 0L),
      n_relapse                = tidyr::replace_na(n_relapse, 0L),
      n_reinfection            = tidyr::replace_na(n_reinfection, 0L),
      n_indeterminate          = tidyr::replace_na(n_indeterminate, 0L),
      pct_relapse              = round(100 * n_relapse      / pmax(classified_pairs, 1L), 1),
      pct_reinfection          = round(100 * n_reinfection  / pmax(classified_pairs, 1L), 1),
      pct_indeterminate        = round(100 * n_indeterminate / pmax(classified_pairs, 1L), 1)
    )

  if (table_type == "mono") {
    summ <- summ %>%
      left_join(no_pairs_by_arm, by = "treatment_arm") %>%
      mutate(n_no_mono_pairs = tidyr::replace_na(n_no_mono_pairs, 0L))
  }

  # Column order
  base_cols <- c("trial", "pair_mode", "table_type", "treatment_arm",
                 "n_d0_genotyped", "n_with_recurrence", "pct_with_recurrence",
                 "n_no_recurrence", "pct_no_recurrence",
                 "classified_pairs")
  if (table_type == "full") {
    base_cols <- c(base_cols, "mean_days", "median_days")
  }
  if (table_type == "mono") {
    base_cols <- c(base_cols, "n_no_mono_pairs")
  }
  class_cols <- c("n_relapse", "pct_relapse",
                  "n_reinfection", "pct_reinfection",
                  "n_indeterminate", "pct_indeterminate")

  summ %>% select(any_of(c(base_cols, class_cols)))
}

# ── run_and_write: run Pv3Rs on a restricted dataset and write outputs ─────────
#
# This function encapsulates the full Pv3Rs run + output writing so it can be
# called once for the full dataset and once per subset analysis (first_rec_only,
# mono_only). Each call operates on an independently filtered dfsmp/dsmp/coi,
# so Pv3Rs sees only the episodes in that subset — this is essential because
# Pv3Rs uses all presented episodes as joint context when computing posteriors.
#
# Arguments:
#   sub_dfsmp   filtered dfsmp (patient/sample index for this run)
#   sub_dsmp    dsmp list subsetted to sub_dfsmp$smp_id
#   sub_coi     COI vector aligned to sub_dfsmp rows
#   sub_pairs   pair_tbl rows for patients in this subset
#   sub_has_geno logical vector aligned to sub_pairs rows
#   sub_shr     pre-computed allele-sharing tibble aligned to sub_pairs rows
#   label       short string for log messages (e.g. "full", "first_rec_only")
#   out_suffix  appended to output filenames: "" for full, "_first_rec_only" etc.
#   subset_note character description written into sanity_check.txt
run_and_write <- function(sub_dfsmp, sub_dsmp, sub_coi,
                          sub_pairs, sub_has_geno, sub_shr,
                          label, out_suffix, subset_note,
                          table_type = "full") {

  message(sprintf("\n>> [%s] n_patients=%d  n_pairs=%d",
                  label, dplyr::n_distinct(sub_dfsmp$id), nrow(sub_pairs)))

  # ── Assemble base output ────────────────────────────────────────────────────
  # IMPORTANT: sub_dfsmp/sub_dsmp/sub_coi (passed to Pv3Rs below) are NOT
  # filtered for substituted baselines — Pv3Rs still sees the full episode
  # sequence, including any substituted first episode, as joint context. Only
  # sub_out (the reporting table written to disk) has substituted-baseline rows
  # suppressed, and only when exclude_substituted_first=on.
  sub_out <- bind_cols(
    sub_pairs %>% select(-.s1, -.s2, -.idx_s1, -.idx_s2),
    sub_shr
  )

  # ── Suppress substituted-baseline rows from reporting (does not affect Pv3Rs) ──
  if (exclude_substituted_first == "on" && "is_substituted_baseline" %in% names(sub_out)) {
    n_excl <- sum(sub_out$is_substituted_baseline, na.rm = TRUE)
    if (n_excl > 0) {
      message(sprintf(">> [%s] excluding %d substituted-baseline row(s) from reporting (exclude_substituted_first=on)",
                      label, n_excl))
      sub_out <- sub_out %>% filter(!is_substituted_baseline)
    }
  }

  # ── Pv3Rs ──────────────────────────────────────────────────────────────────
  # build_pv3rs_inputs and run_pv3rs_one read dfsmp/dsmp/coi/afreq_base from
  # the enclosing environment. We temporarily rebind them to the subset versions
  # so the functions operate on only the episodes in this run.
  # For mono_only: afreq_base and eps_floor are also recomputed from monoclonal
  # D0 samples only — polyclonal samples would distort frequency estimates.
  old_dfsmp <- dfsmp; old_dsmp <- dsmp; old_coi <- coi
  dfsmp <<- sub_dfsmp; dsmp <<- sub_dsmp; coi <<- sub_coi

  # af_run is already set correctly at the top of the script:
  # mono_only on  -> af_mono (monoclonal D0 only)
  # mono_only off -> af_full (all D0)

  if (pv3rs_flag == "on") {
    if (!requireNamespace("Pv3Rs", quietly = TRUE)) {
      message(">> Pv3Rs not installed; skipping.")
    } else {
      message(sprintf(">> [%s] running Pv3Rs", label))
      pids_sub <- unique(sub_dfsmp$id)
      pv_tab   <- dplyr::bind_rows(lapply(pids_sub, run_pv3rs_one,
                                          afreq  = af_run$afreq,
                                          eps_fl = af_run$eps_floor))
      if (nrow(pv_tab)) {
        pv_tab <- pv_tab %>%
          mutate(pv3rs_top = c("C", "L", "I")[max.col(
            as.matrix(select(., C, L, I)), ties.method = "first")])
        rec_col <- if (pair_mode == "baseline") "recurrence_no" else "step_no"
        sub_out <- sub_out %>%
          left_join(
            pv_tab %>% transmute(id = pid, !!rec_col := rec_idx,
                                 pv3rs_C = C, pv3rs_L = L, pv3rs_I = I, pv3rs_top),
            by = c("id", rec_col)
          )
      }
    }
  }

  # Restore global dfsmp/dsmp/coi
  dfsmp <<- old_dfsmp; dsmp <<- old_dsmp; coi <<- old_coi

  # ── Recrudescence reassignment ─────────────────────────────────────────────
  if (pv3rs_flag == "on" && "pv3rs_top" %in% names(sub_out))
    sub_out <- reassign_recrud(
      sub_out, max_recrud_days,
      al_arms     = if (!is.null(drug_params$al)) drug_params$al$arms else character(0),
      cq_arms     = if (!is.null(drug_params$cq)) drug_params$cq$arms else character(0),
      max_days_al = if (!is.null(drug_params$al)) drug_params$al$max_recrud_days else NULL,
      max_days_cq = if (!is.null(drug_params$cq)) drug_params$cq$max_recrud_days else NULL
    )

  # ── Write output TSV ───────────────────────────────────────────────────────
  fname <- paste0("per_infection_classifications", out_suffix, ".tsv")
  fp    <- file.path(out_dir, fname)
  write_tsv(sub_out, fp)
  cat(sprintf(">> [%s] wrote: %s\n", label, fp))

  # ── Sanity summary ─────────────────────────────────────────────────────────
  summ_v <- function(x) {
    if (!length(x) || all(is.na(x))) return(sprintf("n=%d  all-NA", length(x)))
    v <- x[is.finite(x)]
    sprintf("n=%d  NA=%d  min=%.3f  q25=%.3f  med=%.3f  q75=%.3f  max=%.3f",
            length(x), sum(is.na(x)), min(v),
            as.numeric(quantile(v, 0.25)), median(v),
            as.numeric(quantile(v, 0.75)), max(v))
  }
  sc <- c(
    sprintf("TRIAL: %s | pair_mode: %s | subset: %s | site: %s",
            TRIAL_UP, pair_mode, subset_note,
            if (length(site_filter)) paste(site_filter, collapse = ", ") else "all sites"),
    sprintf("n_patients: %d  |  n_pairs: %d  |  pairs_with_geno: %d",
            dplyr::n_distinct(sub_dfsmp$id), nrow(sub_pairs), sum(sub_has_geno)),
    sprintf("r_bg = %.4f  (computed from full-dataset D0 samples, n=%d)", rbg, n_d0),
    sprintf("afreq_base: %s D0 samples used for allele frequency estimation",
            if (table_type == "mono")
              sprintf("monoclonal only (eMOI <= %.2f), n=%d", mono_emoi_th, sum(sub_dfsmp$day == 0, na.rm = TRUE))
            else
              sprintf("all, n=%d", n_d0)),
    sprintf("COI source: MOIRE post_effective_coi_mean (fallback COI=1 for %d samples)",
            n_coi_na),
    sprintf("Pv3Rs=%s | max_recrud_days=%.0f | af_model=%s alpha=%.2f eps=%.3g",
            pv3rs_flag, max_recrud_days, af_model, alpha, eps),
    sprintf("GLOBAL schedule (pL ALWAYS, pC default): w1=%g w2=%g | pC: %s | pL: %s | lambda_foi=%.4f | pq_discount=%.2f",
            w1, w2, paste(pC_w, collapse=" "), paste(pL_w, collapse=" "),
            lambda_foi_w, pq_discount),
    if (!is.null(drug_params$al))
      sprintf("AL-arm pC-ONLY overrides: w1C=%g w2C=%g | pC: %s | max_recrud_days=%g (arms: %s) | pL unaffected, uses global",
              drug_params$al$w1C, drug_params$al$w2C,
              paste(drug_params$al$pC_w, collapse=" "),
              drug_params$al$max_recrud_days, paste(al_arms_flag, collapse=", "))
    else "AL-arm pC overrides: none",
    if (!is.null(drug_params$cq))
      sprintf("CQ-arm pC-ONLY overrides: w1C=%g w2C=%g | pC: %s | max_recrud_days=%g (arms: %s) | pL unaffected, uses global",
              drug_params$cq$w1C, drug_params$cq$w2C,
              paste(drug_params$cq$pC_w, collapse=" "),
              drug_params$cq$max_recrud_days, paste(cq_arms_flag, collapse=", "))
    else "CQ-arm pC overrides: none",
    sprintf("arm_group mapping: %s",
            if (length(arm_group_map)) paste(names(arm_group_map), arm_group_map, sep="->", collapse=", ") else "none"),
    sprintf("exclude_substituted_first: %s (%d substituted-baseline row(s) detected in this subset, %s)",
            exclude_substituted_first,
            sum(sub_pairs$is_substituted_baseline, na.rm = TRUE),
            if (exclude_substituted_first == "on") "excluded from output" else "retained in output"),
    sprintf("input_dir:  %s", input_dir),
    sprintf("geno:       %s", geno_fp),
    sprintf("out_dir:    %s", out_dir),
    sprintf("output_file: %s", fname),
    sprintf("Allele rows dropped (raw genotype cleaning): %d", n_before - n_after),
    paste0("dcifer_r   : ", summ_v(sub_out$dcifer_r)),
    paste0("share_frac : ", summ_v(sub_out$share_frac))
  )
  sc_fname <- paste0("sanity_check", out_suffix, ".txt")
  writeLines(sc, file.path(out_dir, sc_fname))
  cat(sprintf(">> [%s] sanity: %s\n", label, file.path(out_dir, sc_fname)))

  # ── Arm summary TSV ────────────────────────────────────────────────────────
  if (pv3rs_flag == "on" && "pv3rs_L" %in% names(sub_out)) {
    # n_d0_genotyped denominator always comes from ALL genotyped D0 patients,
    # regardless of subset — the subset filter affects which pairs are classified,
    # not who was enrolled and genotyped at baseline.
    meta_d0_full <- meta %>%
      filter(day == 0) %>%
      select(id, treatment_arm) %>%
      distinct(id, .keep_all = TRUE)

    tryCatch({
      summ <- make_arm_summary(
        df           = sub_out,
        meta_d0      = meta_d0_full,
        pair_mode    = pair_mode,
        table_type   = table_type,
        mono_emoi_th = mono_emoi_th
      )
      summ_fname <- paste0("summary_by_arm", out_suffix, ".tsv")
      summ_fp    <- file.path(out_dir, summ_fname)
      write_tsv(summ, summ_fp)
      cat(sprintf(">> [%s] summary: %s\n", label, summ_fp))
    }, error = function(e) {
      message(sprintf(">> [%s] WARNING: arm summary failed: %s", label, conditionMessage(e)))
    })

    # ── Grouped arm summary (reporting only; only written if arm_group supplied) ──
    # Reuses make_arm_summary() unchanged by remapping treatment_arm to the
    # group label before calling it — no effect on classification whatsoever.
    if (length(arm_group_map) > 0) {
      tryCatch({
        meta_d0_grouped <- meta_d0_full %>%
          mutate(treatment_arm = vapply(treatment_arm, map_arm_group, character(1)))
        sub_out_grouped <- sub_out %>%
          mutate(treatment_arm = vapply(treatment_arm, map_arm_group, character(1)))
        summ_grouped <- make_arm_summary(
          df           = sub_out_grouped,
          meta_d0      = meta_d0_grouped,
          pair_mode    = pair_mode,
          table_type   = table_type,
          mono_emoi_th = mono_emoi_th
        )
        summ_g_fname <- paste0("summary_by_arm_group", out_suffix, ".tsv")
        summ_g_fp    <- file.path(out_dir, summ_g_fname)
        write_tsv(summ_grouped, summ_g_fp)
        cat(sprintf(">> [%s] grouped summary: %s\n", label, summ_g_fp))
      }, error = function(e) {
        message(sprintf(">> [%s] WARNING: grouped arm summary failed: %s", label, conditionMessage(e)))
      })
    }
  }

  invisible(sub_out)
}

# ── Helper: subset dfsmp/dsmp/coi/pair_tbl to a given patient ID set ──────────
subset_to_pids <- function(pid_set, src_dfsmp, src_dsmp, src_coi,
                           src_pair_tbl, src_has_geno, src_shr_list) {
  keep_rows <- src_dfsmp$id %in% pid_set
  s_dfsmp   <- src_dfsmp[keep_rows, ]
  s_dsmp    <- src_dsmp[keep_rows]
  s_coi     <- src_coi[keep_rows]

  keep_pairs   <- src_pair_tbl$id %in% pid_set
  s_pairs      <- src_pair_tbl[keep_pairs, ]
  s_has_geno   <- src_has_geno[keep_pairs]
  s_shr        <- bind_rows(src_shr_list[keep_pairs])

  list(dfsmp     = s_dfsmp,
       dsmp      = s_dsmp,
       coi       = s_coi,
       pairs     = s_pairs,
       has_geno  = s_has_geno,
       shr       = s_shr)
}

# ── Pre-compute allele-sharing for all pairs (used by all runs) ────────────────
empty_shr <- tibble(
  shared_loci = NA_integer_, loci_tested = NA_integer_, share_frac = NA_real_
)
shr_list_full <- vector("list", nrow(pair_tbl))
for (k in seq_len(nrow(pair_tbl))) {
  if (has_geno[k]) {
    i   <- pair_tbl$.idx_s1[k]; j <- pair_tbl$.idx_s2[k]
    shr <- allele_share(i, j)
    shr_list_full[[k]] <- tibble(
      shared_loci = shr[["shared"]],
      loci_tested = shr[["tested"]],
      share_frac  = if (shr[["tested"]] > 0L) shr[["shared"]] / shr[["tested"]] else NA_real_
    )
  } else {
    shr_list_full[[k]] <- empty_shr
  }
}

# ── Run 1: Full dataset ────────────────────────────────────────────────────────
message("\n>> ── Full dataset run ──────────────────────────────────────────────")
run_and_write(
  sub_dfsmp    = dfsmp,
  sub_dsmp     = dsmp,
  sub_coi      = coi,
  sub_pairs    = pair_tbl,
  sub_has_geno = has_geno,
  sub_shr      = bind_rows(shr_list_full),
  label        = "full",
  out_suffix   = "",
  subset_note  = "all patients (no subset filter)",
  table_type   = "full"
)

# ── Run 2: First-recurrence-only subset ───────────────────────────────────────
# Re-runs Pv3Rs on only D0 + first recurrence per patient (recurrence_no == 1
# in baseline mode; step_no == 1 in consecutive mode). Pv3Rs sees only this
# two-episode context per patient, so posteriors differ from the full run.
# NOTE: if exclude_substituted_first=on, the (only) interval for an affected
# patient in this subset would be entirely suppressed — such patients simply
# will not appear in per_infection_classifications_first_rec_only.tsv.
if (first_rec_only == "on") {
  message("\n>> ── First-recurrence-only subset run ─────────────────────────────")

  rec_col <- if (pair_mode == "baseline") "recurrence_no" else "step_no"

  if (!rec_col %in% names(pair_tbl)) {
    message(sprintf(">> first_rec_only: column '%s' not found — skipping", rec_col))
  } else {
    # Patients who have a first recurrence in the pair table
    fr_pids <- pair_tbl %>%
      filter(.data[[rec_col]] == 1L) %>%
      pull(id) %>% unique()

    if (!length(fr_pids)) {
      message(">> first_rec_only: no patients with recurrence_no/step_no == 1 — skipping")
    } else {
      # For each patient, keep only D0 and the first recurrence sample.
      # In baseline mode: d0_sample_id and the dr_sample_id where recurrence_no == 1.
      # In consecutive mode: from_sample_id and to_sample_id where step_no == 1.
      fr_pairs <- pair_tbl %>% filter(.data[[rec_col]] == 1L, id %in% fr_pids)

      if (pair_mode == "baseline") {
        fr_sids <- unique(c(fr_pairs$d0_sample_id, fr_pairs$dr_sample_id))
      } else {
        fr_sids <- unique(c(fr_pairs$from_sample_id, fr_pairs$to_sample_id))
      }

      # Restrict dfsmp to only these samples (preserves patient structure)
      fr_dfsmp_rows <- dfsmp$id %in% fr_pids & dfsmp$smp_id %in% fr_sids
      fr_dfsmp      <- dfsmp[fr_dfsmp_rows, ]
      fr_dsmp       <- dsmp[fr_dfsmp_rows]
      fr_coi        <- coi[fr_dfsmp_rows]

      # Restrict pair_tbl and shr_list to first-rec pairs only
      fr_pair_mask  <- pair_tbl$id %in% fr_pids & pair_tbl[[rec_col]] == 1L
      fr_pairs_tbl  <- pair_tbl[fr_pair_mask, ]
      fr_has_geno   <- has_geno[fr_pair_mask]
      fr_shr        <- bind_rows(shr_list_full[fr_pair_mask])

      message(sprintf(">> first_rec_only: %d patients, %d pairs",
                      length(fr_pids), nrow(fr_pairs_tbl)))

      run_and_write(
        sub_dfsmp    = fr_dfsmp,
        sub_dsmp     = fr_dsmp,
        sub_coi      = fr_coi,
        sub_pairs    = fr_pairs_tbl,
        sub_has_geno = fr_has_geno,
        sub_shr      = fr_shr,
        label        = "first_rec_only",
        out_suffix   = "_first_rec_only",
        subset_note  = sprintf(
          "first recurrence only (%s == 1); Pv3Rs run on D0 + first DR per patient",
          rec_col),
        table_type   = "first_rec"
      )
    }
  }
}

# ── Run 3: Monoclonal-D0 subset ───────────────────────────────────────────────
# Re-runs Pv3Rs restricted to patients whose D0 sample has raw eMOI <= mono_emoi_th.
# Default threshold is 1.1 (not 1.0) to account for MOIRE rounding — values just
# above 1.0 may represent truly monoclonal infections.
# The summary TSV additionally filters classified pairs to those where BOTH samples
# have eMOI <= mono_emoi_th, matching the Table Xc definition.
if (mono_only == "on") {
  message("\n>> ── Monoclonal subset run ─────────────────────────────────────────")
  message(sprintf(">> mono_only: eMOI threshold = %.2f (both samples in a pair must qualify)",
                  mono_emoi_th))

  # Filter pair_tbl to pairs where BOTH samples have eMOI <= threshold.
  # emoi_s1 and emoi_s2 carry the raw post_effective_coi_mean for each sample.
  mono_pair_mask <- !is.na(pair_tbl$emoi_s1) & pair_tbl$emoi_s1 <= mono_emoi_th &
                    !is.na(pair_tbl$emoi_s2) & pair_tbl$emoi_s2 <= mono_emoi_th
  mono_pair_tbl  <- pair_tbl[mono_pair_mask, ]

  n_excl_pairs <- nrow(pair_tbl) - nrow(mono_pair_tbl)
  message(sprintf(">> mono_only: %d pairs qualify (both samples eMOI <= %.2f); %d excluded",
                  nrow(mono_pair_tbl), mono_emoi_th, n_excl_pairs))

  # Only patients who have at least one qualifying pair go into Pv3Rs.
  # dfsmp is also restricted to samples that appear in qualifying pairs,
  # so Pv3Rs only sees monoclonal episodes.
  mono_pids <- unique(mono_pair_tbl$id)

  if (!length(mono_pids)) {
    message(">> mono_only: no qualifying pairs — skipping")
  } else {
    # Identify all sample IDs involved in qualifying pairs
    if (pair_mode == "baseline") {
      mono_sids <- unique(c(mono_pair_tbl$d0_sample_id, mono_pair_tbl$dr_sample_id))
    } else {
      mono_sids <- unique(c(mono_pair_tbl$from_sample_id, mono_pair_tbl$to_sample_id))
    }

    # Restrict dfsmp/dsmp/coi to only these samples for these patients
    mono_dfsmp_mask <- dfsmp$id %in% mono_pids & dfsmp$smp_id %in% mono_sids
    mono_dfsmp <- dfsmp[mono_dfsmp_mask, ]
    mono_dsmp  <- dsmp[mono_dfsmp_mask]
    mono_coi   <- coi[mono_dfsmp_mask]

    mono_has_geno <- has_geno[mono_pair_mask]
    mono_shr      <- bind_rows(shr_list_full[mono_pair_mask])

    message(sprintf(">> mono_only: %d patients, %d pairs passed to Pv3Rs",
                    length(mono_pids), nrow(mono_pair_tbl)))

    run_and_write(
      sub_dfsmp    = mono_dfsmp,
      sub_dsmp     = mono_dsmp,
      sub_coi      = mono_coi,
      sub_pairs    = mono_pair_tbl,
      sub_has_geno = mono_has_geno,
      sub_shr      = mono_shr,
      label        = "mono_only",
      out_suffix   = "_mono_only",
      subset_note  = sprintf(
        "both samples monoclonal (eMOI <= %.2f); %d patients, %d pairs; Pv3Rs re-run on subset",
        mono_emoi_th, length(mono_pids), nrow(mono_pair_tbl)),
      table_type   = "mono"
    )
  }
}
