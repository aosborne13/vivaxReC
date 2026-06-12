#!/usr/bin/env Rscript
# =============================================================================
# classify_recurrences.R
#
# Recurrence classification combining DCiFER (pairwise IBD) and Pv3Rs
# (probabilistic joint inference). Allele-sharing statistics are computed
# for all pairs and retained as supplementary descriptive output.
#
# ── Input files ───────────────────────────────────────────────────────────────
# All live under your analysis directory:
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
#   maxmarkers      integer  (Pv3Rs marker cap, default 60)
#   maxrec          integer  (Pv3Rs max recurrences per patient, default 6)
#   max_recrud_days numeric  (post-hoc C->L reassignment threshold, default 45)
#
# ── Three-window prior (primary prior system) ─────────────────────────────────
#   Biologically grounded on P. vivax biology and published radical cure trial data.
#   Defaults are derived from published P. vivax radical cure trial data and should
#   be reviewed and adjusted to reflect the epidemiology of the setting being analysed.
#
#     CQ efficacy example (Getachew et al. 2015): 0.7% early, 8% late treatment failure
#     PQ efficacy example (Taylor et al. 2019 Lancet): ~5.3x reduction in relapse
#
#   window_days     w1 w2  (default: 35 120)
#     Window 1 ≤w1:   recrudescence possible, early relapse possible
#     Window 2 w1-w2: recrudescence very unlikely, primary relapse window
#     Window 3 >w2:   recrudescence implausible, reinfection accumulates
#
#   pC              pC_w1 pC_w2 pC_w3  (default: 0.08 0.01 0.005)
#   pL              pL_w1 pL_w2 pL_w3  (default: 0.25 0.70 0.40)
#     Note: pL values shown are for non-PQ arms. PQ arms receive pL*pq_discount
#     in windows 2 and 3 (window 1 undiscounted — PQ course still being taken)
#
#   lambda_foi_w    daily reinfection rate, uncapped  (default: 0.002/day)
#     pI grows as 1-exp(-lambda_foi_w * t); no ceiling — reinfection possible
#     at any time regardless of treatment arm
#
#   pq_arms         names of PQ treatment arms (space-separated; used for pq_discount)
#                   e.g.: pq_arms "7-day" "14-day"
#                   If omitted, no arms receive the PQ discount.
#
#   armtime         legacy 5-value log-normal prior (back-compat only)
#                   lambda_foi bump_center bump_sdlog pC_base pI_base
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
    cap        = 1,   # kept for back-compat, unused
    maxmarkers = 1,
    maxrec         = 1,
    armtime        = 5,   # legacy: lambda_foi, bump_center, bump_sdlog, pC_base, pI_base
    # ── Three-window prior (replaces armtime log-normal bump) ──────────────────
    # Biologically informed time windows for recurrence state priors.
    # Window 1: 0–w1 days     — recrudescence possible, short-latency relapse possible
    # Window 2: w1+1–w2 days  — recrudescence very unlikely, primary relapse window
    # Window 3: >w2 days      — recrudescence implausible, reinfection accumulates
    window_days    = 2,   # w1 w2  (default: 35 120)
    pC             = 3,   # pC_w1 pC_w2 pC_w3  (default: 0.08 0.01 0.005)
    pL             = 3,   # pL_w1 pL_w2 pL_w3  (default: 0.25 0.70 0.40)
    # pI_cap removed: reinfection risk is uncapped (grows freely with lambda_foi)
    lambda_foi_w   = 1,   # daily force of infection for pI growth (default: 0.002)
    pq_discount    = 1,   # multiplier on pL for PQ arms in windows 2+3 (default: 0.19)
    pq_arms        = -1,  # names of PQ treatment arms (variable length; space-separated)
    first_rec_only = 1,   # on | off — re-run Pv3Rs on first-recurrence patients only
    mono_only      = 1,   # on | off — re-run Pv3Rs on monoclonal-D0 patients only (D0 eMOI <= threshold)
    mono_emoi_th   = 1    # eMOI threshold for monoclonal D0 definition (default: 1.0)
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
pv3rs_flag  <- tolower(cli$pv3rs    %||% "off")
af_model    <- tolower(cli$af_model  %||% "mix")
alpha       <- as.numeric(cli$alpha  %||% 0.15)
eps         <- as.numeric(cli$eps    %||% 0.05)
maxmarkers      <- as.integer(cli$maxmarkers      %||% 60L)
maxrec          <- as.integer(cli$maxrec          %||% 6L)

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

# ── Three-window prior parameters ─────────────────────────────────────────────
# Defaults are grounded in P. vivax biology and published Ethiopian trial data.
# Users should review and adjust these values to reflect their own study setting.
#
# Window 1 (≤w1 days):  recrudescence possible, short-latency relapse possible
# Window 2 (w1-w2 days): recrudescence very unlikely, primary relapse window
# Window 3 (>w2 days):  recrudescence implausible, reinfection accumulates
#
# pq_discount: multiplier on pL for PQ arms in windows 2 & 3
#   default 0.19 = 1/5.3 from published ~5-fold reduction in relapse with PQ
use_window_prior <- !is.null(cli$window_days)
if (use_window_prior) {
  w1           <- as.numeric(cli$window_days[1])  # default 35
  w2           <- as.numeric(cli$window_days[2])  # default 120
  pC_w         <- as.numeric(cli$pC)              # 3 values: w1 w2 w3
  pL_w         <- as.numeric(cli$pL)              # 3 values: w1 w2 w3
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

# Subset analysis flags — each triggers a separate Pv3Rs run on the restricted dataset.
# The full-dataset run always happens first; subsets are additional outputs.
first_rec_only   <- tolower(cli$first_rec_only %||% "off")  # on | off
mono_only        <- tolower(cli$mono_only      %||% "off")  # on | off
mono_emoi_th     <- as.numeric(cli$mono_emoi_th %||% 1.0)  # raw eMOI threshold (not rounded)

stopifnot(first_rec_only %in% c("on", "off"))
stopifnot(mono_only      %in% c("on", "off"))

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

message(sprintf(">> pairs loaded: %d rows (trial=%s, after dropping NA relatedness)",
                nrow(pairs), TRIAL_UP))
if (nrow(pairs) == 0)
  stop("No pairs found for trial ", TRIAL_UP, " in ", pair_mode, " mode.")

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
    id            = as.character(patient_id),   # external_id from patient_infection_sequence
    day           = suppressWarnings(as.integer(day)),
    time_point    = as.character(time_point),
    treatment_arm = as.character(treatment_arm),
    site          = as.character(site)
  )

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

# ── Allele frequencies (π) from D0 samples ───────────────────────────────────
is_d0 <- dfsmp$day == 0
if (!any(is_d0)) stop("No D0 samples found for ", TRIAL_UP, ".")

afreq_base <- dcifer::calcAfreq(dsmp[is_d0], coi[is_d0], tol = 1e-5)
minfr      <- 1 / sum(coi[is_d0])
eps_floor  <- max(minfr * eps, 1e-12)

afreq_base <- lapply(afreq_base, function(af) {
  af[af == 0] <- minfr; af / sum(af)
})

# Augment with every allele observed across all samples (pairwise-stable)
panel_all <- list()
for (s in dsmp) for (L in names(s))
  panel_all[[L]] <- union(panel_all[[L]], unique(names(s[[L]])))

afreq_base <- lapply(names(panel_all), function(L) {
  alleles <- panel_all[[L]]
  v <- afreq_base[[L]]
  if (is.null(v)) {
    v <- rep(eps_floor, length(alleles)); names(v) <- alleles
  } else {
    miss <- setdiff(alleles, names(v))
    if (length(miss)) v <- c(v, setNames(rep(eps_floor, length(miss)), miss))
  }
  v / sum(v)
})
names(afreq_base) <- names(panel_all)

# ── r_bg: trial-specific D0 pairwise dcifer IBD ───────────────────────────────
# Always computed from the current trial's own D0 samples only (no cross-trial mixing).
n_d0 <- sum(is_d0)
message(sprintf(">> computing r_bg from %s D0 samples (n=%d)", TRIAL_UP, n_d0))
rbg <- mean(
  dcifer::ibdDat(dsmp[is_d0], coi[is_d0], afreq_base, pval = FALSE),
  na.rm = TRUE
)
message(sprintf(">> r_bg = %.4f (trial: %s, n_D0 = %d)", rbg, TRIAL_UP, n_d0))

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
  pair_tbl <- pairs %>%
    transmute(
      trial          = TRIAL_UP,
      pair_type      = "D0_to_DR",
      id             = external_id,                  # correct patient identifier
      d0_sample_id   = dplyr::coalesce(baseline_sample_id, sample_id1),
      dr_sample_id   = current_sample_id,
      dr_day         = current_day,
      dr_time_point  = current_time_point,
      treatment_arm,
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
    ungroup()
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
    ungroup()
}

pair_tbl <- pair_tbl %>%
  mutate(
    .idx_s1 = match(.s1, dfsmp$smp_id),
    .idx_s2 = match(.s2, dfsmp$smp_id)
  )

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
message(sprintf(">> PQ arm patients: %d (pq_discount=%.2f applied in windows 2-3)",
                length(pq_arm_pids), pq_discount))
message(sprintf(">> PQ arms specified: %s",
                if (length(pq_arms_flag)) paste(pq_arms_flag, collapse=", ") else "none"))

has_geno <- !is.na(pair_tbl$.idx_s1) & !is.na(pair_tbl$.idx_s2)
message(sprintf(">> pairs with genotypes for allele-sharing: %d / %d",
                sum(has_geno), nrow(pair_tbl)))

# ── Pv3Rs (optional) ─────────────────────────────────────────────────────────
entropy <- function(p) {
  p <- p[p > 0 & is.finite(p)]; if (!length(p)) 0 else -sum(p * log(p))
}
select_markers <- function(fs, k) {
  if (!length(fs)) return(character(0))
  sc <- vapply(fs, \(v) entropy(v / sum(v)), numeric(1))
  head(names(sort(sc, decreasing = TRUE)), k)
}

# ── Window-based prior (primary) ──────────────────────────────────────────────
# Replaces log-normal bump with biologically motivated time windows.
# Falls back to legacy armtime log-normal if window_days not specified but
# armtime is, for backwards compatibility.
tt_prior_row <- function(dt, is_pq_arm = FALSE) {
  if (!is.finite(dt) || dt < 0) return(c(C = 1/3, L = 1/3, I = 1/3))

  # Determine window index
  win <- if (dt <= w1) 1L else if (dt <= w2) 2L else 3L

  pC  <- pC_w[win]
  pL  <- pL_w[win]
  # PQ arms: discount pL in windows 2 and 3 (hypnozoite reservoir cleared)
  if (is_pq_arm && win >= 2L) pL <- pL * pq_discount
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

build_pv3rs_inputs <- function(pid) {
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
  fs_base     <- afreq_base[names(afreq_base) %in% all_markers]
  keep_m      <- select_markers(fs_base, min(maxmarkers, length(fs_base)))
  all_markers <- intersect(all_markers, keep_m)
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
    base <- afreq_base[[m]]
    if (is.null(base)) {
      if (!length(alle_epi)) return(setNames(numeric(0), character(0)))
      vec <- rep(eps_floor, length(alle_epi)); names(vec) <- alle_epi; return(vec / sum(vec))
    }
    miss <- setdiff(alle_epi, names(base))
    if (length(miss)) base <- c(base, setNames(rep(eps_floor, length(miss)), miss))
    base / sum(base)
  })
  names(fs) <- all_markers
  keep_m2 <- names(fs)[vapply(fs, length, integer(1)) > 0]
  if (!length(keep_m2))
    return(list(y = NULL, fs = NULL, sids = sids, note = "empty fs"))
  list(y = lapply(y, \(ep) ep[keep_m2]), fs = fs[keep_m2], sids = sids, note = NULL)
}

run_pv3rs_one <- function(pid) {
  if (!requireNamespace("Pv3Rs", quietly = TRUE)) return(NULL)
  inp <- build_pv3rs_inputs(pid); if (is.null(inp)) return(NULL)
  n_rec <- max(0L, length(inp$y) - 1L)
  if (!is.null(inp$note)) {
    if (n_rec < 1L) return(NULL)
    return(tibble(pid = pid, rec_idx = seq_len(n_rec),
                  C = NA_real_, L = NA_real_, I = NA_real_, pv3rs_note = inp$note))
  }
  if (n_rec < 1L) return(NULL)
  # Infer MOI from allele data using Pv3Rs::determine_MOIs.
  # This is what Pv3Rs was designed for and gives lower, more appropriate
  # MOI estimates than MOIRE eMOI for the sum(MOIs) <= 8 feasibility check.
  # MOIRE eMOI (via the coi vector) is still used for dcifer.
  MOIs_min <- tryCatch(
    Pv3Rs::determine_MOIs(inp$y),
    error = function(e) rep(1L, length(inp$y))
  )
  if (sum(MOIs_min) > 8L) {
    # Stepwise fallback: run each consecutive pair independently.
    # This loses the multi-episode joint inference but recovers classifications
    # for patients who would otherwise be entirely NA.
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
      dt_k     <- if (length(dt) >= k) dt[k] else NA_real_
      is_pq    <- pid %in% pq_arm_pids
      pr_k     <- matrix(tt_prior_row(dt_k, is_pq_arm = is_pq),
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
        # marg has one row per recurrence in the pair; take the recurrence row (row 1)
        rows[[k]] <- tibble(pid = pid, rec_idx = k,
                            C = m[nrow(m), "C"],
                            L = m[nrow(m), "L"],
                            I = m[nrow(m), "I"],
                            pv3rs_note = "stepwise_fallback")
      }
    }
    return(dplyr::bind_rows(rows))
  }
  idx_ord  <- which(dfsmp$id == pid)
  idx_ord  <- idx_ord[order(dfsmp$day[idx_ord], dfsmp$smp_id[idx_ord])]
  dt       <- diff(dfsmp$day[idx_ord])
  is_pq    <- pid %in% pq_arm_pids
  prior    <- matrix(1/3, nrow = n_rec, ncol = 3, dimnames = list(NULL, c("C", "L", "I")))
  if (length(dt) == n_rec) {
    # Window-based prior is always used when available (default)
    # Falls back to legacy log-normal if --armtime supplied without --window_days
    for (k in seq_len(n_rec)) {
      prior[k, ] <- if (use_armtime && !use_window_prior)
        tt_prior_row_lognormal(dt[k])
      else
        tt_prior_row(dt[k], is_pq_arm = is_pq)
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
# ~45 days given strong chloroquine efficacy. When Pv3Rs calls C (recrudescence)
# but the interval between the paired samples exceeds max_recrud_days, the call
# is reassigned to the next most probable state (L or I) based on the pv3rs_L
# and pv3rs_I posteriors. The raw pv3rs_C/L/I posteriors are retained unchanged;
# only pv3rs_top and the derived agreement columns are updated.
reassign_recrud <- function(df, max_days) {
  if (!"pv3rs_top" %in% names(df)) return(df)
  interval <- if ("dr_day" %in% names(df)) {
    df$dr_day
  } else if (all(c("from_day", "to_day") %in% names(df))) {
    df$to_day - df$from_day
  } else {
    rep(NA_real_, nrow(df))
  }
  reassign_mask <- !is.na(df$pv3rs_top) &
                   df$pv3rs_top == "C" &
                   !is.na(interval) &
                   interval > max_days
  n_reassigned <- sum(reassign_mask)
  if (n_reassigned > 0) {
    message(sprintf(
      ">> Pv3Rs recrudescence reassignment: %d pair(s) with interval > %g days reassigned",
      n_reassigned, max_days))
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
                          label, out_suffix, subset_note) {

  message(sprintf("\n>> [%s] n_patients=%d  n_pairs=%d",
                  label, dplyr::n_distinct(sub_dfsmp$id), nrow(sub_pairs)))

  # ── Assemble base output ────────────────────────────────────────────────────
  sub_out <- bind_cols(
    sub_pairs %>% select(-.s1, -.s2, -.idx_s1, -.idx_s2),
    sub_shr
  )

  # ── Pv3Rs ──────────────────────────────────────────────────────────────────
  # build_pv3rs_inputs and run_pv3rs_one read dfsmp/dsmp/coi/afreq_base from
  # the enclosing environment. We temporarily rebind them to the subset versions
  # so the functions operate on only the episodes in this run.
  old_dfsmp <- dfsmp; old_dsmp <- dsmp; old_coi <- coi
  dfsmp <<- sub_dfsmp; dsmp <<- sub_dsmp; coi <<- sub_coi

  if (pv3rs_flag == "on") {
    if (!requireNamespace("Pv3Rs", quietly = TRUE)) {
      message(">> Pv3Rs not installed; skipping.")
    } else {
      message(sprintf(">> [%s] running Pv3Rs", label))
      pids_sub <- unique(sub_dfsmp$id)
      pv_tab   <- dplyr::bind_rows(lapply(pids_sub, run_pv3rs_one))
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
    sub_out <- reassign_recrud(sub_out, max_recrud_days)

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
    sprintf("TRIAL: %s | pair_mode: %s | subset: %s", TRIAL_UP, pair_mode, subset_note),
    sprintf("n_patients: %d  |  n_pairs: %d  |  pairs_with_geno: %d",
            dplyr::n_distinct(sub_dfsmp$id), nrow(sub_pairs), sum(sub_has_geno)),
    sprintf("r_bg = %.4f  (computed from full-dataset D0 samples, n=%d)", rbg, n_d0),
    sprintf("COI source: MOIRE post_effective_coi_mean (fallback COI=1 for %d samples)",
            n_coi_na),
    sprintf("Pv3Rs=%s | max_recrud_days=%.0f | af_model=%s alpha=%.2f eps=%.3g",
            pv3rs_flag, max_recrud_days, af_model, alpha, eps),
    sprintf("window_days: w1=%g w2=%g | pC: %s | pL: %s | lambda_foi=%.4f | pq_discount=%.2f",
            w1, w2, paste(pC_w, collapse=" "), paste(pL_w, collapse=" "),
            lambda_foi_w, pq_discount),
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
  subset_note  = "all patients (no subset filter)"
)

# ── Run 2: First-recurrence-only subset ───────────────────────────────────────
# Re-runs Pv3Rs on only D0 + first recurrence per patient (recurrence_no == 1
# in baseline mode; step_no == 1 in consecutive mode). Pv3Rs sees only this
# two-episode context per patient, so posteriors differ from the full run.
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
          rec_col)
      )
    }
  }
}

# ── Run 3: Monoclonal-D0 subset ───────────────────────────────────────────────
# Re-runs Pv3Rs restricted to patients whose D0 sample has raw eMOI <= mono_emoi_th.
# "Raw" means the unrounded post_effective_coi_mean from MOIRE, not the integer
# used internally by dcifer. Only the D0 sample's eMOI is used for inclusion;
# polyclonal recurrence samples within an included patient are still analysed.
# This matches the sensitivity analysis denominator definition agreed analytically.
if (mono_only == "on") {
  message("\n>> ── Monoclonal-D0 subset run ──────────────────────────────────────")
  message(sprintf(">> mono_only: eMOI threshold = %.2f (raw post_effective_coi_mean)",
                  mono_emoi_th))

  # Identify D0 sample per patient from dfsmp
  d0_rows  <- dfsmp[dfsmp$day == 0, ]
  # Get raw (unrounded) eMOI for D0 samples from the lookup
  d0_emoi  <- get_emoi(d0_rows$smp_id, fallback = NA_real_)

  # Patients whose D0 eMOI is <= threshold AND is not missing
  mono_pids <- d0_rows$id[!is.na(d0_emoi) & d0_emoi <= mono_emoi_th]
  mono_pids <- unique(mono_pids)

  n_excl_missing <- sum(is.na(d0_emoi))
  n_excl_poly    <- sum(!is.na(d0_emoi) & d0_emoi > mono_emoi_th)

  message(sprintf(">> mono_only: %d patients qualify (eMOI <= %.2f)",
                  length(mono_pids), mono_emoi_th))
  message(sprintf(">> mono_only: excluded %d polyclonal D0 (eMOI > %.2f), %d missing eMOI",
                  n_excl_poly, mono_emoi_th, n_excl_missing))

  if (!length(mono_pids)) {
    message(">> mono_only: no qualifying patients — skipping")
  } else {
    sub <- subset_to_pids(mono_pids, dfsmp, dsmp, coi,
                          pair_tbl, has_geno, shr_list_full)

    run_and_write(
      sub_dfsmp    = sub$dfsmp,
      sub_dsmp     = sub$dsmp,
      sub_coi      = sub$coi,
      sub_pairs    = sub$pairs,
      sub_has_geno = sub$has_geno,
      sub_shr      = sub$shr,
      label        = "mono_only",
      out_suffix   = "_mono_only",
      subset_note  = sprintf(
        "monoclonal D0 only (raw eMOI <= %.2f); n_patients=%d; Pv3Rs re-run on subset",
        mono_emoi_th, length(mono_pids))
    )
  }
}