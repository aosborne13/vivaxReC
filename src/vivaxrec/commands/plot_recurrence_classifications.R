#!/usr/bin/env Rscript
# =============================================================================
# plot_recurrence_classifications.R
#
# Plots outputs from classify_recurrences.R. Supports both pair modes;
# reads whichever classification files exist or restricts to the mode(s)
# requested via pair_mode.
#
# ── Expected file layout ──────────────────────────────────────────────────────
#   output-analysis/
#     baseline-to-recurrence/
#       per_infection_classifications.tsv    (pair_mode=baseline)
#     consecutive-pairs/
#       per_infection_classifications.tsv                  (pair_mode=consecutive, all recurrences)
#       per_infection_classifications_first_rec_only.tsv    (pair_mode=consecutive, first_rec_only=on)
#
#   input-files/patient_level_tables/
#     patient_infection_sequence.tsv         (for eMOI plots and individual outcomes)
#
# ── CLI flags ─────────────────────────────────────────────────────────────────
#   trial       your trial name (must match trial column in patient_infection_sequence.tsv)
#   pair_mode   baseline | consecutive | both   (default: both)
#   bl_dir      path to baseline-to-recurrence/ output dir
#   con_dir     path to consecutive-pairs/ output dir
#   bl_file     path to a specific baseline TSV (overrides bl_dir search)
#   con_file    path to a specific consecutive TSV (overrides con_dir search)
#               Use this to point directly at first_rec_only TSVs, e.g.:
#               con_file .../consecutive-pairs/per_infection_classifications_first_rec_only.tsv
#   figs_dir    path for figures output
#   pis         path to patient_infection_sequence.tsv
#   min_rec     minimum recurrences for patient timeline plot (default: 2)
#
#   arm_group   "raw_arm:group_label" pairs, space-separated — reporting-only
#               arm relabelling applied to bl, con, AND patient_infection_sequence.tsv
#               (eMOI-over-time, eMOI-spaghetti, individual-outcomes plots).
#               Does not change any classification, only which facet/legend
#               label a row's treatment_arm renders under.
#               e.g. OPRA:   arm_group "ALPQ:PQ" "CQPQ:PQ"      (AL, CQ stay separate)
#               e.g. IMPROV: arm_group "7-day:PQ" "14-day:PQ"   (Control stays separate)
#
#   first_rec_as_baseline   on | off (default: off)
#               When on, the data loaded via con_dir/con_file is expected to
#               be a first-rec-only CONSECUTIVE-schema file (step_no, to_day,
#               from_day, to_sample_id, from_sample_id — i.e. what
#               classify_recurrences.R writes when run with pair_mode
#               consecutive first_rec_only on). Its columns are renamed to the
#               baseline schema (recurrence_no, dr_day, dr_sample_id) and the
#               result is stored in the same internal slot the script has
#               always used for "one D0-anchored observation per patient"
#               plots. This lets those plots work correctly WITHOUT ever
#               running classify_recurrences.R with pair_mode baseline. With
#               this flag on, con is NOT also populated.
#
# ── Headline outputs ──────────────────────────────────────────────────────────
#   individual_outcomes_stacked_<TRIAL>_<MODE>.png
#     One bar per (grouped) treatment arm; every Day 0-genotyped individual is
#     placed into Relapse / Reinfection / Indeterminate / Recurred-unclassifiable
#     / Did not recur. Generated for WHICHEVER of bl/con is loaded — works for
#     first-rec-only runs (first_rec_as_baseline=on, uses bl, <MODE>=first_rec)
#     and all-recurrences runs (uses con, <MODE>=all_rec) alike.
#   emoi_spaghetti_by_episode_<TRIAL>.png
#     eMOI trajectory per individual, x-axis = infection episode index (0=D0,
#     1=first recurrence, ...) rather than calendar day. Bold line = arm
#     median, shown only where N>=5 contribute. Panel N = individuals who
#     actually recurred (>=1 episode beyond D0), not all Day 0-genotyped.
#     Panels are forced to a single row (nrow=1) with square-ish dimensions.
#
# NOTE: relapse_connected_dots_first_rec_<TRIAL>.png and
# relapse_averted_pv3rs_<TRIAL>.png are DISABLED — the functions remain
# defined below but are not called. Re-enable by uncommenting the two call
# sites near the end of the script (search "DISABLED").
#
# ── Example ───────────────────────────────────────────────────────────────────
#   BASE=/path/to/your/analysis
#
#   # All recurrences (consecutive chain)
#   Rscript plot_recurrence_classifications.R \
#     trial    mytrial pair_mode consecutive \
#     con_dir  $BASE/output-analysis/consecutive-pairs \
#     figs_dir $BASE/output-analysis/figures \
#     pis      $BASE/input-files/patient_level_tables/patient_infection_sequence.tsv \
#     arm_group "ALPQ:PQ" "CQPQ:PQ"
#
#   # First recurrence only (same underlying consecutive schema, aliased to bl)
#   Rscript plot_recurrence_classifications.R \
#     trial    mytrial pair_mode consecutive \
#     con_file $BASE/output-analysis/consecutive-pairs/per_infection_classifications_first_rec_only.tsv \
#     figs_dir $BASE/output-analysis/figures \
#     pis      $BASE/input-files/patient_level_tables/patient_infection_sequence.tsv \
#     arm_group "ALPQ:PQ" "CQPQ:PQ" \
#     first_rec_as_baseline on
# =============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(forcats); library(ggplot2); library(scales); library(grid)
  if (requireNamespace("ggnewscale", quietly = TRUE))
    library(ggnewscale)
})

# ── Helpers ───────────────────────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b
nice       <- function(s) gsub("_", " ", s)
has_cols   <- function(df, cols) all(cols %in% names(df))
n_label <- function(df, id_col = "id") {
  n_pairs <- nrow(df)
  n_indiv <- if (id_col %in% names(df)) dplyr::n_distinct(df[[id_col]]) else NA_integer_
  if (is.na(n_indiv)) sprintf("N = %d pairs", n_pairs)
  else sprintf("N = %d pairs (%d individuals)", n_pairs, n_indiv)
}
safe_quant <- function(x, p) {
  v <- x[is.finite(x) & !is.na(x)]; if (!length(v)) return(NA_real_)
  as.numeric(quantile(v, p))
}
# IMPORTANT: OPRA (and possibly other trials) use purely numeric patient IDs
# (e.g. "202135"). readr auto-detects such columns as double, while IDs
# sourced from patient_infection_sequence.tsv are explicitly cast to
# character. ANY join on "id" between a classification data frame and a
# pis-derived data frame MUST coerce both sides to character first, or dplyr
# will error with "incompatible types". Use this helper everywhere such a
# join happens.
as_id_chr <- function(x) as.character(x)

# ── CLI parsing ───────────────────────────────────────────────────────────────
parse_cli <- function(args) {
  known_flags <- c("trial","pair_mode","bl_dir","con_dir","bl_file","con_file",
                    "figs_dir","pis","min_rec","arm_group","first_rec_as_baseline")
  out <- list(trial = NULL, pair_mode = NULL,
              bl_dir = NULL, con_dir = NULL,
              bl_file = NULL, con_file = NULL,
              figs_dir = NULL, pis = NULL,
              min_rec = NULL,
              arm_group = NULL,
              first_rec_as_baseline = NULL)
  i <- 1L
  while (i <= length(args)) {
    a <- tolower(args[i])
    if (a == "trial"     && i+1 <= length(args)) { out$trial     <- tolower(args[i+1]); i <- i+2; next }
    if (a == "pair_mode" && i+1 <= length(args)) { out$pair_mode <- tolower(args[i+1]); i <- i+2; next }
    if (a == "bl_dir"    && i+1 <= length(args)) { out$bl_dir    <- args[i+1];          i <- i+2; next }
    if (a == "con_dir"   && i+1 <= length(args)) { out$con_dir   <- args[i+1];          i <- i+2; next }
    if (a == "bl_file"   && i+1 <= length(args)) { out$bl_file   <- args[i+1];          i <- i+2; next }
    if (a == "con_file"  && i+1 <= length(args)) { out$con_file  <- args[i+1];          i <- i+2; next }
    if (a == "figs_dir"  && i+1 <= length(args)) { out$figs_dir  <- args[i+1];          i <- i+2; next }
    if (a == "pis"       && i+1 <= length(args)) { out$pis       <- args[i+1];          i <- i+2; next }
    if (a == "min_rec"   && i+1 <= length(args)) { out$min_rec   <- args[i+1];          i <- i+2; next }
    if (a == "first_rec_as_baseline" && i+1 <= length(args)) { out$first_rec_as_baseline <- tolower(args[i+1]); i <- i+2; next }
    if (a == "arm_group") {
      vals <- character(0); j <- i + 1L
      while (j <= length(args) && !tolower(args[[j]]) %in% known_flags) {
        vals <- c(vals, args[[j]]); j <- j + 1L
      }
      out$arm_group <- vals; i <- j; next
    }
    if (i == 1L && !a %in% known_flags)  {
      out$trial <- a; i <- i+1; next
    }
    i <- i+1L
  }
  out
}

cli          <- parse_cli(commandArgs(trailingOnly = TRUE))
trial        <- cli$trial     %||% stop("'trial' argument is required")
trial        <- tolower(trial)
pair_mode    <- cli$pair_mode %||% "both";   stopifnot(pair_mode %in% c("baseline","consecutive","both"))
min_rec      <- as.integer(cli$min_rec %||% 2L)
TRIAL_UP     <- toupper(trial)

first_rec_as_baseline <- cli$first_rec_as_baseline %||% "off"
stopifnot(first_rec_as_baseline %in% c("on","off"))

# ── Arm grouping (reporting only) ───────────────────────────────────────────────
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
  if (is.na(raw_arm)) return(raw_arm)
  if (raw_arm %in% names(arm_group_map)) arm_group_map[[raw_arm]] else raw_arm
}
apply_arm_group <- function(df) {
  if (!length(arm_group_map) || is.null(df) || !"treatment_arm" %in% names(df)) return(df)
  df$treatment_arm <- vapply(as.character(df$treatment_arm), map_arm_group, character(1))
  df
}

bl_dir   <- cli$bl_dir   %||% file.path("output-analysis", "baseline-to-recurrence")
con_dir  <- cli$con_dir  %||% file.path("output-analysis", "consecutive-pairs")
figs_dir <- cli$figs_dir %||% file.path("output-analysis", "figures")
fp_pis   <- cli$pis      %||% file.path("input-files", "patient_level_tables",
                                         "patient_infection_sequence.tsv")

dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)

find_classifications <- function(dir) {
  for (ext in c("tsv", "csv")) {
    fp <- file.path(dir, paste0("per_infection_classifications.", ext))
    if (file.exists(fp)) return(fp)
  }
  file.path(dir, "per_infection_classifications.tsv")
}

fp_bl  <- if (!is.null(cli$bl_file))  cli$bl_file  else find_classifications(bl_dir)
fp_con <- if (!is.null(cli$con_file)) cli$con_file else find_classifications(con_dir)

message(">> trial: ", TRIAL_UP, " | pair_mode: ", pair_mode)
message(">> bl:   ", fp_bl)
message(">> con:  ", fp_con)
message(">> figs_dir: ", figs_dir)
if (first_rec_as_baseline == "on")
  message(">> first_rec_as_baseline: on — con data will be aliased into the baseline-schema slot")

# ── Classification mutators ───────────────────────────────────────────────────
std_tolerant <- function(v) {
  v2 <- forcats::fct_na_value_to_level(
    factor(v, levels = c("relapse","reinfection","indeterminate")),
    level = "indeterminate"
  )
  forcats::fct_drop(v2)
}

DCIFER_BREAKS  <- c(0, 0.2275, 0.4375, 0.95, 1)
DCIFER_LABELS  <- c("Unrelated (\u22640.2275)",
                    "Low relatedness (0.2275\u20130.4375)",
                    "Related (0.4375\u20130.95)",
                    "Highly related (>0.95)")
DCIFER_PAL_BIN <- c(
  "Unrelated (\u22640.2275)"            = "#D9EEF8",
  "Low relatedness (0.2275\u20130.4375)" = "#85C1E9",
  "Related (0.4375\u20130.95)"           = "#2E86C1",
  "Highly related (>0.95)"               = "#0A3D62"
)

dcifer_ibd_bin <- function(r) {
  cut(suppressWarnings(as.numeric(r)),
      breaks         = DCIFER_BREAKS,
      labels         = DCIFER_LABELS,
      include.lowest = TRUE)
}

dcifer_to_class <- function(r, hi = 0.60, lo = 0.20) {
  dplyr::case_when(
    is.na(r) ~ NA_character_,
    r >= hi  ~ "relapse",
    r <= lo  ~ "reinfection",
    TRUE     ~ "indeterminate"
  )
}

pv3rs_map_three <- function(code) dplyr::case_when(
  code == "I" ~ "reinfection",
  code == "C" ~ "recrudescence",
  code == "L" ~ "relapse",
  TRUE        ~ "indeterminate"
)
pv3rs_map_two <- function(code) dplyr::case_when(
  code == "I"            ~ "reinfection",
  code %in% c("C","L")  ~ "relapse",
  TRUE                   ~ NA_character_
)

prep_classifications <- function(df) {
  df %>%
    mutate(
      id           = as_id_chr(id),
      loci_tested  = suppressWarnings(as.numeric(loci_tested)),
      share_frac   = suppressWarnings(as.numeric(share_frac)),
      dcifer_r     = suppressWarnings(as.numeric(dcifer_r)),
      emoi_s1      = suppressWarnings(as.numeric(emoi_s1)),
      emoi_s2      = suppressWarnings(as.numeric(emoi_s2)),
      dcifer_class = dcifer_to_class(dcifer_r),
      dcifer_ibd   = dcifer_ibd_bin(dcifer_r),
      pv3rs_C      = suppressWarnings(as.numeric(pv3rs_C)),
      pv3rs_L      = suppressWarnings(as.numeric(pv3rs_L)),
      pv3rs_I      = suppressWarnings(as.numeric(pv3rs_I)),
      pv3rs_top    = dplyr::if_else(is.na(pv3rs_top) | as.character(pv3rs_top) == "NA",
                                "indeterminate", as.character(pv3rs_top)),
      pv3rs_three  = pv3rs_map_three(pv3rs_top),
      pv3rs        = pv3rs_map_two(pv3rs_top),
      pv3rs_postR  = suppressWarnings(rowSums(cbind(pv3rs_C, pv3rs_L), na.rm = FALSE)),
      pv3rs_call   = dplyr::case_when(
        is.na(pv3rs_L)        ~ "uncertain",
        pv3rs_L >= 0.75       ~ "relapse",
        pv3rs_L <= 0.25       ~ "reinfection",
        TRUE                  ~ "uncertain"
      )
    )
}

load_mode <- function(fp, mode_label) {
  if (!file.exists(fp)) {
    message(">> ", mode_label, " file not found: ", fp)
    return(NULL)
  }
  sep <- if (grepl("\\.csv$", fp)) "," else "\t"
  df <- read_delim(fp, delim = sep, show_col_types = FALSE)
  message(sprintf(">> loaded %s: %d rows", mode_label, nrow(df)))
  prep_classifications(df)
}

alias_consecutive_as_baseline <- function(df) {
  if (is.null(df)) return(NULL)
  need <- c("from_sample_id","to_sample_id","from_day","to_day","step_no")
  missing <- setdiff(need, names(df))
  if (length(missing)) {
    warning("first_rec_as_baseline=on but consecutive data is missing columns: ",
            paste(missing, collapse=", "), " — cannot alias.")
    return(NULL)
  }
  df %>%
    rename(
      d0_sample_id  = from_sample_id,
      dr_sample_id  = to_sample_id,
      dr_day        = to_day,
      recurrence_no = step_no
    ) %>%
    { if ("from_time_point" %in% names(.)) rename(., d0_time_point = from_time_point) else . } %>%
    { if ("to_time_point"   %in% names(.)) rename(., dr_time_point = to_time_point)   else . }
}

if (first_rec_as_baseline == "on") {
  con_raw <- load_mode(fp_con, "consecutive (first-rec, to be aliased)")
  bl  <- alias_consecutive_as_baseline(con_raw)
  con <- NULL
  if (!is.null(bl))
    message(sprintf(">> aliased %d first-rec row(s) from consecutive schema into baseline slot", nrow(bl)))
} else {
  bl  <- if (pair_mode %in% c("baseline",    "both")) load_mode(fp_bl,  "baseline")    else NULL
  con <- if (pair_mode %in% c("consecutive", "both")) load_mode(fp_con, "consecutive") else NULL
}

bl  <- apply_arm_group(bl)
con <- apply_arm_group(con)

if (is.null(bl) && is.null(con))
  stop("No classification files found. Check bl_dir / con_dir / bl_file / con_file.")

is_first_rec <- function(df, fp) {
  if (is.null(df)) return(FALSE)
  if (grepl("first_rec_only", basename(fp), fixed = TRUE)) return(TRUE)
  if ("recurrence_no" %in% names(df)) return(all(df$recurrence_no == 1L, na.rm = TRUE))
  FALSE
}
bl_is_first_rec <- is_first_rec(bl, fp_bl) || (first_rec_as_baseline == "on" && !is.null(bl))
if (bl_is_first_rec)
  message(">> first-rec-only data detected — reinfection rank plot will be generated")

# ── Colour palettes ───────────────────────────────────────────────────────────
PAL_CLASS <- c(
  relapse       = "#6aa9ff",
  recrudescence = "#C0392B",
  reinfection   = "#2fbf71",
  indeterminate = "#f08a8a"
)

METHOD_PAL <- c("Pv3Rs" = "#9B59B6")

MOI_SHAPES <- c(monoclonal = 21L, polyclonal = 24L, unknown = 23L)

LVL_PV3RS_2 <- c("relapse","reinfection")
LVL_PV3RS_3 <- c("reinfection","recrudescence","relapse","indeterminate")
LVL_DCIFER  <- DCIFER_LABELS

get_arm_pal <- function(arm_vec) {
  arms <- sort(unique(arm_vec[!is.na(arm_vec)]))
  okabe <- c("#E69F00","#56B4E9","#009E73","#F0E442",
             "#0072B2","#D55E00","#CC79A7","#999999")
  pal <- setNames(rep(okabe, length.out = length(arms)), arms)
  pal
}

theme_pub <- theme_minimal(base_size = 12) +
  theme(plot.title   = element_text(size = 13, face = "bold"),
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8),
        strip.text   = element_text(face = "bold"))

theme_cd <- theme_pub +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 8),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(0.8, "lines")
  )

hm_theme <- theme_minimal(base_size = 11) +
  theme(plot.title   = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8),
        legend.key.height = unit(3.5,"mm"),
        legend.key.width  = unit(4.0,"mm"),
        axis.text.x  = element_text(size = 9),
        axis.text.y  = element_text(size = 9),
        panel.grid   = element_blank())

NICE_LABELS <- c(
  reinfection     = "Reinfection",
  indeterminate   = "Indeterminate",
  relapse         = "Relapse",
  recrudescence   = "Recrudescence",
  pv3rs           = "Pv3Rs",
  dcifer_class    = "DCiFER"
)
nice_lbl <- function(x) ifelse(x %in% names(NICE_LABELS), NICE_LABELS[x], nice(x))

# ── Heatmap builder ───────────────────────────────────────────────────────────
mk_heat_pair <- function(df, xvar, yvar, x_levels, y_levels,
                          title, file, width = 6.4, height = 4.2) {
  if (is.null(df) || !has_cols(df, c(xvar, yvar))) return(invisible(NULL))
  df2 <- df %>%
    transmute(x = factor(.data[[xvar]], levels = x_levels),
              y = factor(.data[[yvar]], levels = y_levels)) %>%
    tidyr::drop_na(x, y)
  if (!nrow(df2)) return(invisible(NULL))
  cnt  <- df2 %>% count(y, x, name = "n")
  tot  <- cnt %>% group_by(y) %>% summarise(tot = sum(n), .groups = "drop")
  grid <- tidyr::expand_grid(y = factor(y_levels, levels = y_levels),
                              x = factor(x_levels, levels = x_levels))
  tab  <- grid %>%
    left_join(cnt, by = c("y","x")) %>%
    left_join(tot, by = "y") %>%
    mutate(n = replace_na(n, 0L), tot = replace_na(tot, 0L),
           frac = ifelse(tot > 0, n/tot, 0))
  p <- ggplot(tab, aes(x, y, fill = frac)) +
    geom_tile(color = "white") +
    geom_text(aes(label = scales::percent(frac, 1)), color = "white", size = 3) +
    scale_x_discrete(drop = FALSE, labels = nice_lbl) +
    scale_y_discrete(drop = FALSE, labels = nice_lbl) +
    scale_fill_gradient(limits = c(0,1), low = "#0b2239", high = "#59A5FF",
                        breaks = c(0,.25,.5,.75,1), labels = percent,
                        guide = guide_colorbar(barheight = unit(20,"mm"),
                                               barwidth  = unit(4,"mm"),
                                               ticks = FALSE,
                                               title.position = "right",
                                               title = "row share")) +
    labs(x = nice_lbl(xvar), y = nice_lbl(yvar), title = title) + hm_theme
  ggsave(file.path(figs_dir, file), p, width = width, height = height, dpi = 300)
  invisible(p)
}

# ── N-count labels on x axis ─────────────────────────────────────────────────
apply_n_labels <- function(p, df, xvar, levels, facet_arms = FALSE) {
  df2 <- df %>%
    mutate(grp = as.character(.data[[xvar]])) %>%
    filter(!is.na(.data[[xvar]]), as.character(.data[[xvar]]) != "NA")
  if (facet_arms && "treatment_arm" %in% names(df2)) {
    cnt <- df2 %>% count(treatment_arm, grp, name = "n") %>%
      mutate(.x = factor(grp, levels = levels))
    p <- p +
      scale_x_discrete(limits = levels, drop = FALSE) +
      geom_text(data = cnt,
                aes(x = .x, y = -Inf, label = paste0("n=", n),
                    group = interaction(grp, treatment_arm)),
                vjust = -0.3, size = 2.6, color = "grey40",
                inherit.aes = FALSE)
  } else {
    cnt   <- df2 %>% count(grp, name = "n")
    n_for <- function(lv) { v <- cnt$n[cnt$grp == lv]; if (!length(v)) 0L else v }
    labvec <- setNames(
      sapply(levels, \(lv) paste0(lv, "\n(n=", n_for(lv), ")")),
      levels
    )
    p <- p + scale_x_discrete(limits = levels, labels = labvec, drop = FALSE)
  }
  p
}

# ── Violin + boxplot + jitter ─────────────────────────────────────────────────
vbj_plot <- function(df, xvar, yvar, x_levels, x_lab, y_lab, title,
                     facet_arms = FALSE) {
  if (is.null(df) || !has_cols(df, c(xvar, yvar, "dcifer_r"))) return(NULL)
  d <- df %>%
    filter(is.finite(.data[[yvar]]),
           .data[[yvar]] >= 0,
           !is.na(.data[[xvar]])) %>%
    mutate(
      .x         = factor(.data[[xvar]], levels = x_levels),
      dcifer_bin = dcifer_ibd_bin(dcifer_r)
    )
  if (!nrow(d)) return(NULL)
  p <- ggplot(d, aes(x = .x, y = .data[[yvar]], fill = .x)) +
    geom_violin(trim = FALSE, alpha = 0.55, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(aes(color = dcifer_bin), width = 0.15, height = 0,
                alpha = 0.70, size = 1.8) +
    scale_fill_manual(values = PAL_CLASS, guide = "none") +
    scale_color_manual(values = DCIFER_PAL_BIN, name = "Relatedness\n(DCiFER r)",
                       drop = FALSE) +
    coord_cartesian(ylim = c(0, NA), clip = "off") +
    labs(x = x_lab, y = y_lab, title = title,
         subtitle = if (facet_arms && "treatment_arm" %in% names(d)) {
           d %>% group_by(treatment_arm) %>%
             summarise(lbl = sprintf("%s: N=%d pairs (%d indiv)",
                                     first(treatment_arm), dplyr::n(),
                                     dplyr::n_distinct(id)),
                       .groups = "drop") %>%
             pull(lbl) %>% paste(collapse = "  |  ")
         } else {
           n_label(d)
         }) +
    theme_pub + theme(legend.position = "right")
  if (facet_arms && "treatment_arm" %in% names(d))
    p <- p + facet_wrap(~treatment_arm, nrow = 1, scales = "fixed") +
             theme(strip.text = element_text(face = "bold", size = 9))
  apply_n_labels(p, d, ".x", x_levels, facet_arms = facet_arms)
}
vbj_plot_byarm <- function(...) vbj_plot(..., facet_arms = TRUE)

# ── Evidence maps ─────────────────────────────────────────────────────────────
evidence_map_pv3rs <- function(df, mode_label) {
  if (is.null(df) || !("pv3rs_postR" %in% names(df))) return(NULL)

  df <- df %>%
    mutate(
      y_val       = suppressWarnings(as.numeric(dcifer_r)),
      pv3rs_L_val = suppressWarnings(as.numeric(pv3rs_L)),
      pv3rs_C_val = suppressWarnings(as.numeric(pv3rs_C)),
      pv3rs_I_val = suppressWarnings(as.numeric(pv3rs_I))
    )

  if (!any(is.finite(df$y_val))) return(NULL)

  thresh_call_pL <- function(pl) {
    dplyr::case_when(
      !is.finite(pl) ~ "uncertain",
      pl >= 0.75     ~ "relapse",
      pl <= 0.25     ~ "reinfection",
      TRUE           ~ "uncertain"
    )
  }

  lvls <- c("reinfection", "relapse", "uncertain")

  d_classified <- df %>%
    filter(is.finite(y_val), is.finite(pv3rs_L_val)) %>%
    mutate(
      thresh_plot = factor(thresh_call_pL(pv3rs_L_val), levels = lvls),
      x_val       = pv3rs_L_val
    )

  n_failed <- df %>% filter(is.finite(y_val), !is.finite(pv3rs_L_val)) %>% nrow()
  indet_note <- if (n_failed > 0)
    sprintf("  |  %d Pv3Rs failed (excluded)", n_failed)
  else ""

  d <- d_classified
  if (!nrow(d)) return(NULL)

  pv3rs_pal    <- c(reinfection = "#2fbf71",
                    relapse     = "#6aa9ff",
                    uncertain   = "#f08a8a")
  pv3rs_shapes <- c(reinfection = 16L,
                    relapse     = 15L,
                    uncertain   = 4L)

  ggplot(d, aes(x_val, y_val)) +
    annotate("rect", xmin = 0.75, xmax = 1.01, ymin = -Inf, ymax = Inf,
             fill = "#6aa9ff", alpha = 0.07) +
    annotate("rect", xmin = -0.01, xmax = 0.25, ymin = -Inf, ymax = Inf,
             fill = "#2fbf71", alpha = 0.07) +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "#f08a8a",
               linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(shape = thresh_plot, color = thresh_plot),
               alpha = 0.72, size = 3.0) +
    geom_vline(xintercept = 0.75, linetype = 2, color = "grey40") +
    geom_vline(xintercept = 0.25, linetype = 2, color = "grey40") +
    annotate("text", x = 0.875, y = Inf, label = "Relapse",
             hjust = 0.5, vjust = 1.3, size = 2.8, color = "#3a7fd5",
             fontface = "italic") +
    annotate("text", x = 0.125, y = Inf, label = "Reinfection",
             hjust = 0.5, vjust = 1.3, size = 2.8, color = "#1a9e5a",
             fontface = "italic") +
    annotate("text", x = 0.5, y = Inf, label = "Uncertain",
             hjust = 0.5, vjust = 1.3, size = 2.5, color = "#c0392b",
             fontface = "italic") +
    scale_color_manual(values = pv3rs_pal,    name = "Pv3Rs call", drop = TRUE) +
    scale_shape_manual(values = pv3rs_shapes, name = "Pv3Rs call", drop = TRUE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Probability of relapse [P(L)] (Pv3Rs)",
         y = "Pairwise relatedness (DCiFER r)",
         title = sprintf("Evidence map — Pv3Rs (%s)", mode_label),
         subtitle = paste0(n_label(d), indet_note)) +
    theme_pub
}

# ── Relapse probability vs day of recurrence scatter plot ────────────────────
plot_relapse_prob_vs_day <- function(df, day_col, mode_label, x_label) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message(">> patchwork required for relapse_prob_vs_day grid — skipping")
    return(invisible(NULL))
  }
  if (is.null(df)) return(invisible(NULL))
  has_pv  <- all(c("pv3rs_C","pv3rs_L","pv3rs_top") %in% names(df))
  has_dcf <- "dcifer_r" %in% names(df)
  if (!has_pv && !has_dcf) {
    message(">> relapse_prob_vs_day: no classifier columns found — skipping")
    return(invisible(NULL))
  }
  if (!all(c(day_col, "treatment_arm", "id") %in% names(df))) {
    message(">> relapse_prob_vs_day: missing required columns — skipping")
    return(invisible(NULL))
  }
  arm_levels <- sort(unique(na.omit(as.character(df$treatment_arm))))
  pv3rs_call_pal <- c(relapse = "#6aa9ff", reinfection = "#2fbf71",
                       indeterminate = "#f08a8a", "NA" = "grey70")
  dcifer_call_pal <- DCIFER_PAL_BIN

  theme_scatter <- theme_pub +
    theme(strip.text       = element_text(face = "bold", size = 9),
          legend.position  = "right",
          panel.grid.minor = element_blank())

  build_scatter <- function(d_plot, y_col, y_label, call_col, call_levels,
                             call_pal, method_title, indet_at_half = FALSE) {
    if (!all(c(y_col, call_col) %in% names(d_plot))) return(NULL)
    indet_at_half <- FALSE

    d2 <- d_plot %>%
      filter(is.finite(.data[[day_col]])) %>%
      mutate(
        .day      = .data[[day_col]],
        .y_raw    = suppressWarnings(as.numeric(.data[[y_col]])),
        .call_raw = as.character(.data[[call_col]]),
        .uncertain = (!is.na(.call_raw) & .call_raw == "indeterminate"),
        .failed    = !is.finite(.y_raw),
        .y    = .y_raw,
        .call = factor(
          dplyr::case_when(
            .uncertain | .failed  ~ "indeterminate",
            !is.na(.call_raw)     ~ .call_raw,
            TRUE                  ~ "indeterminate"
          ),
          levels = c(call_levels,
                     if ("indeterminate" %in% call_levels) character(0)
                     else "indeterminate")
        ),
        treatment_arm = factor(treatment_arm, levels = arm_levels)
      ) %>%
      filter(is.finite(.y))

    if (!nrow(d2)) return(NULL)

    arm_n <- d2 %>%
      group_by(treatment_arm) %>%
      summarise(n_pairs = dplyr::n(), n_indiv = dplyr::n_distinct(id),
                .groups = "drop") %>%
      mutate(label = sprintf("N=%d (%d indiv)", n_pairs, n_indiv),
             .y = -Inf, .day = Inf)

    p <- ggplot(d2, aes(x = .day, y = .y)) +
      geom_hline(yintercept = 0.5, linetype = "dotted",
                 colour = "#f08a8a", linewidth = 0.5, alpha = 0.7) +
      geom_point(aes(colour = .call,
                     alpha  = ifelse(.uncertain, 0.55, 0.72),
                     shape  = ifelse(.uncertain, 4L, 16L)),
                 size = 2.0) +
      scale_alpha_identity() +
      scale_shape_identity() +
      geom_smooth(data = d2 %>% filter(!.uncertain, !.failed, is.finite(.y)),
                  method = "loess", formula = y ~ x, se = TRUE,
                  colour = "grey30", fill = "grey80",
                  linewidth = 0.7, alpha = 0.22) +
      geom_text(data = arm_n,
                aes(x = .day, y = .y, label = label),
                hjust = 1.05, vjust = -0.4, size = 2.3, colour = "grey45",
                inherit.aes = FALSE) +
      facet_wrap(~treatment_arm, nrow = 1, scales = "fixed") +
      scale_colour_manual(values = call_pal, name = "Call",
                          drop = TRUE, na.value = "#f08a8a") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      labs(x = x_label, y = y_label, title = method_title) +
      theme_scatter
    p
  }

  panels <- list()
  if (has_pv) {
    d_pv <- df %>% mutate(
      pv3rs_L_num = suppressWarnings(as.numeric(pv3rs_L)),
      pv3rs_postR = pv3rs_L_num,
      pv3rs_two   = dplyr::case_when(
        !is.finite(pv3rs_L_num)  ~ "indeterminate",
        pv3rs_L_num >= 0.75      ~ "relapse",
        pv3rs_L_num <= 0.25      ~ "reinfection",
        TRUE                     ~ "indeterminate"
      )
    )
    pv3rs_scatter_pal <- c(relapse       = "#6aa9ff",
                            reinfection   = "#2fbf71",
                            indeterminate = "#f08a8a")
    p_pv <- build_scatter(d_pv, "pv3rs_postR", "P(relapse) [Pv3Rs]",
                           "pv3rs_two",
                           c("relapse", "reinfection", "indeterminate"),
                           pv3rs_scatter_pal,
                           sprintf("Pv3Rs  |  %s %s", TRIAL_UP, mode_label),
                           indet_at_half = TRUE)
    if (!is.null(p_pv)) panels[["pv3rs"]] <- p_pv
  }
  if (has_dcf) {
    p_dcf <- build_scatter(df, "dcifer_r",
                            "Pairwise relatedness r  [DCiFER]",
                            "dcifer_ibd", DCIFER_LABELS, dcifer_call_pal,
                            sprintf("DCiFER r  |  %s %s", TRIAL_UP, mode_label),
                            indet_at_half = FALSE)
    if (!is.null(p_dcf)) panels[["dcifer"]] <- p_dcf
  }
  if (!length(panels)) return(invisible(NULL))

  n_arms <- max(1, length(arm_levels))
  n_meth <- length(panels)
  pw <- n_arms * 2.9 + 1.8
  ph <- n_meth * 3.3 + 1.0

  grid <- patchwork::wrap_plots(panels, ncol = 1) +
    patchwork::plot_annotation(
      title    = sprintf("%s — Relapse probability vs day of recurrence (%s)",
                         TRIAL_UP, mode_label),
      subtitle = "Each row = one classifier  |  Facets = treatment arms  |  Loess smoother \u00b1 95% CI",
      theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                       plot.subtitle = element_text(size = 9, colour = "grey40"))
    )
  fname <- sprintf("relapse_prob_vs_day_grid_%s.png", mode_label)
  ggsave(file.path(figs_dir, fname), grid,
         width = pw, height = ph, dpi = 300, bg = "white", limitsize = FALSE)
  message(">> wrote: ", fname)
  invisible(grid)
}

# ── Patient ribbon helpers ────────────────────────────────────────────────────
build_ribbon_edges <- function(bl, con) {
  rows <- list()
  if (!is.null(bl) && has_cols(bl, c("id","dr_day")))
    rows[["bl"]] <- bl %>%
      transmute(id = as_id_chr(id), from_day = 0, to_day = dr_day,
                pv3rs_three, dcifer_class, pair_type,
                emoi_from = emoi_s1, emoi_to = emoi_s2)
  if (!is.null(con) && has_cols(con, c("id","from_day","to_day")))
    rows[["con"]] <- con %>%
      transmute(id = as_id_chr(id), from_day, to_day,
                pv3rs_three, dcifer_class, pair_type,
                emoi_from = emoi_s1, emoi_to = emoi_s2)
  if (!length(rows)) return(NULL)
  bind_rows(rows) %>%
    filter(is.finite(from_day), is.finite(to_day)) %>%
    arrange(id, from_day, to_day)
}

MOI_CUTOFF <- 1.1

build_ribbon_nodes <- function(edges) {
  bind_rows(
    edges %>% transmute(id = as_id_chr(id), day = from_day, emoi = emoi_from),
    edges %>% transmute(id = as_id_chr(id), day = to_day,   emoi = emoi_to)
  ) %>%
    filter(is.finite(day)) %>%
    distinct(id, day, .keep_all = TRUE) %>%
    mutate(
      emoi        = dplyr::coalesce(suppressWarnings(as.numeric(emoi)), 1.0),
      node_shape  = ifelse(emoi > MOI_CUTOFF, 23L, 21L),
      shape_label = ifelse(emoi > MOI_CUTOFF,
                           paste0("eMOI > ", MOI_CUTOFF, " (polyclonal)"),
                           paste0("eMOI \u2264 ", MOI_CUTOFF, " (single-clone)"))
    )
}

plot_ribbons <- function(edges, color_by, top_n = 30, fname) {
  if (is.null(edges) || !nrow(edges) || !color_by %in% names(edges))
    return(invisible(NULL))

  make_ribbon_plot <- function(edge_sub, title_prefix, n_shown) {
    if (!nrow(edge_sub)) return(NULL)
    edge_sub <- edge_sub %>% mutate(id = as_id_chr(id))
    ord <- edge_sub %>%
      group_by(id) %>%
      summarise(last_day = max(to_day, na.rm = TRUE),
                n_edges  = dplyr::n(), .groups = "drop") %>%
      arrange(desc(last_day), desc(n_edges)) %>%
      mutate(row = dplyr::row_number())
    dat   <- edge_sub %>% inner_join(ord, by = "id") %>% filter(row <= n_shown)
    if (!nrow(dat)) return(NULL)
    nodes <- build_ribbon_nodes(dat) %>% inner_join(ord, by = "id")

    ggplot(dat, aes(x = from_day, y = row)) +
      geom_segment(aes(xend = to_day, yend = row, color = .data[[color_by]]),
                   linewidth = 2, lineend = "round") +
      geom_point(data = nodes,
                 aes(x = day, y = row, shape = shape_label),
                 fill = "white", color = "grey20", size = 2.8, stroke = 0.6) +
      scale_shape_manual(
        values = setNames(c(21L, 23L),
                          c(paste0("eMOI \u2264 ", MOI_CUTOFF, " (single-clone)"),
                            paste0("eMOI > ",    MOI_CUTOFF, " (polyclonal)"))),
        name = "Infection marker"
      ) +
      scale_y_reverse(breaks = dat$row, labels = dat$id,
                      expand = expansion(mult = c(.02, .05))) +
      scale_color_manual(values = PAL_CLASS, na.value = "grey70",
                         name = nice(color_by)) +
      labs(x = "Days since D0", y = NULL,
           title = sprintf("%s — %s", title_prefix, nice(color_by))) +
      theme_pub +
      theme(panel.grid.major.y = element_blank(),
            axis.text.y        = element_text(size = 7))
  }

  ord_all <- edges %>%
    mutate(id = as_id_chr(id)) %>%
    group_by(id) %>%
    summarise(last_day = max(to_day, na.rm = TRUE),
              n_edges  = dplyr::n(), .groups = "drop") %>%
    arrange(desc(last_day), desc(n_edges)) %>%
    mutate(row = dplyr::row_number())
  n_all <- min(top_n, nrow(ord_all))
  p_all <- make_ribbon_plot(edges, sprintf("All arms (top %d)", n_all), top_n)
  if (!is.null(p_all)) {
    ggsave(file.path(figs_dir, fname), p_all,
           width = 10, height = max(4, n_all * 0.32 + 2), dpi = 300,
           limitsize = FALSE)
    message(">> wrote: ", fname)
  }

  arm_col <- if ("treatment_arm" %in% names(edges)) "treatment_arm" else NULL
  if (is.null(arm_col)) return(invisible(p_all))

  arm_levels <- sort(unique(na.omit(as.character(edges[[arm_col]]))))
  arm_plots  <- list()

  for (arm in arm_levels) {
    arm_edges <- edges %>% filter(.data[[arm_col]] == arm)
    n_arm     <- min(top_n, dplyr::n_distinct(arm_edges$id))
    p_arm     <- make_ribbon_plot(arm_edges,
                                  sprintf("%s — %s", arm, TRIAL_UP),
                                  top_n)
    if (is.null(p_arm)) next
    arm_plots[[arm]] <- p_arm
    arm_fname <- sub("\\.png$",
                     sprintf("_%s.png", gsub("[^A-Za-z0-9]", "_", arm)),
                     fname)
    ggsave(file.path(figs_dir, arm_fname), p_arm,
           width = 10, height = max(4, n_arm * 0.32 + 2),
           dpi = 300, limitsize = FALSE)
    message(">> wrote: ", arm_fname)
  }

  if (length(arm_plots) > 1 && requireNamespace("patchwork", quietly = TRUE)) {
    n_cols  <- min(length(arm_plots), 2L)
    n_rows  <- ceiling(length(arm_plots) / n_cols)
    grid_pw <- n_cols * 10
    grid_ph <- n_rows * max(4, min(top_n, 15) * 0.32 + 2)
    combined <- patchwork::wrap_plots(arm_plots, ncol = n_cols) +
      patchwork::plot_annotation(
        title = sprintf("%s — ribbons by arm — %s", TRIAL_UP, nice(color_by)),
        theme = theme(plot.title = element_text(size = 13, face = "bold"))
      )
    grid_fname <- sub("\\.png$", "_by_arm_grid.png", fname)
    ggsave(file.path(figs_dir, grid_fname), combined,
           width = grid_pw, height = grid_ph,
           dpi = 300, bg = "white", limitsize = FALSE)
    message(">> wrote: ", grid_fname)
  }
  invisible(p_all)
}

# ── eMOI over time ────────────────────────────────────────────────────────────
plot_emoi_over_time <- function(fp_pis, trial_up) {
  if (!file.exists(fp_pis)) {
    message(">> patient_infection_sequence.tsv not found at: ", fp_pis,
            " — skipping eMOI trend. Pass pis <path> to enable.")
    return(invisible(NULL))
  }
  pis <- read_tsv(fp_pis, show_col_types = FALSE) %>%
    filter(toupper(trial) == trial_up) %>%
    transmute(
      individual_id = as_id_chr(external_id),
      treatment_arm = as.factor(vapply(as.character(treatment_arm), map_arm_group, character(1))),
      day           = suppressWarnings(as.numeric(day)),
      emoi          = suppressWarnings(as.numeric(post_effective_coi_mean))
    ) %>%
    filter(is.finite(day), is.finite(emoi), !is.na(treatment_arm))
  if (!nrow(pis)) { message(">> no eMOI data for ", trial_up); return(invisible(NULL)) }
  BIN  <- 30L
  summ <- pis %>%
    mutate(day_bin = pmax(0L, floor(day / BIN) * BIN)) %>%
    group_by(treatment_arm, day_bin) %>%
    summarise(med = median(emoi), lo = quantile(emoi, 0.25),
              hi = quantile(emoi, 0.75), n = dplyr::n(), .groups = "drop") %>%
    mutate(day_mid = day_bin + BIN / 2)
  ylim_top <- safe_quant(summ$hi, 0.99)
  ylim_top <- if (is.finite(ylim_top)) ceiling(ylim_top) + 1 else NA_real_
  arm_n_labels <- pis %>%
    group_by(treatment_arm) %>%
    summarise(lbl = sprintf("%s\nN=%d (%d indiv)",
                            first(as.character(treatment_arm)),
                            dplyr::n(),
                            dplyr::n_distinct(individual_id)),
              .groups = "drop")
  arm_lbl_vec <- setNames(arm_n_labels$lbl, arm_n_labels$treatment_arm)
  p <- ggplot(summ, aes(x = day_mid, group = treatment_arm,
                         color = treatment_arm, fill = treatment_arm)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
    geom_line(aes(y = med), linewidth = 1.3) +
    geom_point(aes(y = med), size = 1.6) +
    scale_colour_discrete(labels = arm_lbl_vec) +
    scale_fill_discrete(labels = arm_lbl_vec) +
    labs(x = "Days since D0 (30-day bins)", y = "MOIRE eMOI (posterior mean)",
         color = "Treatment arm", fill = "Treatment arm",
         title = sprintf("eMOI over time by treatment arm — %s", trial_up)) +
    theme_pub +
    theme(plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA))
  if (is.finite(ylim_top))
    p <- p + coord_cartesian(ylim = c(1, ylim_top))
  ggsave(file.path(figs_dir, "emoi_over_time_by_arm.png"), p,
         width = 7.5, height = 4.8, dpi = 300, bg = "white")
  invisible(p)
}

# =============================================================================
# eMOI trajectory by infection episode — "spaghetti" plot
# X-axis = infection episode index (0 = D0, 1 = first recurrence, ...), derived
# purely from chronological rank within each patient (day order). Bold line =
# arm median at each episode, drawn only where N >= min_n_trend individuals
# contribute. Panel N = individuals who actually recurred (>=1 episode beyond
# D0), not all Day 0-genotyped individuals. Panels forced to a single row
# (nrow=1) with dimensions chosen to keep them roughly square rather than
# vertically stretched.
# =============================================================================
plot_emoi_spaghetti <- function(fp_pis, trial_up, min_n_trend = 5L, max_episode = NULL) {
  if (!file.exists(fp_pis)) {
    message(">> emoi_spaghetti: patient_infection_sequence.tsv not found — skipping")
    return(invisible(NULL))
  }
  pis <- read_tsv(fp_pis, show_col_types = FALSE) %>%
    filter(toupper(trial) == trial_up) %>%
    transmute(
      id            = as_id_chr(external_id),
      treatment_arm = vapply(as.character(treatment_arm), map_arm_group, character(1)),
      day           = suppressWarnings(as.numeric(day)),
      emoi          = suppressWarnings(as.numeric(post_effective_coi_mean))
    ) %>%
    filter(is.finite(day), is.finite(emoi))
  if (!nrow(pis)) { message(">> emoi_spaghetti: no data for ", trial_up); return(invisible(NULL)) }

  pis <- pis %>%
    arrange(id, day) %>%
    group_by(id) %>%
    mutate(episode = dplyr::row_number() - 1L) %>%
    ungroup()

  if (!is.null(max_episode)) pis <- pis %>% filter(episode <= max_episode)

  arm_levels <- sort(unique(na.omit(pis$treatment_arm)))
  arm_pal    <- get_arm_pal(arm_levels)
  n_arms_ep  <- length(arm_levels)

  n_recurred <- pis %>%
    group_by(id, treatment_arm) %>%
    summarise(max_ep = max(episode), .groups = "drop") %>%
    filter(max_ep >= 1) %>%
    count(treatment_arm, name = "n_recurred")
  n_lbl_vec <- setNames(n_recurred$n_recurred, n_recurred$treatment_arm)

  trend <- pis %>%
    group_by(treatment_arm, episode) %>%
    summarise(med = median(emoi), n = dplyr::n(), .groups = "drop") %>%
    filter(n >= min_n_trend)

  p <- ggplot() +
    geom_line(data = pis, aes(x = episode, y = emoi, group = id, colour = treatment_arm),
              alpha = 0.12, linewidth = 0.4) +
    geom_line(data = trend, aes(x = episode, y = med, colour = treatment_arm),
              linewidth = 1.4) +
    geom_point(data = trend, aes(x = episode, y = med, colour = treatment_arm),
               size = 2.2) +
    facet_wrap(~treatment_arm, nrow = 1, scales = "free_y",
               labeller = as_labeller(function(x)
                 sprintf("%s  (N=%d recurred)", x, ifelse(x %in% names(n_lbl_vec), n_lbl_vec[x], 0L)))) +
    scale_colour_manual(values = arm_pal, guide = "none") +
    scale_x_continuous(breaks = function(lims) seq(0, floor(lims[2]), 1)) +
    labs(x = "Infection episode (0 = D0)", y = "MOIRE eMOI",
         title = sprintf("%s — eMOI trajectory by infection episode", trial_up),
         subtitle = sprintf("Thin lines = individuals  |  Bold = arm median (shown where N\u2265%d)", min_n_trend)) +
    theme_pub +
    theme(strip.text = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey40"))

  fname <- sprintf("emoi_spaghetti_by_episode_%s.png", trial_up)
  ggsave(file.path(figs_dir, fname), p,
         width = n_arms_ep * 3.6 + 1, height = 4.4, dpi = 300, bg = "white", limitsize = FALSE)
  message(">> wrote: ", fname)
  invisible(p)
}

# =============================================================================
# Individual-level classification outcome — stacked bar (headline figure)
# One bar per treatment arm. Every Day 0-genotyped individual is placed into
# exactly one of five mutually exclusive categories:
#   Relapse / Reinfection / Indeterminate / Recurred-unclassifiable
#   (no valid pair) / Did not recur
#
# "Did not recur" and "Recurred, unclassifiable" are derived purely from
# patient_infection_sequence.tsv (day > 0 presence/absence) and do NOT depend
# on which classification data is supplied — correct whether df is
# first-rec-only (bl, aliased) or all-recurrences (con) data.
#
# df: EITHER bl (first-rec-only, aliased) OR con (all-recurrences). Both
# valid. mode_label distinguishes the two in the output filename/title.
#
# NOTE: OPRA (and possibly other trials) use purely numeric patient IDs,
# which readr auto-detects as double in the classification TSV while
# patient_infection_sequence.tsv's id is character — both sides are
# explicitly coerced with as_id_chr() before the join to avoid a type
# mismatch error.
# =============================================================================
plot_individual_outcomes_stacked <- function(df, fp_pis, trial_up, mode_label = "all_rec") {
  if (is.null(df) || !has_cols(df, c("id","treatment_arm","pv3rs_L"))) {
    message(">> individual_outcomes_stacked (", mode_label, "): no classification data — skipping")
    return(invisible(NULL))
  }
  if (!file.exists(fp_pis)) {
    message(">> individual_outcomes_stacked: patient_infection_sequence.tsv not found — skipping")
    return(invisible(NULL))
  }

  pis <- read_tsv(fp_pis, show_col_types = FALSE) %>%
    filter(toupper(trial) == trial_up) %>%
    transmute(
      id            = as_id_chr(external_id),
      treatment_arm = vapply(as.character(treatment_arm), map_arm_group, character(1)),
      day           = suppressWarnings(as.numeric(day))
    ) %>%
    filter(is.finite(day))
  if (!nrow(pis)) { message(">> individual_outcomes_stacked: no pis rows for ", trial_up); return(invisible(NULL)) }

  d0_pop <- pis %>% filter(day == 0) %>% distinct(id, treatment_arm)
  recurred_ids <- pis %>% filter(day > 0) %>% pull(id) %>% unique()

  patient_status <- df %>%
    mutate(id = as_id_chr(id)) %>%
    filter(id %in% d0_pop$id) %>%
    group_by(id) %>%
    summarise(
      any_relapse       = any(pv3rs_L >= 0.75, na.rm = TRUE),
      any_reinfection   = any(pv3rs_L <= 0.25, na.rm = TRUE),
      any_indeterminate = any(is.finite(pv3rs_L) & pv3rs_L > 0.25 & pv3rs_L < 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  outcomes <- d0_pop %>%
    left_join(patient_status, by = "id") %>%
    mutate(
      any_relapse       = tidyr::replace_na(any_relapse, FALSE),
      any_reinfection   = tidyr::replace_na(any_reinfection, FALSE),
      any_indeterminate = tidyr::replace_na(any_indeterminate, FALSE),
      recurred          = id %in% recurred_ids,
      category = dplyr::case_when(
        any_relapse       ~ "Relapse",
        any_reinfection   ~ "Reinfection",
        any_indeterminate ~ "Indeterminate",
        recurred          ~ "Recurred, unclassifiable\n(no valid pair)",
        TRUE              ~ "Did not recur"
      )
    )

  lvls <- c("Relapse","Reinfection","Indeterminate",
            "Recurred, unclassifiable\n(no valid pair)","Did not recur")
  cat_pal <- setNames(c("#6aa9ff","#2fbf71","#f08a8a","#b39ddb","#d9d9d9"), lvls)

  arm_levels <- sort(unique(outcomes$treatment_arm))
  dat <- outcomes %>%
    mutate(treatment_arm = factor(treatment_arm, levels = arm_levels),
           category = factor(category, levels = rev(lvls))) %>%
    count(treatment_arm, category, name = "n") %>%
    group_by(treatment_arm) %>%
    mutate(pct = 100 * n / sum(n), n_arm = sum(n)) %>%
    ungroup()

  arm_n_lbl   <- dat %>% distinct(treatment_arm, n_arm) %>%
    mutate(lbl = sprintf("%s\n(N=%d)", treatment_arm, n_arm))
  arm_lbl_vec <- setNames(arm_n_lbl$lbl, arm_n_lbl$treatment_arm)

  mode_title <- ifelse(mode_label == "first_rec", "first recurrence only",
                       ifelse(mode_label == "baseline", "D0\u2192DR", "all recurrences"))

  p <- ggplot(dat, aes(x = treatment_arm, y = pct, fill = category)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(pct >= 3, sprintf("%.0f%%\n(n=%d)", pct, n), "")),
              position = position_stack(vjust = 0.5),
              size = 2.9, colour = "grey15", lineheight = 0.9) +
    scale_fill_manual(values = cat_pal, breaks = lvls, name = NULL) +
    scale_x_discrete(labels = arm_lbl_vec) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "% of Day 0-genotyped individuals",
         title = sprintf("%s — individual-level recurrence outcome (%s)", trial_up, mode_title)) +
    theme_pub +
    theme(legend.position = "right", panel.grid.major.x = element_blank())

  fname <- sprintf("individual_outcomes_stacked_%s_%s.png", trial_up, mode_label)
  ggsave(file.path(figs_dir, fname), p,
         width = max(6, length(arm_levels) * 1.8 + 2), height = 6,
         dpi = 300, bg = "white")
  message(">> wrote: ", fname)
  invisible(p)
}

# ── DCiFER distribution plots ─────────────────────────────────────────────────
plot_dcifer_by_arm <- function(df, mode_label, fname) {
  if (is.null(df) || !has_cols(df, c("treatment_arm","dcifer_r","id"))) return(invisible(NULL))
  d <- df %>%
    filter(!is.na(dcifer_r), !is.na(treatment_arm)) %>%
    mutate(treatment_arm = as.factor(treatment_arm))
  if (!nrow(d)) return(invisible(NULL))
  arm_summ <- d %>%
    group_by(treatment_arm) %>%
    summarise(n_pairs = dplyr::n(), n_indiv = dplyr::n_distinct(id),
              mean_r  = mean(dcifer_r, na.rm = TRUE), .groups = "drop")
  sub_lbl <- arm_summ %>%
    mutate(lbl = sprintf("%s: N=%d (%d indiv), mean r=%.2f",
                         treatment_arm, n_pairs, n_indiv, mean_r)) %>%
    pull(lbl) %>% paste(collapse = "\n")
  arm_levels <- levels(d$treatment_arm)
  arm_pal    <- get_arm_pal(arm_levels)
  p <- ggplot(d, aes(x = dcifer_r, fill = treatment_arm, colour = treatment_arm)) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    geom_rug(alpha = 0.25, linewidth = 0.4) +
    geom_vline(data = arm_summ,
               aes(xintercept = mean_r, colour = treatment_arm),
               linetype = "dashed", linewidth = 0.8, alpha = 0.8) +
    facet_wrap(~treatment_arm, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = arm_pal, guide = "none") +
    scale_colour_manual(values = arm_pal, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    geom_vline(xintercept = 0.5, linetype = "dotted",
               colour = "grey50", linewidth = 0.5) +
    labs(x = "Pairwise relatedness (DCiFER r)", y = "Density",
         title = sprintf("DCiFER r distribution by treatment arm — %s (%s)", TRIAL_UP, mode_label),
         subtitle = sub_lbl,
         caption = "Dashed lines = arm mean  |  Dotted line = r=0.5") +
    theme_pub +
    theme(strip.text    = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 7.5, colour = "grey40", lineheight = 1.4),
          plot.caption  = element_text(size = 8, colour = "grey50"))
  ggsave(file.path(figs_dir, fname), p,
         width = 10, height = 6, dpi = 300, bg = "white")
  message(">> wrote: ", fname)
  invisible(p)
}

plot_dcifer_by_arm_overlap <- function(df, mode_label, fname) {
  if (is.null(df) || !has_cols(df, c("treatment_arm","dcifer_r","id"))) return(invisible(NULL))
  d <- df %>%
    filter(!is.na(dcifer_r), !is.na(treatment_arm)) %>%
    mutate(treatment_arm = as.factor(treatment_arm))
  if (!nrow(d)) return(invisible(NULL))
  arm_summ <- d %>%
    group_by(treatment_arm) %>%
    summarise(n_pairs = dplyr::n(), n_indiv = dplyr::n_distinct(id),
              mean_r  = mean(dcifer_r, na.rm = TRUE), .groups = "drop")
  arm_levels <- levels(d$treatment_arm)
  arm_pal    <- get_arm_pal(arm_levels)
  lbl_vec <- arm_summ %>%
    mutate(lbl = sprintf("%s  (N=%d, %d indiv, mean r=%.2f)",
                         treatment_arm, n_pairs, n_indiv, mean_r)) %>%
    { setNames(.$lbl, .$treatment_arm) }
  p <- ggplot(d, aes(x = dcifer_r, fill = treatment_arm, colour = treatment_arm)) +
    geom_density(alpha = 0.25, linewidth = 0.9) +
    geom_vline(data = arm_summ,
               aes(xintercept = mean_r, colour = treatment_arm),
               linetype = "dashed", linewidth = 0.8, alpha = 0.85) +
    geom_vline(xintercept = 0.5, linetype = "dotted",
               colour = "grey50", linewidth = 0.6) +
    scale_fill_manual(values = arm_pal, labels = lbl_vec, name = "Treatment arm") +
    scale_colour_manual(values = arm_pal, labels = lbl_vec, name = "Treatment arm") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Pairwise relatedness (DCiFER r)", y = "Density",
         title = sprintf("DCiFER r distribution by treatment arm — %s (%s)", TRIAL_UP, mode_label),
         caption = "Dashed lines = arm mean  |  Dotted line = r=0.5") +
    theme_pub +
    theme(legend.position = "right",
          legend.text     = element_text(size = 8),
          plot.caption    = element_text(size = 8, colour = "grey50"))
  ggsave(file.path(figs_dir, fname), p,
         width = 9, height = 5, dpi = 300, bg = "white")
  message(">> wrote: ", fname)
  invisible(p)
}

# =============================================================================
# Generate all plots
# =============================================================================

heatmaps_bl  <- list()
heatmaps_con <- list()
if (!is.null(bl)) {}
if (!is.null(con)) {}

ev_bl_pv  <- if (!is.null(bl))  evidence_map_pv3rs(bl,  "baseline")    else NULL
ev_con_pv <- if (!is.null(con)) evidence_map_pv3rs(con, "consecutive") else NULL

for (nm in list(list(ev_bl_pv,  "evidence_map_baseline_pv3rs.png"),
                list(ev_con_pv, "evidence_map_consecutive_pv3rs.png"))) {
  if (!is.null(nm[[1]]))
    ggsave(file.path(figs_dir, nm[[2]]), nm[[1]], width = 6, height = 5, dpi = 300)
}

p_bl_time_pv  <- vbj_plot(bl,"pv3rs_three","dr_day",
  c("reinfection","recrudescence","relapse"),
  "Pv3Rs call (D0\u2192DR)","Days since D0","Time to recurrence \u2014 Pv3Rs (baseline)")
con_dt <- if (!is.null(con) && has_cols(con, c("from_day","to_day")))
  con %>% mutate(dt_days = to_day - from_day) else NULL
p_con_time_pv  <- vbj_plot(con_dt,"pv3rs_three","dt_days",
  c("reinfection","recrudescence","relapse"),
  "Pv3Rs call (DRk\u2192DRk+1)","Days between recurrences",
  "Time between consecutive recurrences \u2014 Pv3Rs")
for (nm in list(
  list(p_bl_time_pv,   "timing_baseline_pv3rs.png"),
  list(p_con_time_pv,  "timing_consecutive_pv3rs.png"))) {
  if (!is.null(nm[[1]]))
    ggsave(file.path(figs_dir, nm[[2]]), nm[[1]], width = 6, height = 5, dpi = 300)
}

p_bl_time_pv_arm   <- vbj_plot_byarm(bl,"pv3rs_three","dr_day",
  c("reinfection","recrudescence","relapse"),
  "Pv3Rs call (D0→DR)","Days since D0",
  "Time to recurrence — Pv3Rs by treatment arm (baseline)")
p_con_time_pv_arm  <- vbj_plot_byarm(con_dt,"pv3rs_three","dt_days",
  c("reinfection","recrudescence","relapse"),
  "Pv3Rs call (DRk→DRk+1)","Days between recurrences",
  "Time between recurrences — Pv3Rs by treatment arm (consecutive)")
for (nm in list(
  list(p_bl_time_pv_arm,  "timing_baseline_pv3rs_by_arm.png",    14, 5),
  list(p_con_time_pv_arm, "timing_consecutive_pv3rs_by_arm.png", 14, 5))) {
  if (!is.null(nm[[1]]))
    ggsave(file.path(figs_dir, nm[[2]]), nm[[1]], width = nm[[3]], height = nm[[4]], dpi = 300)
}

ribbon_edges <- build_ribbon_edges(bl, con)

plot_dcifer_by_arm(bl,  "baseline",    "dcifer_r_by_arm_baseline.png")
plot_dcifer_by_arm(con, "consecutive", "dcifer_r_by_arm_consecutive.png")
plot_dcifer_by_arm_overlap(bl,  "baseline",    "dcifer_r_by_arm_baseline_overlap.png")
plot_dcifer_by_arm_overlap(con, "consecutive", "dcifer_r_by_arm_consecutive_overlap.png")

plot_ibd_combined <- function(df, mode_label, fname,
                               day_col    = NULL,
                               day_split  = FALSE,
                               day_cut    = 120) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message(">> patchwork required for plot_ibd_combined — skipping ", fname)
    return(invisible(NULL))
  }
  if (is.null(df) || !has_cols(df, c("treatment_arm","dcifer_r","id"))) return(invisible(NULL))

  local_arm_levels <- if (exists("arm_levels") && length(arm_levels) > 0)
    arm_levels
  else
    sort(unique(df$treatment_arm[!is.na(df$treatment_arm)]))

  d_base <- df %>%
    filter(!is.na(dcifer_r), !is.na(treatment_arm)) %>%
    mutate(treatment_arm = factor(treatment_arm, levels = local_arm_levels))
  if (!nrow(d_base)) return(invisible(NULL))

  arm_pal <- get_arm_pal(local_arm_levels)

  make_box <- function(d, title_str) {
    summ <- d %>% group_by(treatment_arm) %>%
      summarise(n_pairs = dplyr::n(), n_indiv = dplyr::n_distinct(id),
                .groups = "drop") %>%
      mutate(label = sprintf("N=%d\n(%d indiv)", n_pairs, n_indiv),
             y = Inf)
    ggplot(d, aes(x = treatment_arm, y = dcifer_r,
                  fill = treatment_arm, colour = treatment_arm)) +
      geom_boxplot(width = 0.35, outlier.shape = NA, alpha = 0.6, colour = "grey40") +
      geom_jitter(width = 0.18, height = 0, alpha = 0.65, size = 1.6) +
      geom_text(data = summ, aes(x = treatment_arm, y = y, label = label),
                vjust = 1.4, size = 2.5, colour = "grey35", inherit.aes = FALSE) +
      scale_fill_manual(values = arm_pal, guide = "none") +
      scale_colour_manual(values = arm_pal, guide = "none") +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                         expand = expansion(mult = c(0.02, 0.12))) +
      labs(x = "Treatment arm", y = "Pairwise relatedness (DCiFER r)",
           title = title_str) +
      theme_pub +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }

  make_density <- function(d, title_str, show_legend = FALSE) {
    summ <- d %>% group_by(treatment_arm) %>%
      summarise(n_pairs = dplyr::n(), n_indiv = dplyr::n_distinct(id),
                mean_r  = mean(dcifer_r, na.rm = TRUE), .groups = "drop")
    lbl_vec <- setNames(as.character(summ$treatment_arm), summ$treatment_arm)
    p <- ggplot(d, aes(x = dcifer_r, fill = treatment_arm, colour = treatment_arm)) +
      geom_density(alpha = 0.25, linewidth = 0.9) +
      geom_vline(data = summ,
                 aes(xintercept = mean_r, colour = treatment_arm),
                 linetype = "dashed", linewidth = 0.8, alpha = 0.85) +
      geom_text(data = summ,
                aes(x = mean_r, label = sprintf("%.2f", mean_r), colour = treatment_arm),
                y = Inf, vjust = -0.4, hjust = 0.5, size = 2.5,
                fontface = "bold", inherit.aes = FALSE,
                show.legend = FALSE) +
      geom_vline(xintercept = 0.5, linetype = "dotted",
                 colour = "grey50", linewidth = 0.6) +
      scale_fill_manual(values = arm_pal, labels = lbl_vec, name = "Treatment arm") +
      scale_colour_manual(values = arm_pal, labels = lbl_vec, name = "Treatment arm") +
      scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      coord_cartesian(clip = "off") +
      labs(x = "Pairwise relatedness (DCiFER r)", y = "Density",
           title = title_str,
           caption = "Dashed lines = arm mean  |  Dotted line = r=0.5") +
      theme_pub +
      theme(plot.caption  = element_text(size = 7.5, colour = "grey50"),
            plot.margin   = margin(t = 18, r = 5, b = 5, l = 5))
    if (show_legend)
      p <- p + theme(legend.position = "right",
                     legend.text     = element_text(size = 8))
    else
      p <- p + theme(legend.position = "none")
    p
  }

  all_title <- sprintf("%s %s — all recurrences", TRIAL_UP, mode_label)

  if (!day_split || is.null(day_col) || !day_col %in% names(d_base)) {
    pA <- make_box(d_base, all_title)
    pB <- make_density(d_base, all_title, show_legend = TRUE)
    p  <- patchwork::wrap_plots(pA, pB, ncol = 1, heights = c(1, 1))
    w  <- 8; h <- 10
  } else {
    d_lo <- d_base %>% filter(.data[[day_col]] <= day_cut)
    d_hi <- d_base %>% filter(.data[[day_col]] >  day_cut)

    pA1 <- make_box(d_base, sprintf("All recurrences"))
    pA2 <- make_box(d_lo,   sprintf("\u2264 Day %d", day_cut))
    pA3 <- make_box(d_hi,   sprintf("> Day %d",  day_cut))

    pB1 <- make_density(d_base, "All recurrences",              show_legend = FALSE)
    pB2 <- make_density(d_lo,   sprintf("\u2264 Day %d", day_cut), show_legend = FALSE)
    pB3 <- make_density(d_hi,   sprintf("> Day %d",  day_cut),  show_legend = TRUE)

    p <- patchwork::wrap_plots(
      pA1, pA2, pA3,
      pB1, pB2, pB3,
      ncol = 3
    ) + patchwork::plot_annotation(
      title    = sprintf("%s %s — IBD by treatment arm", TRIAL_UP, mode_label),
      subtitle = sprintf("Left: all pairs  |  Centre: recurrences \u2264 Day %d  |  Right: recurrences > Day %d",
                         day_cut, day_cut),
      theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                       plot.subtitle = element_text(size = 9,  colour = "grey40"))
    )
    w  <- 18; h <- 10
  }

  ggsave(file.path(figs_dir, fname), p,
         width = w, height = h, dpi = 300, bg = "white")
  message(">> wrote: ", fname)
  invisible(p)
}

con_day_col <- if (!is.null(con) && "to_day" %in% names(con)) "to_day" else NULL
bl_day_col  <- if (!is.null(bl)  && "dr_day" %in% names(bl))  "dr_day" else NULL

plot_ibd_combined(con, "consecutive", "ibd_by_arm_consecutive_combined.png",
                  day_col = con_day_col, day_split = FALSE)
plot_ibd_combined(bl,  "baseline",    "ibd_by_arm_baseline_combined.png",
                  day_col = bl_day_col,  day_split = FALSE)

plot_ibd_combined(con, "consecutive", "ibd_by_arm_consecutive_day120_split.png",
                  day_col = con_day_col, day_split = TRUE, day_cut = 120)
plot_ibd_combined(bl,  "baseline",    "ibd_by_arm_baseline_day120_split.png",
                  day_col = bl_day_col,  day_split = TRUE, day_cut = 120)

con_dt_day <- if (!is.null(con) && has_cols(con, c("from_day","to_day")))
  con %>% mutate(dt_days = to_day - from_day) else con
plot_relapse_prob_vs_day(bl,         "dr_day",  "baseline",    "Days since D0")
plot_relapse_prob_vs_day(con_dt_day, "to_day",  "consecutive", "Days since D0")

# ── 5. eMOI plots + individual outcomes ───────────────────────────────────────
plot_emoi_over_time(fp_pis, TRIAL_UP)
plot_emoi_spaghetti(fp_pis, TRIAL_UP)

outcome_df   <- if (!is.null(con)) con else bl
outcome_mode <- if (!is.null(con)) "all_rec" else if (bl_is_first_rec) "first_rec" else "baseline"
plot_individual_outcomes_stacked(outcome_df, fp_pis, TRIAL_UP, outcome_mode)

# ── 6. Multi-panel summaries ──────────────────────────────────────────────────
if (requireNamespace("patchwork", quietly = TRUE)) {
  suppressPackageStartupMessages(library(patchwork))

  all_hm <- Filter(Negate(is.null), c(heatmaps_bl, heatmaps_con))
  if (length(all_hm) > 0) {
    hm_ncol <- min(length(all_hm), 2L)
    hm_nrow <- ceiling(length(all_hm) / hm_ncol)
    p_hm_all <- patchwork::wrap_plots(all_hm, ncol = hm_ncol, byrow = TRUE) +
      patchwork::plot_annotation(
        title = sprintf("%s — method agreement heatmaps", TRIAL_UP),
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      )
    ggsave(file.path(figs_dir, "panel_heatmaps_all.png"), p_hm_all,
           width = hm_ncol * 5.5, height = max(4, 4.2 * hm_nrow), dpi = 300)
  }

}

if (!is.null(p_bl_time_pv) && !is.null(p_con_time_pv)) {
  ggsave(file.path(figs_dir, "panel_timing_comparison_pv3rs.png"),
         (p_bl_time_pv | p_con_time_pv) +
           patchwork::plot_annotation(
             title = sprintf("%s — Pv3Rs timing: baseline vs consecutive", TRIAL_UP),
             theme = theme(plot.title = element_text(size = 13, face = "bold"))),
         width = 12, height = 5.5, dpi = 300)
}

# ── 7. Summary figures ────────────────────────────────────────────────────────
tryCatch({
  all_df <- dplyr::bind_rows(
    if (!is.null(bl))  bl  %>% mutate(pair_mode = "Baseline (D0\u2192DR)"),
    if (!is.null(con)) con %>% mutate(pair_mode = "Consecutive (DR\u2192DR)")
  )

  class_by_arm_plot <- function(df_in, class_col, lvls, method_label, panel_lbl) {
    if (is.null(df_in) ||
        !has_cols(df_in, c("treatment_arm", class_col, "pair_mode"))) return(NULL)
    dat <- df_in %>%
      filter(!is.na(.data[[class_col]]), !is.na(treatment_arm)) %>%
      mutate(.cls = factor(.data[[class_col]], levels = lvls)) %>%
      count(treatment_arm, pair_mode, .cls) %>%
      group_by(treatment_arm, pair_mode) %>%
      mutate(prop = n / sum(n)) %>% ungroup()
    sub_lbl <- df_in %>%
      filter(!is.na(.data[[class_col]]), !is.na(treatment_arm)) %>%
      group_by(pair_mode) %>%
      summarise(lbl = sprintf("%s: N=%d pairs (%d indiv)",
                              first(pair_mode), dplyr::n(),
                              dplyr::n_distinct(id)),
                .groups = "drop") %>%
      pull(lbl) %>% paste(collapse = "  |  ")
    ggplot(dat, aes(treatment_arm, prop, fill = .cls)) +
      geom_col(width = 0.7, color = "white", linewidth = 0.3) +
      geom_text(aes(label = ifelse(prop >= 0.04,
                                   paste0(scales::percent(prop, 1), "\n(n=", n, ")"), "")),
                position = position_stack(vjust = 0.5),
                size = 2.5, color = "white", fontface = "bold") +
      facet_wrap(~pair_mode) +
      scale_fill_manual(values = PAL_CLASS, name = "Classification", drop = FALSE) +
      scale_y_continuous(labels = percent) +
      labs(x = "Treatment arm", y = "Proportion of recurrences",
           title = sprintf("%s  Classification by treatment arm (%s)", panel_lbl, method_label),
           subtitle = sub_lbl) +
      theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }

  pA_pv3rs <- class_by_arm_plot(all_df,
                                 "pv3rs_call",
                                 c("reinfection", "uncertain", "relapse"),
                                 "Pv3Rs", "A")
  if (!is.null(pA_pv3rs))
    ggsave(file.path(figs_dir, "panel_A_classification_pv3rs.png"),
           pA_pv3rs, width = 10, height = 5.5, dpi = 300)
  message(">> wrote per-method Panel A classification figures")

}, error = function(e) {
  message(">> WARNING: summary figures failed: ", conditionMessage(e))
})

build_cd_points <- function(df, day_col) {
  if (is.null(df)) return(NULL)
  need <- c("id", day_col, "treatment_arm", "pv3rs_L", "pv3rs_call")
  if (!all(need %in% names(df))) return(NULL)

  emoi_col   <- if ("emoi_s2"  %in% names(df)) "emoi_s2"  else NULL
  dcifer_col <- if ("dcifer_r" %in% names(df)) "dcifer_r" else NULL

  df %>%
    transmute(
      id            = as_id_chr(id),
      day           = .data[[day_col]],
      treatment_arm = treatment_arm,
      value         = suppressWarnings(as.numeric(pv3rs_L)),
      uncertain     = is.na(pv3rs_L) | pv3rs_call == "uncertain",
      dcifer_r      = if (!is.null(dcifer_col))
                        suppressWarnings(as.numeric(.data[[dcifer_col]]))
                      else NA_real_,
      shape_lbl     = if (!is.null(emoi_col))
                        dplyr::case_when(
                          is.na(.data[[emoi_col]])  ~ "unknown",
                          .data[[emoi_col]] <= 1.1  ~ "monoclonal",
                          TRUE                      ~ "polyclonal"
                        )
                      else "unknown"
    )
}

# NOTE: plot_connected_dots_first_rec() remains defined but is not called
# below (see "DISABLED" near the end of the script). Left intact in case you
# want to re-enable it.
plot_connected_dots_first_rec <- function(bl) {
  if (is.null(bl)) return(invisible(NULL))
  rec_col <- if ("recurrence_no" %in% names(bl)) "recurrence_no" else NULL
  df1 <- if (!is.null(rec_col)) bl %>% filter(.data[[rec_col]] == 1L) else bl
  if (!nrow(df1)) { message(">> connected_dots_first_rec: no recurrence_no==1 rows — skipping"); return(invisible(NULL)) }

  pts <- build_cd_points(df1, day_col = "dr_day")
  if (is.null(pts) || !nrow(pts)) return(invisible(NULL))
  arm_levels <- sort(unique(na.omit(pts$treatment_arm)))
  n_arms     <- length(arm_levels)

  connectors <- pts %>%
    group_by(id, day, treatment_arm) %>%
    summarise(ymin = min(value), ymax = max(value),
              any_uncertain = any(uncertain), .groups = "drop")

  arm_plots <- lapply(seq_along(arm_levels), function(i) {
    arm      <- arm_levels[i]
    show_leg <- i == n_arms
    pd  <- pts        %>% filter(treatment_arm == arm)
    cd  <- connectors %>% filter(treatment_arm == arm)
    n_p <- dplyr::n_distinct(pd$id)
    p <- ggplot() +
      geom_segment(data = cd,
                   aes(x = day, xend = day, y = ymin, yend = ymax,
                       linetype = any_uncertain),
                   colour = "grey60", linewidth = 0.7) +
      scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dashed"),
                            labels = c("FALSE" = "Certain", "TRUE" = "Uncertain"),
                            name   = NULL,
                            guide  = if (show_leg) guide_legend() else "none") +
      geom_point(data = pd,
                 aes(x = day, y = value, fill = dcifer_r,
                     alpha = uncertain, shape = shape_lbl),
                 colour = "grey30", size = 2.5, stroke = 0.3) +
      scale_fill_gradientn(
        colours  = c("#deebf7","#9ecae1","#3182bd","#08519c"),
        limits   = c(0, 1),
        name     = "DCiFER r",
        na.value = "grey70",
        guide    = if (show_leg) guide_colourbar(barwidth = 0.8, barheight = 5) else "none"
      ) +
      scale_alpha_manual(values = c("FALSE" = 0.85, "TRUE" = 0.35), guide = "none") +
      scale_shape_manual(values = MOI_SHAPES, name = "eMOI",
                         guide = if (show_leg) guide_legend() else "none") +
      scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
      scale_x_continuous(expand = expansion(mult = 0.08)) +
      labs(title = sprintf("%s  (N=%d)", arm, n_p),
           x = "Days since D0",
           y = if (i == 1) "P(relapse)" else NULL) +
      annotate("text", x = -Inf, y = 0.5, label = "uncertain\n(0.5)",
               hjust = -0.1, size = 2, colour = "grey55", fontface = "italic") +
      theme_cd +
      theme(legend.position = if (show_leg) "right" else "none") +
      guides(fill   = guide_legend(override.aes = list(shape = 21, size = 3, colour = "grey20")),
             shape  = guide_legend(override.aes = list(fill  = "grey50", colour = "grey20", size = 3)),
             colour = guide_legend(override.aes = list(linewidth = 2)))
    p
  })
  names(arm_plots) <- arm_levels

  combined <- patchwork::wrap_plots(arm_plots, nrow = 1) +
    patchwork::plot_annotation(
      title    = sprintf("%s \u2014 Relapse probability: D0 \u2192 first recurrence", TRIAL_UP),
      subtitle = paste0("Point fill = DCiFER r  |  Shape = eMOI  |  ",
                        "Dashed connector = indeterminate  |  ",
                        "Faded points = uncertain Pv3Rs call"),
      theme = theme(plot.title    = element_text(size = 13, face = "bold"),
                    plot.subtitle = element_text(size = 8, colour = "grey40"))
    )
  fname <- sprintf("relapse_connected_dots_first_rec_%s.png", TRIAL_UP)
  ggsave(file.path(figs_dir, fname), combined,
         width = max(5, n_arms * 3.8), height = 4.8,
         dpi = 300, bg = "white", limitsize = FALSE)
  message(">> wrote: ", fname)
  invisible(combined)
}

# DISABLED — see note at top of script:
# plot_connected_dots_first_rec(bl)

# ── 8. Patient trajectory timeline ────────────────────────────────────────────
plot_patient_timelines <- function(con, min_steps = 2L) {
  if (is.null(con)) {
    message(">> patient_timelines: no consecutive data — skipping")
    return(invisible(NULL))
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message(">> patchwork required for patient timelines — skipping")
    return(invisible(NULL))
  }
  if (!all(c("from_day","to_day","pv3rs_L","id","treatment_arm") %in% names(con))) {
    message(">> patient_timelines: missing required columns — skipping")
    return(invisible(NULL))
  }

  df <- con %>%
    mutate(
      pv3rs_L   = suppressWarnings(as.numeric(pv3rs_L)),
      from_day  = suppressWarnings(as.numeric(from_day)),
      to_day    = suppressWarnings(as.numeric(to_day)),
      emoi_s1   = suppressWarnings(as.numeric(emoi_s1)),
      emoi_s2   = suppressWarnings(as.numeric(emoi_s2)),
      step_call = dplyr::case_when(
        is.na(pv3rs_L)    ~ "indeterminate",
        pv3rs_L >= 0.75   ~ "relapse",
        pv3rs_L <= 0.25   ~ "reinfection",
        TRUE              ~ "indeterminate"
      ),
      shape_s1 = dplyr::case_when(
        is.na(emoi_s1)   ~ "unknown",
        emoi_s1 <= 1.1   ~ "monoclonal",
        TRUE             ~ "polyclonal"
      ),
      shape_s2 = dplyr::case_when(
        is.na(emoi_s2)   ~ "unknown",
        emoi_s2 <= 1.1   ~ "monoclonal",
        TRUE             ~ "polyclonal"
      )
    ) %>%
    filter(is.finite(from_day), is.finite(to_day))

  step_counts <- df %>%
    group_by(id) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    filter(n >= min_steps)

  df <- df %>% filter(id %in% step_counts$id)
  if (!nrow(df)) {
    message(sprintf(">> patient_timelines: no patients with >= %d steps — skipping", min_steps))
    return(invisible(NULL))
  }

  arm_levels <- sort(unique(na.omit(as.character(df$treatment_arm))))

  step_pal <- c(
    relapse       = "#6aa9ff",
    reinfection   = "#2fbf71",
    indeterminate = "#e74c3c"
  )
  step_labels <- c(
    relapse       = "Relapse  p(L) \u22650.75",
    reinfection   = "Reinfection  p(L) \u22640.25",
    indeterminate = "Indeterminate"
  )
  step_lty <- c(
    relapse       = "solid",
    reinfection   = "solid",
    indeterminate = "dashed"
  )

  arm_plots <- lapply(arm_levels, function(arm) {
    d <- df %>%
      filter(treatment_arm == arm) %>%
      mutate(step_call = factor(step_call, levels = c("relapse","reinfection","indeterminate")))

    pat_order <- d %>%
      group_by(id) %>%
      summarise(first = min(from_day), last = max(to_day), .groups = "drop") %>%
      arrange(desc(last), desc(first)) %>%
      pull(id)
    d <- d %>% mutate(id = factor(id, levels = pat_order))

    nodes_from <- d %>%
      transmute(id, day = from_day, shape_lbl = shape_s1) %>%
      distinct()
    nodes_to <- d %>%
      transmute(id, day = to_day, shape_lbl = shape_s2) %>%
      distinct()
    nodes <- bind_rows(nodes_from, nodes_to) %>%
      distinct(id, day, .keep_all = TRUE) %>%
      mutate(shape_lbl = factor(shape_lbl, levels = c("monoclonal","polyclonal","unknown")))

    n_pts <- dplyr::n_distinct(d$id)

    ggplot() +
      geom_segment(data = d,
                   aes(x = from_day, xend = to_day,
                       y = id,       yend = id,
                       colour   = step_call,
                       linetype = step_call),
                   linewidth = 1.2, lineend = "round") +
      geom_point(data = nodes,
                 aes(x = day, y = id, shape = shape_lbl),
                 fill = "white", colour = "grey30",
                 size = 2.5, stroke = 0.7) +
      scale_colour_manual(
        values = step_pal, labels = step_labels,
        name   = "Classification",
        drop   = FALSE,
        guide  = guide_legend(
          title.theme  = element_text(size = 11, face = "bold"),
          label.theme  = element_text(size = 10),
          override.aes = list(linewidth = 2, shape = NA)
        )
      ) +
      scale_linetype_manual(values = step_lty, guide = "none") +
      scale_shape_manual(
        values = c(monoclonal = 21L, polyclonal = 24L, unknown = 23L),
        name   = "eMOI",
        guide  = guide_legend(
          title.theme  = element_text(size = 11, face = "bold"),
          label.theme  = element_text(size = 10),
          override.aes = list(size = 3, fill = "white", colour = "grey30")
        )
      ) +
      scale_x_continuous(
        name   = if (arm == arm_levels[length(arm_levels)]) "Days since D0" else NULL,
        expand = expansion(mult = 0.05)
      ) +
      labs(title = sprintf("%s  (N=%d)", arm, n_pts)) +
      theme_pub +
      theme(
        axis.title.y    = element_blank(),
        axis.text.y     = element_text(size = 7),
        plot.title      = element_text(size = 11, face = "bold"),
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
        panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
        panel.grid.minor   = element_blank(),
        legend.position = "none"
      )
  })

  arm_plots[[length(arm_plots)]] <- arm_plots[[length(arm_plots)]] +
    theme(legend.position = "right")

  n_per_arm <- df %>%
    group_by(treatment_arm) %>%
    summarise(n = dplyr::n_distinct(id), .groups = "drop") %>%
    arrange(match(treatment_arm, arm_levels)) %>%
    pull(n)
  heights <- pmax(n_per_arm, 1)

  combined <- patchwork::wrap_plots(arm_plots, ncol = 1, heights = heights) +
    patchwork::plot_annotation(
      title = sprintf("%s — Patient recurrence trajectories (consecutive pairs)", TRIAL_UP),
      subtitle = sprintf(
        "Each row = one patient | Line colour = Pv3Rs classification | Shape = eMOI | min. %d steps",
        min_steps),
      theme = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9,  colour = "grey40")
      )
    )

  pw <- 8
  ph <- max(6, sum(heights) * 0.35 + 2)

  fname <- sprintf("patient_timelines_%s.png", TRIAL_UP)
  ggsave(file.path(figs_dir, fname), combined,
         width = pw, height = ph, dpi = 300, bg = "white", limitsize = FALSE)
  message(">> wrote: ", fname)
  invisible(combined)
}

plot_patient_timelines(con, min_steps = 2L)

# ── 9. Reinfection rank plot — first-rec-only baseline data only ─────────────
plot_reinfection_rank <- function(df) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message(">> patchwork required for reinfection rank plot — skipping")
    return(invisible(NULL))
  }
  if (is.null(df)) return(invisible(NULL))
  need <- c("treatment_arm", "pv3rs_I")
  missing_cols <- setdiff(need, names(df))
  if (length(missing_cols)) {
    message(">> reinfection rank plot: missing columns — ",
            paste(missing_cols, collapse = ", "), " — skipping")
    return(invisible(NULL))
  }
  series_pal <- c("Pv3Rs" = "#377EB8")
  d <- df %>%
    filter(!is.na(treatment_arm)) %>%
    mutate(
      p_reinf_pv3rs = suppressWarnings(as.numeric(pv3rs_I)),
      p_reinf_mid   = dplyr::case_when(
        is.finite(p_reinf_pv3rs) ~ p_reinf_pv3rs,
        TRUE ~ NA_real_
      )
    ) %>%
    group_by(treatment_arm) %>%
    mutate(rank = rank(p_reinf_mid, ties.method = "first", na.last = "keep")) %>%
    ungroup() %>%
    filter(is.finite(rank))
  if (!nrow(d)) return(invisible(NULL))

  arm_levels <- sort(unique(na.omit(as.character(d$treatment_arm))))
  segs <- d %>%
    group_by(treatment_arm, rank) %>%
    summarise(xmin = min(rank, na.rm = TRUE),
              xmax = max(rank, na.rm = TRUE), .groups = "drop") %>%
    mutate(treatment_arm = factor(treatment_arm, levels = arm_levels))
  d_long <- d %>%
    tidyr::pivot_longer(cols = starts_with("p_reinf_pv3rs"),
                        names_to = "series", values_to = "p_reinf") %>%
    mutate(series        = factor("Pv3Rs", levels = names(series_pal)),
           treatment_arm = factor(treatment_arm, levels = arm_levels)) %>%
    filter(is.finite(p_reinf))
  if (!nrow(d_long)) return(invisible(NULL))

  arm_n <- d %>%
    group_by(treatment_arm) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    mutate(treatment_arm = factor(treatment_arm, levels = arm_levels))

  arm_plots <- lapply(seq_along(arm_levels), function(i) {
    arm     <- arm_levels[i]
    show_lg <- i == length(arm_levels)
    pd <- d_long %>% filter(treatment_arm == arm)
    sg <- segs   %>% filter(treatment_arm == arm)
    an <- arm_n  %>% filter(treatment_arm == arm)
    ggplot() +
      geom_segment(data = sg,
                   aes(x = xmin, xend = xmax, y = rank, yend = rank),
                   colour = "grey60", linetype = "dotted", linewidth = 0.5) +
      geom_point(data = pd,
                 aes(x = p_reinf, y = rank, colour = series),
                 size = 1.8, alpha = 0.8) +
      geom_text(data = an,
                aes(x = Inf, y = -Inf, label = sprintf("N=%d", n)),
                hjust = 1.1, vjust = -0.5, size = 2.5, colour = "grey45",
                inherit.aes = FALSE) +
      scale_colour_manual(values = series_pal, name = NULL) +
      scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
      labs(title = arm,
           x = "Reinfection probability",
           y = if (i == 1) "First recurrence ranked\nby mean reinfection estimate" else NULL) +
      theme_pub +
      theme(legend.position  = if (show_lg) "right" else "none",
            axis.text.y      = element_text(size = 6),
            axis.title.y     = if (i == 1) element_text() else element_blank(),
            plot.title       = element_text(face = "bold", size = 10, hjust = 0.5))
  })

  combined <- patchwork::wrap_plots(arm_plots, nrow = 1) +
    patchwork::plot_annotation(
      title    = sprintf("%s \u2014 First recurrence: reinfection probability (Pv3Rs)", TRIAL_UP),
      subtitle = "Patients ranked within each arm by reinfection probability",
      theme = theme(plot.title    = element_text(size = 13, face = "bold"),
                    plot.subtitle = element_text(size = 8,  colour = "grey40"))
    )
  fname <- sprintf("reinfection_rank_first_rec_%s.png", TRIAL_UP)
  ggsave(file.path(figs_dir, fname), combined,
         width = length(arm_levels) * 3.2 + 1.5, height = 6,
         dpi = 300, bg = "white", limitsize = FALSE)
  message(">> wrote: ", fname)
  invisible(combined)
}

if (bl_is_first_rec) {
  plot_reinfection_rank(bl)
}

# ── 10. Relapse-averted estimate plot ────────────────────────────────────────
# NOTE: plot_relapse_averted() remains defined but is not called below (see
# "DISABLED"). Left intact in case you want to re-enable it.
plot_relapse_averted <- function(df,
                                  class_col     = "pv3rs_call",
                                  relapse_vals  = c("relapse"),
                                  method_label  = "Pv3Rs",
                                  ref_arm       = NULL,
                                  pq_arms       = NULL,
                                  arm_enrolled  = NULL) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message(">> patchwork required for relapse_averted plot — skipping")
    return(invisible(NULL))
  }
  if (is.null(df) || !has_cols(df, c("treatment_arm", class_col, "id")))
    return(invisible(NULL))

  if (is.null(pq_arms)) {
    message(">> relapse_averted: pq_arms not supplied — all arms used in Panel B")
  }

  dat <- df %>%
    filter(!is.na(.data[[class_col]]),
           !is.na(treatment_arm),
           .data[[class_col]] %in% c(relapse_vals, "reinfection")) %>%
    group_by(treatment_arm) %>%
    summarise(
      n_relapse    = sum(.data[[class_col]] %in% relapse_vals),
      n_classified = dplyr::n(),
      n_indiv      = dplyr::n_distinct(id),
      prop         = n_relapse / n_classified,
      .groups      = "drop"
    )

  if (!nrow(dat) || nrow(dat) < 2) {
    message(">> relapse_averted: fewer than 2 arms with data — skipping")
    return(invisible(NULL))
  }

  if (is.null(ref_arm)) {
    non_pq_dat <- if (!is.null(pq_arms))
      dat %>% filter(!treatment_arm %in% pq_arms)
    else
      dat
    if (!nrow(non_pq_dat)) non_pq_dat <- dat
    ref_arm <- non_pq_dat$treatment_arm[which.max(non_pq_dat$prop)]
    message(sprintf(">> relapse_averted: reference arm (non-PQ, highest relapse) = '%s'",
                    ref_arm))
  }
  if (!ref_arm %in% dat$treatment_arm) {
    message(">> relapse_averted: ref_arm '", ref_arm, "' not in data — skipping")
    return(invisible(NULL))
  }

  p_ref <- dat$prop[dat$treatment_arm == ref_arm]

  pq_in_data <- if (!is.null(pq_arms))
    dat %>% filter(treatment_arm %in% pq_arms)
  else
    dat %>% filter(treatment_arm != ref_arm)

  if (!nrow(pq_in_data)) {
    message(">> relapse_averted: no PQ arms found in data — skipping Panel B")
    return(invisible(NULL))
  }

  get_enrolled <- function(arm) {
    if (!is.null(arm_enrolled) && arm %in% names(arm_enrolled))
      as.integer(arm_enrolled[[arm]])
    else
      dat$n_classified[dat$treatment_arm == arm]
  }

  averted <- pq_in_data %>%
    rowwise() %>%
    mutate(
      p_diff        = p_ref - prop,
      n_enrolled    = get_enrolled(treatment_arm),
      n_averted_est = p_diff * n_enrolled
    ) %>%
    ungroup()

  extra_note <- if (!is.null(arm_enrolled))
    "Enrolled N supplied — extrapolated to full arm"
  else
    "No enrolled N supplied — estimate uses classified pairs only (conservative)"

  arm_levels <- c(ref_arm, setdiff(dat$treatment_arm, ref_arm))
  arm_pal    <- get_arm_pal(arm_levels)
  pA <- ggplot(dat %>% mutate(treatment_arm = factor(treatment_arm, arm_levels)),
               aes(treatment_arm, prop, colour = treatment_arm)) +
    geom_point(size = 4) +
    geom_text(aes(label = sprintf("%d/%d", n_relapse, n_classified),
                  vjust = ifelse(prop >= 0.92, 2.2, -1.5)),
              size = 2.8, colour = "grey35") +
    geom_hline(yintercept = p_ref, linetype = "dashed",
               colour = arm_pal[[ref_arm]], linewidth = 0.7, alpha = 0.6) +
    scale_colour_manual(values = arm_pal, guide = "none") +
    scale_y_continuous(labels = percent, limits = c(0, 1.07),
                       breaks = seq(0, 1, 0.25)) +
    labs(x = "Treatment arm",
         y = "Proportion classified as relapse",
         title = sprintf("A  Relapse proportion by arm — %s", method_label),
         subtitle = sprintf("Reference arm: %s (dashed line)", ref_arm)) +
    theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1))

  pq_arm_levels <- intersect(arm_levels, averted$treatment_arm)
  pB <- ggplot(averted %>% mutate(treatment_arm = factor(treatment_arm, pq_arm_levels)),
               aes(treatment_arm, n_averted_est, colour = treatment_arm)) +
    geom_point(size = 4) +
    geom_text(aes(label = sprintf("~%.0f", n_averted_est)),
              vjust = ifelse(averted$n_averted_est >= 0, -1.5, 2.2),
              size = 2.8, colour = "grey35") +
    geom_hline(yintercept = 0, linetype = "solid",
               colour = "grey60", linewidth = 0.4) +
    scale_colour_manual(values = arm_pal[pq_arm_levels], guide = "none") +
    scale_y_continuous(labels = function(x) round(x),
                       expand = expansion(mult = c(0.15, 0.15))) +
    labs(x = "Treatment arm (PQ only)",
         y = "Estimated relapses averted\n(vs reference arm)",
         title = "B  Estimated relapses averted (PQ arms only)",
         subtitle = paste0(extra_note, "\nNon-PQ arms excluded — cannot clear hypnozoites")) +
    theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1))

  combined <- (pA | pB) +
    patchwork::plot_annotation(
      title    = sprintf("%s \u2014 Relapse burden and averted relapses (%s)",
                         TRIAL_UP, method_label),
      subtitle = paste0(
        "\u26a0\ufe0f  EXPLORATORY ESTIMATE — interpret with caution:\n",
        "(1) Classified subset is biased toward PQ non-responders ",
        "(most PQ patients never recurred and have no pair).\n",
        "(2) Assumes classification proportions in classified pairs represent ",
        "the full enrolled arm — this assumption is not validated.\n",
        "(3) CIs are binomial on classified pairs only, not on enrolled population."
      ),
      theme = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "#8B4513",
                                     lineheight = 1.3)
      )
    )

  fname <- sprintf("relapse_averted_%s_%s.png", tolower(method_label), TRIAL_UP)
  ggsave(file.path(figs_dir, fname), combined,
         width = 10, height = 5.5, dpi = 300, bg = "white")
  message(">> wrote: ", fname)
  invisible(combined)
}

# DISABLED — see note at top of script:
# arm_n_enrolled <- NULL
# if (bl_is_first_rec && !is.null(bl)) {
#   if (is.null(arm_n_enrolled) && "treatment_arm" %in% names(bl)) {
#     d0_enrolled <- bl %>%
#       filter(if ("recurrence_no" %in% names(.)) recurrence_no == 1L else TRUE) %>%
#       count(treatment_arm) %>%
#       tibble::deframe()
#     arm_n_enrolled <- d0_enrolled
#   }
#   plot_relapse_averted(
#     df           = bl,
#     class_col    = "pv3rs_call",
#     relapse_vals = c("relapse"),
#     method_label = "Pv3Rs",
#     ref_arm      = NULL,
#     pq_arms      = if ("PQ" %in% arm_group_map) "PQ" else NULL,
#     arm_enrolled = arm_n_enrolled
#   )
# }

# ── END ───────────────────────────────────────────────────────────────────────
cat(sprintf(">> %s figures written to:\n   %s\n", TRIAL_UP, figs_dir))
