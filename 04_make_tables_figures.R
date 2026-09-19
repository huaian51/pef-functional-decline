# ============================================================================
# 04_make_tables_figures.R
# Exact manuscript-facing tables and figures for the public release.
#
# This orchestrator:
#   1) builds the locked manuscript Table 1 and Table 2;
#   2) creates a temporary compatibility layer whose filenames/columns match
#      the locked final plotting scripts;
#   3) executes public-path versions of the final locked Figure 1-3 and
#      Supplementary Figure S1-S3 scripts.
#
# The figure modules in R/locked_figures/ are derived directly from the final
# manuscript plotting scripts. Their analytical/plotting logic is unchanged;
# only local-drive paths and input filenames are redirected to this repository.
# ============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, warn = 1, width = 320, scipen = 999)

suppressPackageStartupMessages(library(data.table))

COHORTS <- c("HRS", "ELSA", "SHARE", "CHARLS")

DATA_DIR <- file.path("output", "derived_data")
TABLE_DIR <- file.path("output", "tables")
QC_DIR <- file.path("output", "QC")
FIG_DIR <- file.path("output", "figures")
COMPAT_ROOT <- file.path("output", "compat")
U30_ROOT <- file.path(COMPAT_ROOT, "unified_burden_check")
U31_ROOT <- file.path(COMPAT_ROOT, "unified_burden_final")

for (d in c(
  TABLE_DIR, QC_DIR, FIG_DIR,
  file.path(U30_ROOT, "data"),
  file.path(U30_ROOT, "tables"),
  file.path(U30_ROOT, "QC"),
  file.path(U31_ROOT, "tables"),
  file.path(U31_ROOT, "QC")
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

write_bom <- function(x, path) {
  fwrite(x, path, bom = TRUE, na = "")
  message("Saved: ", path)
}

fmt_n <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE)
}

fmt_n_pct <- function(n, denom, digits = 1L) {
  paste0(
    fmt_n(n), " (",
    formatC(100 * n / denom, format = "f", digits = digits),
    "%)"
  )
}

fmt_or_ci <- function(est, low, high, digits = 2L) {
  paste0(
    formatC(est, format = "f", digits = digits),
    " (",
    formatC(low, format = "f", digits = digits),
    "-",
    formatC(high, format = "f", digits = digits),
    ")"
  )
}

detect_female <- function(x, cohort) {
  sx <- trimws(tolower(as.character(x)))

  if (any(sx %in% c("female", "woman", "women", "f"), na.rm = TRUE)) {
    return(as.integer(sx %in% c("female", "woman", "women", "f")))
  }

  lev <- sort(unique(sx[!is.na(sx)]))

  # Locked SHARE exception only: gender 1=Male, 2=Female.
  if (identical(cohort, "SHARE") && setequal(lev, c("1", "2"))) {
    return(as.integer(sx == "2"))
  }

  stop(
    "[", cohort,
    "] Female level cannot be identified safely. Observed sex values: ",
    paste(lev, collapse = ", "),
    call. = FALSE
  )
}

get_first_person_row <- function(d, cohort) {
  z <- copy(d)
  z[, .source_row__ := .I]

  if ("origin_wave" %in% names(z)) {
    suppressWarnings(
      z[, origin_wave_num__ := as.numeric(as.character(origin_wave))]
    )
    if (all(is.na(z$origin_wave_num__))) {
      z[, origin_wave_num__ := .source_row__]
    }
    setorder(z, id, origin_wave_num__, .source_row__, na.last = TRUE)
  } else if ("interval_f" %in% names(z)) {
    if (is.factor(z$interval_f)) {
      z[, interval_order__ := as.integer(interval_f)]
    } else {
      z[, interval_order__ := match(
        as.character(interval_f),
        unique(as.character(interval_f))
      )]
    }
    setorder(z, id, interval_order__, .source_row__, na.last = TRUE)
  } else {
    setorder(z, id, .source_row__)
  }

  out <- z[, .SD[1L], by = id]
  if (nrow(out) != uniqueN(d$id)) {
    stop("[", cohort, "] participant deduplication failed.", call. = FALSE)
  }
  out[]
}

# ============================================================================
# 1. LOCKED PRIMARY INPUTS
# ============================================================================

PRIMARY_FILES <- setNames(
  file.path(DATA_DIR, paste0("01_", COHORTS, "_primary_analysis_dataset.rds")),
  COHORTS
)

required_primary <- c(
  PRIMARY_FILES,
  file.path(TABLE_DIR, "02B_any_limitation.csv"),
  file.path(TABLE_DIR, "02C_multinomial.csv"),
  file.path(TABLE_DIR, "02E_direct_severe_vs_mild.csv"),
  file.path(TABLE_DIR, "02F_primary_meta.csv"),
  file.path(TABLE_DIR, "31J_low_count_supportive_meta.csv"),
  file.path(TABLE_DIR, "31F_smoking_extended_adjustment_meta.csv"),
  file.path(TABLE_DIR, "31O_IPW_observed_followup_meta.csv"),
  file.path(TABLE_DIR, "31K_severe_vs_nonsevere_results.csv"),
  file.path(TABLE_DIR, "31L_severe_vs_nonsevere_meta.csv"),
  file.path(TABLE_DIR, "31T_leave_one_cohort_out_meta.csv"),
  file.path(QC_DIR, "31U_FINAL_UNIFIED_SENSITIVITY_GATE.csv")
)

missing <- required_primary[!file.exists(required_primary)]
if (length(missing)) {
  stop(
    "Run scripts 01-03 first. Missing required file(s):\n",
    paste(missing, collapse = "\n"),
    call. = FALSE
  )
}

# Strict outcome is part of final Figure S3 and therefore required for the exact
# manuscript figure set.
strict_files <- c(
  file.path(TABLE_DIR, "31R_ge2_deterioration_results.csv"),
  file.path(TABLE_DIR, "31S_ge2_deterioration_meta.csv")
)
if (any(!file.exists(strict_files))) {
  stop(
    "The exact manuscript Figure S3 requires the strict >=2 deterioration module.\n",
    "Supply data_harmonized/<COHORT>_ge2_sensitivity.rds for all four cohorts, ",
    "rerun 03_sensitivity_analysis.R, then rerun this script.",
    call. = FALSE
  )
}

# ============================================================================
# 2. TABLE 1 — SAME SUMMARY RULES USED FOR THE FINAL MANUSCRIPT
# ============================================================================

EXPECTED <- fread(file.path("metadata", "expected_sample_counts.csv"))

summary_rows <- list()
qc_rows <- list()

for (coh in COHORTS) {
  d <- as.data.table(readRDS(PRIMARY_FILES[[coh]]))
  d[, destination := as.character(destination)]

  required <- c("id", "destination", "burden_z", "age", "sex", "bmi")
  miss <- setdiff(required, names(d))
  if (length(miss)) {
    stop("[", coh, "] missing Table 1 variable(s): ", paste(miss, collapse = ", "))
  }

  interval_N <- nrow(d)
  person_N <- uniqueN(d$id)
  independent_N <- sum(d$destination == "Independent", na.rm = TRUE)
  mild_N <- sum(d$destination == "Mild limitation", na.rm = TRUE)
  severe_N <- sum(d$destination == "Severe limitation", na.rm = TRUE)
  any_limitation_N <- mild_N + severe_N

  p <- get_first_person_row(d, coh)
  p[, female__ := detect_female(sex, coh)]

  summary_rows[[coh]] <- data.table(
    cohort = coh,
    interval_N = interval_N,
    person_N = person_N,
    age_mean = mean(as.numeric(p$age), na.rm = TRUE),
    age_sd = sd(as.numeric(p$age), na.rm = TRUE),
    female_N = sum(p$female__ == 1L, na.rm = TRUE),
    female_pct = 100 * mean(p$female__ == 1L, na.rm = TRUE),
    bmi_mean = mean(as.numeric(p$bmi), na.rm = TRUE),
    bmi_sd = sd(as.numeric(p$bmi), na.rm = TRUE),
    burden_mean_interval = mean(as.numeric(d$burden_z), na.rm = TRUE),
    burden_sd_interval = sd(as.numeric(d$burden_z), na.rm = TRUE),
    independent_N = independent_N,
    mild_N = mild_N,
    severe_N = severe_N,
    any_limitation_N = any_limitation_N
  )

  qc_rows[[coh]] <- data.table(
    cohort = coh,
    participant_rows = nrow(p),
    participant_unique_ids = uniqueN(p$id),
    first_interval_rule = if ("origin_wave" %in% names(d)) {
      "Earliest origin_wave"
    } else {
      "Earliest interval_f / row order"
    }
  )
}

T1_SUMMARY <- rbindlist(summary_rows, fill = TRUE)
T1_QC <- rbindlist(qc_rows, fill = TRUE)

lock_qc <- merge(
  T1_SUMMARY[, .(
    cohort,
    interval_N,
    person_N,
    independent_N,
    mild_N,
    severe_N
  )],
  EXPECTED,
  by = "cohort",
  all.x = TRUE
)

lock_qc[, exact_match :=
  interval_N == interval_N_expected &
  person_N == person_N_expected &
  independent_N == independent_N_expected &
  mild_N == mild_N_expected &
  severe_N == severe_N_expected
]

if (!all(lock_qc$exact_match)) {
  print(lock_qc)
  stop("Locked primary-sample QC failed. Table 1 generation stopped.")
}

T1_LONG <- rbindlist(
  lapply(COHORTS, function(coh) {
    z <- T1_SUMMARY[cohort == coh]
    data.table(
      characteristic = c(
        "Participants, n",
        "Person-intervals, n",
        "Age, years",
        "Female, n (%)",
        "BMI, kg/m2",
        "Unified burden z",
        "Independent, n (%)",
        "Mild, n (%)",
        "Severe, n (%)",
        "Any limitation, n (%)"
      ),
      cohort = coh,
      value = c(
        fmt_n(z$person_N),
        fmt_n(z$interval_N),
        paste0(
          formatC(z$age_mean, format = "f", digits = 1),
          " (", formatC(z$age_sd, format = "f", digits = 1), ")"
        ),
        fmt_n_pct(z$female_N, z$person_N),
        paste0(
          formatC(z$bmi_mean, format = "f", digits = 1),
          " (", formatC(z$bmi_sd, format = "f", digits = 1), ")"
        ),
        paste0(
          formatC(z$burden_mean_interval, format = "f", digits = 3),
          " (", formatC(z$burden_sd_interval, format = "f", digits = 3), ")"
        ),
        fmt_n_pct(z$independent_N, z$interval_N),
        fmt_n_pct(z$mild_N, z$interval_N),
        fmt_n_pct(z$severe_N, z$interval_N),
        fmt_n_pct(z$any_limitation_N, z$interval_N)
      )
    )
  }),
  fill = TRUE
)

T1_DISPLAY <- dcast(
  T1_LONG,
  characteristic ~ cohort,
  value.var = "value"
)

T1_DISPLAY[, characteristic := factor(
  characteristic,
  levels = c(
    "Participants, n",
    "Person-intervals, n",
    "Age, years",
    "Female, n (%)",
    "BMI, kg/m2",
    "Unified burden z",
    "Independent, n (%)",
    "Mild, n (%)",
    "Severe, n (%)",
    "Any limitation, n (%)"
  )
)]
setorder(T1_DISPLAY, characteristic)
T1_DISPLAY[, characteristic := as.character(characteristic)]

write_bom(T1_SUMMARY, file.path(TABLE_DIR, "04A_Table1_summary.csv"))
write_bom(T1_DISPLAY, file.path(TABLE_DIR, "04B_Table1_display.csv"))
write_bom(T1_QC, file.path(QC_DIR, "04A_Table1_participant_QC.csv"))
write_bom(lock_qc, file.path(QC_DIR, "04B_Table1_locked_sample_QC.csv"))

# ============================================================================
# 3. TABLE 2 — PRIMARY ASSOCIATIONS
# ============================================================================

A <- fread(file.path(TABLE_DIR, "02B_any_limitation.csv"))
M <- fread(file.path(TABLE_DIR, "02C_multinomial.csv"))
D <- fread(file.path(TABLE_DIR, "02E_direct_severe_vs_mild.csv"))
META <- fread(file.path(TABLE_DIR, "02F_primary_meta.csv"))

cohort_rows <- rbindlist(
  lapply(COHORTS, function(coh) {
    data.table(
      Cohort = coh,
      `Any limitation vs Independent` = fmt_or_ci(
        A[cohort == coh, OR],
        A[cohort == coh, CI_low],
        A[cohort == coh, CI_high]
      ),
      `Mild vs Independent` = fmt_or_ci(
        M[cohort == coh & estimand == "Mild vs Independent", effect],
        M[cohort == coh & estimand == "Mild vs Independent", CI_low_percentile],
        M[cohort == coh & estimand == "Mild vs Independent", CI_high_percentile]
      ),
      `Severe vs Independent` = fmt_or_ci(
        M[cohort == coh & estimand == "Severe vs Independent", effect],
        M[cohort == coh & estimand == "Severe vs Independent", CI_low_percentile],
        M[cohort == coh & estimand == "Severe vs Independent", CI_high_percentile]
      ),
      `Severe vs Mild` = fmt_or_ci(
        D[cohort == coh, OR],
        D[cohort == coh, CI_low],
        D[cohort == coh, CI_high]
      ),
      `P for Severe vs Mild` = ifelse(
        D[cohort == coh, P] < 0.001,
        "<0.001",
        formatC(D[cohort == coh, P], format = "f", digits = 3)
      )
    )
  }),
  fill = TRUE
)

pooled <- data.table(
  Cohort = "Pooled",
  `Any limitation vs Independent` = with(
    META[estimand == "Any limitation vs Independent"],
    fmt_or_ci(pooled_OR, CI_low_mKH, CI_high_mKH)
  ),
  `Mild vs Independent` = with(
    META[estimand == "Mild vs Independent"],
    fmt_or_ci(pooled_OR, CI_low_mKH, CI_high_mKH)
  ),
  `Severe vs Independent` = with(
    META[estimand == "Severe vs Independent"],
    fmt_or_ci(pooled_OR, CI_low_mKH, CI_high_mKH)
  ),
  `Severe vs Mild` = with(
    META[estimand == "Severe vs Mild"],
    fmt_or_ci(pooled_OR, CI_low_mKH, CI_high_mKH)
  ),
  `P for Severe vs Mild` = with(
    META[estimand == "Severe vs Mild"],
    ifelse(
      P_mKH < 0.001,
      "<0.001",
      formatC(P_mKH, format = "f", digits = 3)
    )
  )
)

T2 <- rbind(cohort_rows, pooled, fill = TRUE)
write_bom(T2, file.path(TABLE_DIR, "04C_Table2_display.csv"))

# ============================================================================
# 4. BUILD COMPATIBILITY INPUTS FOR THE LOCKED FINAL FIGURE SCRIPTS
# ============================================================================

# 4.1 Primary RDS aliases (same objects, manuscript-final filenames).
for (coh in COHORTS) {
  target <- file.path(
    U30_ROOT, "data",
    paste0("30_", coh, "_primary_analysis_UNIFIED_burden.rds")
  )
  file.copy(PRIMARY_FILES[[coh]], target, overwrite = TRUE)
}

# 4.2 Script-30-style primary result files.
A30 <- copy(A)
A30[, effect := OR]
write_bom(A30, file.path(U30_ROOT, "tables", "30C_any_limitation_unified_burden.csv"))

C30 <- copy(M)
C30[, point_beta := beta]
C30[, boot_SE := SE_boot]
write_bom(C30, file.path(U30_ROOT, "tables", "30D_multinomial_unified_burden.csv"))

D30 <- copy(D)
D30[, effect := OR]
write_bom(D30, file.path(U30_ROOT, "tables", "30F_direct_severe_vs_mild_unified_burden.csv"))

M30 <- copy(META)
M30[, estimand := fcase(
  estimand == "Any limitation vs Independent", "A: Any limitation vs Independent",
  estimand == "Mild vs Independent", "C: Mild vs Independent",
  estimand == "Severe vs Independent", "C: Severe vs Independent",
  estimand == "Severe vs Mild", "Severity contrast: Severe vs Mild",
  default = estimand
)]
write_bom(M30, file.path(U30_ROOT, "tables", "30G_meta_unified_burden.csv"))

# 4.3 Primary locked-result gate.
expected_primary <- fread(file.path("metadata", "expected_primary_results.csv"))

observed_primary <- rbindlist(
  list(
    A[, .(type = "cohort", cohort, estimand = "Any limitation vs Independent", observed_OR = OR)],
    M[estimand == "Mild vs Independent",
      .(type = "cohort", cohort, estimand, observed_OR = effect)],
    M[estimand == "Severe vs Independent",
      .(type = "cohort", cohort, estimand, observed_OR = effect)],
    D[, .(type = "cohort", cohort, estimand = "Severe vs Mild", observed_OR = OR)],
    META[, .(type = "pooled", cohort = "ALL", estimand, observed_OR = pooled_OR)]
  ),
  fill = TRUE
)

primary_audit <- merge(
  expected_primary,
  observed_primary,
  by = c("type", "cohort", "estimand"),
  all.x = TRUE
)
primary_audit[, abs_diff := abs(observed_OR - expected_OR)]
primary_result_ok <- all(primary_audit$abs_diff <= 5e-5, na.rm = FALSE)

primary_gate <- data.table(
  check = c(
    "Locked primary sample counts reproduced",
    "Locked cohort-specific and pooled primary ORs reproduced"
  ),
  pass = c(
    all(lock_qc$exact_match),
    primary_result_ok
  ),
  value = c(
    paste0(sum(lock_qc$exact_match), "/4 cohort samples exact"),
    paste0(
      "max |OR difference| = ",
      formatC(max(primary_audit$abs_diff, na.rm = TRUE), format = "f", digits = 6)
    )
  )
)
write_bom(primary_gate, file.path(U30_ROOT, "QC", "30J_unified_burden_stability_gate.csv"))
write_bom(primary_audit, file.path(QC_DIR, "04C_primary_locked_result_audit.csv"))

# 4.4 Script-31 final sensitivity aliases.
u31_table_files <- c(
  "31E_smoking_extended_adjustment_results.csv",
  "31F_smoking_extended_adjustment_meta.csv",
  "31I_low_count_supportive_results.csv",
  "31J_low_count_supportive_meta.csv",
  "31K_severe_vs_nonsevere_results.csv",
  "31L_severe_vs_nonsevere_meta.csv",
  "31N_IPW_observed_followup_results.csv",
  "31O_IPW_observed_followup_meta.csv",
  "31R_ge2_deterioration_results.csv",
  "31S_ge2_deterioration_meta.csv",
  "31T_leave_one_cohort_out_meta.csv"
)

for (f in u31_table_files) {
  file.copy(
    file.path(TABLE_DIR, f),
    file.path(U31_ROOT, "tables", f),
    overwrite = TRUE
  )
}

u31_qc_files <- c(
  "31U_FINAL_UNIFIED_SENSITIVITY_GATE.csv"
)
for (f in u31_qc_files) {
  file.copy(
    file.path(QC_DIR, f),
    file.path(U31_ROOT, "QC", f),
    overwrite = TRUE
  )
}

# ============================================================================
# 5. RUN PUBLIC-PATH COPIES OF THE LOCKED FINAL FIGURE SCRIPTS
# ============================================================================

locked_scripts <- c(
  "04_Fig1_Study_Design_v5_UNIFIED_FINAL_PUBLIC.R",
  "05_Fig2_Primary_Results_v4_UNIFIED_FINAL_PUBLIC.R",
  "06_Fig3_Repeated_Low_PEF_Supportive_v3_UNIFIED_FINAL_PUBLIC.R",
  "07_FigS1_Leave_One_Cohort_Out_Sensitivity_v3_UNIFIED_FINAL_PUBLIC.R",
  "08_FigS2_Pooled_Sensitivity_Analyses_v3_UNIFIED_FINAL_PUBLIC.R",
  "09_FigS3_Alternative_Outcome_Robustness_v2_UNIFIED_FINAL_PUBLIC.R"
)

for (s in locked_scripts) {
  f <- file.path("R", "locked_figures", s)
  if (!file.exists(f)) stop("Missing locked figure module: ", f)

  message("\n===============================================================================")
  message("Running locked figure module: ", s)
  message("===============================================================================")

  # New environment isolates each original plotting script, which begins by
  # clearing its own workspace.
  sys.source(
    f,
    envir = new.env(parent = globalenv()),
    keep.source = FALSE
  )
}

cat("\n")
cat("===============================================================================\n")
cat("TABLES AND LOCKED FINAL FIGURES COMPLETE\n")
cat("===============================================================================\n")
cat("Tables: ", normalizePath(TABLE_DIR, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Figures: ", normalizePath(FIG_DIR, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("===============================================================================\n")
