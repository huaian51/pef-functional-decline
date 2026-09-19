# ============================================================================
# 02_primary_analysis.R
# Primary statistical analysis for the public release.
# ============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("nnet", quietly = TRUE)) stop("Package 'nnet' is required.")
if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' is required.")
source(file.path("R", "helpers.R"))

COHORTS <- c("HRS", "ELSA", "SHARE", "CHARLS")
CORE0 <- c("age", "sex", "interval_f", "height", "bmi", "education")
B_BOOT <- 1000L
SEED_BASE <- 20260919L

DATA_DIR <- file.path("output", "derived_data")
TABLE_DIR <- file.path("output", "tables")
QC_DIR <- file.path("output", "QC")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

FILES <- setNames(
  file.path(DATA_DIR, paste0("01_", COHORTS, "_primary_analysis_dataset.rds")),
  COHORTS
)
missing <- FILES[!file.exists(FILES)]
if (length(missing)) stop("Run 01_build_analysis_datasets.R first. Missing:\n", paste(missing, collapse = "\n"))

D <- setNames(lapply(FILES, function(f) as.data.table(readRDS(f))), COHORTS)

# ---- Sample QC ----------------------------------------------------------------
sample_qc <- rbindlist(lapply(COHORTS, function(coh) {
  d <- D[[coh]]
  data.table(
    cohort = coh,
    interval_N = nrow(d),
    person_N = uniqueN(d$id),
    Independent_N = sum(d$destination == "Independent", na.rm = TRUE),
    Mild_N = sum(d$destination == "Mild limitation", na.rm = TRUE),
    Severe_N = sum(d$destination == "Severe limitation", na.rm = TRUE)
  )
}))
write_csv_utf8(sample_qc, file.path(QC_DIR, "02A_primary_sample_QC.csv"))

expected_sample_file <- file.path("metadata", "expected_sample_counts.csv")
if (file.exists(expected_sample_file)) {
  expected_sample <- fread(expected_sample_file)
  sample_compare <- merge(expected_sample, sample_qc, by = "cohort", all.x = TRUE, suffixes = c("_expected", "_observed"))
  sample_compare[, exact_match :=
    interval_N_expected == interval_N_observed &
    person_N_expected == person_N_observed &
    independent_N_expected == Independent_N &
    mild_N_expected == Mild_N &
    severe_N_expected == Severe_N
  ]
  write_csv_utf8(sample_compare, file.path(QC_DIR, "02A2_locked_sample_comparison.csv"))
}

# ---- A. Any limitation vs Independent ----------------------------------------
any_rows <- lapply(COHORTS, function(coh) {
  d <- D[[coh]]
  fit_cluster_logit(
    d, outcome = "any_limitation", exposure = "burden_z", covars = CORE0,
    cohort = coh, analysis_label = "A_ANY_LIMITATION", id_col = "id"
  )
})
any_results <- rbindlist(any_rows, fill = TRUE)
write_csv_utf8(any_results, file.path(TABLE_DIR, "02B_any_limitation.csv"))

# ---- B. Multinomial Mild / Severe vs Independent -----------------------------
multinom_rows <- list()
boot_qc_rows <- list()

for (i in seq_along(COHORTS)) {
  coh <- COHORTS[i]
  d <- D[[coh]]
  check_required(d, c("id", "dest3", "burden_z", CORE0), paste0(coh, " primary data"))

  fit <- fit_multinom(d, exposure = "burden_z", covars = CORE0)
  b0 <- extract_multinom_exposure(fit, exposure = "burden_z")

  boot <- bootstrap_multinom_person(
    d, exposure = "burden_z", covars = CORE0,
    B = B_BOOT, seed = SEED_BASE + i, progress_every = 100L
  )

  valid_mild <- is.finite(boot[, "beta_mild"])
  valid_severe <- is.finite(boot[, "beta_severe"])
  valid_contrast <- is.finite(boot[, "beta_severe_minus_mild"])

  boot_qc_rows[[coh]] <- data.table(
    cohort = coh,
    B_requested = B_BOOT,
    valid_mild_B = sum(valid_mild),
    valid_severe_B = sum(valid_severe),
    valid_contrast_B = sum(valid_contrast)
  )

  make_row <- function(label, beta0, boot_vec) {
    z <- boot_vec[is.finite(boot_vec)]
    se <- sd(z)
    qs <- quantile(z, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    data.table(
      cohort = coh,
      estimand = label,
      interval_N = nrow(d),
      person_N = uniqueN(d$id),
      beta = beta0,
      SE_boot = se,
      effect = exp(beta0),
      CI_low_percentile = exp(qs[1]),
      CI_high_percentile = exp(qs[2]),
      P_boot_normal = 2 * pnorm(abs(beta0 / se), lower.tail = FALSE),
      valid_B = length(z)
    )
  }

  multinom_rows[[paste0(coh, "_mild")]] <- make_row(
    "Mild vs Independent", b0["mild"], boot[, "beta_mild"]
  )
  multinom_rows[[paste0(coh, "_severe")]] <- make_row(
    "Severe vs Independent", b0["severe"], boot[, "beta_severe"]
  )
  multinom_rows[[paste0(coh, "_contrast")]] <- make_row(
    "Multinomial Severe minus Mild", b0["contrast"], boot[, "beta_severe_minus_mild"]
  )
}

multinom_results <- rbindlist(multinom_rows, fill = TRUE)
boot_qc <- rbindlist(boot_qc_rows, fill = TRUE)
write_csv_utf8(multinom_results, file.path(TABLE_DIR, "02C_multinomial.csv"))
write_csv_utf8(boot_qc, file.path(QC_DIR, "02D_multinomial_bootstrap_QC.csv"))

# ---- C. Formal direct Severe vs Mild -----------------------------------------
direct_rows <- lapply(COHORTS, function(coh) {
  d <- build_estimand_data(D[[coh]], "SEVERE_VS_MILD")
  fit_cluster_logit(
    d, outcome = "y", exposure = "burden_z", covars = CORE0,
    cohort = coh, analysis_label = "SEVERE_VS_MILD", id_col = "id"
  )
})
direct_results <- rbindlist(direct_rows, fill = TRUE)
write_csv_utf8(direct_results, file.path(TABLE_DIR, "02E_direct_severe_vs_mild.csv"))

# ---- D. Four-cohort meta-analysis --------------------------------------------
meta_inputs <- list(
  "Any limitation vs Independent" = any_results[, .(cohort, beta, SE)],
  "Mild vs Independent" = multinom_results[estimand == "Mild vs Independent", .(
    cohort, beta, SE = SE_boot
  )],
  "Severe vs Independent" = multinom_results[estimand == "Severe vs Independent", .(
    cohort, beta, SE = SE_boot
  )],
  "Severe vs Mild" = direct_results[, .(cohort, beta, SE)]
)

meta_rows <- lapply(names(meta_inputs), function(label) {
  z <- meta_inputs[[label]]
  m <- meta_reml_mkh(z$beta, z$SE)
  m[, estimand := label]
  m[, positive_cohort_N := sum(exp(z$beta) > 1, na.rm = TRUE)]
  m
})
meta_results <- rbindlist(meta_rows, fill = TRUE)
setcolorder(meta_results, c("estimand", setdiff(names(meta_results), "estimand")))
write_csv_utf8(meta_results, file.path(TABLE_DIR, "02F_primary_meta.csv"))

# ---- E. Leave-one-cohort-out -------------------------------------------------
loo_rows <- list()
for (label in c("Any limitation vs Independent", "Severe vs Independent", "Severe vs Mild")) {
  z0 <- meta_inputs[[label]]
  full <- meta_reml_mkh(z0$beta, z0$SE)
  full[, `:=`(estimand = label, omitted = "None", analysis = "Primary pooled")]
  loo_rows[[length(loo_rows) + 1L]] <- full

  for (omit in COHORTS) {
    z <- z0[cohort != omit]
    m <- meta_reml_mkh(z$beta, z$SE)
    m[, `:=`(estimand = label, omitted = omit, analysis = "Leave-one-out")]
    loo_rows[[length(loo_rows) + 1L]] <- m
  }
}
loo <- rbindlist(loo_rows, fill = TRUE)
write_csv_utf8(loo, file.path(TABLE_DIR, "02G_leave_one_cohort_out.csv"))

# ---- F. Locked-result comparison (informational, never used in fitting) ------
expected_file <- file.path("metadata", "expected_primary_results.csv")
if (file.exists(expected_file)) {
  expected <- fread(expected_file)
  observed <- rbindlist(list(
    any_results[, .(cohort, estimand = "Any limitation vs Independent", OR = OR)],
    multinom_results[estimand == "Mild vs Independent", .(cohort, estimand, OR = effect)],
    multinom_results[estimand == "Severe vs Independent", .(cohort, estimand, OR = effect)],
    direct_results[, .(cohort, estimand = "Severe vs Mild", OR = OR)]
  ), fill = TRUE)

  qc_cohort <- merge(expected[type == "cohort"], observed, by = c("cohort", "estimand"), all.x = TRUE)
  qc_cohort[, abs_difference := abs(OR - expected_OR)]
  qc_cohort[, within_0_01 := abs_difference <= 0.01]

  pooled_obs <- meta_results[, .(type = "pooled", cohort = "ALL", estimand, OR = pooled_OR)]
  qc_pooled <- merge(expected[type == "pooled"], pooled_obs, by = c("type", "cohort", "estimand"), all.x = TRUE)
  qc_pooled[, abs_difference := abs(OR - expected_OR)]
  qc_pooled[, within_0_01 := abs_difference <= 0.01]

  qc_cohort[, type := "cohort"]
  qc <- rbindlist(list(qc_cohort, qc_pooled), fill = TRUE)
  write_csv_utf8(qc, file.path(QC_DIR, "02H_locked_result_comparison.csv"))
}

cat("\nPrimary analysis completed.\n")
print(meta_results[, .(estimand, pooled_OR, CI_low_mKH, CI_high_mKH, P_mKH, I2)])
