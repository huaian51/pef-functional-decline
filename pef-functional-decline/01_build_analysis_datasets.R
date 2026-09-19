# ============================================================================
# 01_build_analysis_datasets.R
# Public release for:
# "Cumulative low peak expiratory flow burden and severity of subsequent
# functional decline: a four-cohort longitudinal study"
#
# PURPOSE
# -------
# Build the unified repeated low-PEF exposure and cohort-specific analysis
# datasets from harmonized, access-controlled cohort extracts.
#
# IMPORTANT
# ---------
# This repository does not redistribute HRS, ELSA, SHARE, or CHARLS data.
# Users must obtain the source data and create the two harmonized input files
# per cohort described in data_harmonized/README.md.
# ============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

suppressPackageStartupMessages(library(data.table))
source(file.path("R", "helpers.R"))

COHORTS <- c("HRS", "ELSA", "SHARE", "CHARLS")
CORE0 <- c("age", "sex", "interval_f", "height", "bmi", "education")
LOW_Z_CUTOFF <- -1.0

INPUT_DIR <- file.path("data_harmonized")
OUT_DATA <- file.path("output", "derived_data")
OUT_QC <- file.path("output", "QC")
dir.create(OUT_DATA, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_QC, recursive = TRUE, showWarnings = FALSE)

PEF_REQUIRED <- c(
  "id",
  "pef_1", "age_1", "sex_1", "height_1",
  "pef_2", "age_2", "sex_2", "height_2"
)

INTERVAL_REQUIRED <- c(
  "id", "interval_f", "destination",
  "age", "sex", "height", "bmi", "education"
)

OPTIONAL_INTERVAL <- c(
  "origin_wave", "next_wave", "smoking", "current_smoking", "ever_smoking",
  "lung_disease", "diabetes", "hypertension", "heart_disease", "stroke",
  "cancer", "depression", "observed_next", "adl_iadl_t0", "adl_iadl_t1"
)

normalize_destination <- function(x) {
  y <- trimws(as.character(x))
  y[tolower(y) %in% c("independent", "independence")] <- "Independent"
  y[tolower(y) %in% c("mild", "mild limitation", "mild functional limitation")] <- "Mild limitation"
  y[tolower(y) %in% c("severe", "severe limitation", "severe functional limitation")] <- "Severe limitation"
  allowed <- c("Independent", "Mild limitation", "Severe limitation")
  bad <- unique(y[!is.na(y) & !(y %in% allowed)])
  if (length(bad)) stop("Unexpected destination value(s): ", paste(bad, collapse = ", "))
  y
}

build_unified_exposure <- function(pef_wide, cohort) {
  d <- as.data.table(copy(pef_wide))
  check_required(d, PEF_REQUIRED, paste0(cohort, " PEF input"))
  d[, id := normalize_id(id)]
  if (anyNA(d$id)) stop(cohort, ": missing participant ID in PEF input.")
  if (anyDuplicated(d$id)) stop(cohort, ": PEF input must contain one row per participant.")

  long1 <- d[, .(
    id,
    assessment = 1L,
    pef = suppressWarnings(as.numeric(pef_1)),
    age = suppressWarnings(as.numeric(age_1)),
    sex = as.factor(sex_1),
    height = suppressWarnings(as.numeric(height_1))
  )]
  long2 <- d[, .(
    id,
    assessment = 2L,
    pef = suppressWarnings(as.numeric(pef_2)),
    age = suppressWarnings(as.numeric(age_2)),
    sex = as.factor(sex_2),
    height = suppressWarnings(as.numeric(height_2))
  )]
  long <- rbindlist(list(long1, long2), use.names = TRUE)

  fit <- lm(pef ~ age + sex + height, data = as.data.frame(long), na.action = na.exclude)
  rr <- residuals(fit)
  resid_sd <- sd(rr, na.rm = TRUE)
  if (!is.finite(resid_sd) || resid_sd <= 0) stop(cohort, ": invalid residual SD in expected-PEF model.")

  long[, residual_z := as.numeric(rr) / resid_sd]
  wide_z <- dcast(long, id ~ assessment, value.var = "residual_z")
  setnames(wide_z, c("1", "2"), c("z1", "z2"))

  wide_z[, deficit1 := pmax(-z1, 0)]
  wide_z[, deficit2 := pmax(-z2, 0)]
  wide_z[, burden_raw_unified := rowMeans(cbind(deficit1, deficit2), na.rm = FALSE)]
  wide_z[, low1 := fifelse(!is.na(z1), as.integer(z1 <= LOW_Z_CUTOFF), NA_integer_)]
  wide_z[, low2 := fifelse(!is.na(z2), as.integer(z2 <= LOW_Z_CUTOFF), NA_integer_)]
  wide_z[, low_pef_count := fifelse(
    !is.na(low1) & !is.na(low2), low1 + low2, NA_integer_
  )]

  # Optional human-readable two-assessment state.
  wide_z[, pef_state := fifelse(
    low_pef_count == 0, "Stable preserved",
    fifelse(
      low_pef_count == 2, "Persistent low",
      fifelse(low1 == 1 & low2 == 0, "Low-to-preserved", "Incident low")
    )
  )]

  coef_dt <- data.table(
    cohort = cohort,
    term = names(coef(fit)),
    coefficient = as.numeric(coef(fit)),
    pooled_residual_sd = resid_sd,
    model_N = nobs(fit)
  )

  list(exposure = wide_z, model = fit, model_summary = coef_dt)
}

sample_qc <- list()
model_qc <- list()

for (coh in COHORTS) {
  message("\n==================== ", coh, " ====================")

  f_pef <- file.path(INPUT_DIR, paste0(coh, "_pef_exposure.rds"))
  f_int <- file.path(INPUT_DIR, paste0(coh, "_functional_intervals.rds"))
  if (!file.exists(f_pef) || !file.exists(f_int)) {
    stop(
      "Missing harmonized input(s) for ", coh, ":\n",
      f_pef, "\n", f_int,
      "\nSee data_harmonized/README.md.", call. = FALSE
    )
  }

  pef <- as.data.table(readRDS(f_pef))
  intervals <- as.data.table(readRDS(f_int))
  check_required(intervals, INTERVAL_REQUIRED, paste0(coh, " interval input"))

  exp_obj <- build_unified_exposure(pef, coh)
  e <- exp_obj$exposure
  model_qc[[coh]] <- exp_obj$model_summary

  d <- copy(intervals)
  d[, id := normalize_id(id)]
  d[, destination := normalize_destination(destination)]
  d[, interval_f := factor(interval_f)]
  d[, sex := factor(sex)]
  d[, age := suppressWarnings(as.numeric(age))]
  d[, height := suppressWarnings(as.numeric(height))]
  d[, bmi := suppressWarnings(as.numeric(bmi))]
  # Education is intentionally preserved in the harmonized class. In the
  # submitted analysis it is numeric in some cohorts and categorical in others.

  if (!("observed_next" %in% names(d))) {
    d[, observed_next := as.integer(destination %in% c(
      "Independent", "Mild limitation", "Severe limitation"
    ))]
  } else {
    d[, observed_next := as.integer(observed_next)]
  }

  d[, .row_order__ := .I]
  d <- merge(d, e, by = "id", all.x = TRUE, sort = FALSE)
  setorder(d, .row_order__)
  d[, .row_order__ := NULL]

  d[, dest3 := factor(
    destination,
    levels = c("Independent", "Mild limitation", "Severe limitation")
  )]
  d[, any_limitation := fifelse(
    destination %in% c("Mild limitation", "Severe limitation"), 1L,
    fifelse(destination == "Independent", 0L, NA_integer_)
  )]
  d[, severe_vs_nonsevere := fifelse(
    destination == "Severe limitation", 1L,
    fifelse(destination %in% c("Independent", "Mild limitation"), 0L, NA_integer_)
  )]
  d[, severe_vs_mild := fifelse(
    destination == "Severe limitation", 1L,
    fifelse(destination == "Mild limitation", 0L, NA_integer_)
  )]

  primary_needed <- c("id", "burden_raw_unified", "destination", CORE0)
  eligible_primary <- d$observed_next == 1L & complete.cases(d[, ..primary_needed])

  # Scale on unique participants in the locked primary sample, avoiding
  # interval-count weighting.
  person_burden <- unique(d[eligible_primary, .(id, burden_raw_unified)])
  if (anyDuplicated(person_burden$id)) {
    chk <- person_burden[, .(n_burden = uniqueN(burden_raw_unified)), by = id]
    if (any(chk$n_burden > 1L)) {
      stop(coh, ": burden_raw_unified is not constant within participant.")
    }
    person_burden <- person_burden[, .SD[1], by = id]
  }

  m <- mean(person_burden$burden_raw_unified, na.rm = TRUE)
  s <- sd(person_burden$burden_raw_unified, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) stop(coh, ": invalid primary-sample burden SD.")
  d[, burden_z := (burden_raw_unified - m) / s]
  d[, burden_z_unified := burden_z]

  primary_needed2 <- c("id", "burden_z", "destination", CORE0)
  eligible_primary <- d$observed_next == 1L & complete.cases(d[, ..primary_needed2])
  primary <- copy(d[eligible_primary])

  saveRDS(d, file.path(OUT_DATA, paste0("01_", coh, "_analysis_dataset.rds")), compress = "xz")
  saveRDS(primary, file.path(OUT_DATA, paste0("01_", coh, "_primary_analysis_dataset.rds")), compress = "xz")

  sample_qc[[coh]] <- data.table(
    cohort = coh,
    full_interval_N = nrow(d),
    primary_interval_N = nrow(primary),
    primary_person_N = uniqueN(primary$id),
    independent_N = sum(primary$destination == "Independent", na.rm = TRUE),
    mild_N = sum(primary$destination == "Mild limitation", na.rm = TRUE),
    severe_N = sum(primary$destination == "Severe limitation", na.rm = TRUE),
    burden_raw_mean_person = m,
    burden_raw_sd_person = s,
    exposure_coverage_primary_pct = 100 * mean(!is.na(primary$burden_z)),
    low0_N = sum(primary$low_pef_count == 0, na.rm = TRUE),
    low1_N = sum(primary$low_pef_count == 1, na.rm = TRUE),
    low2_N = sum(primary$low_pef_count == 2, na.rm = TRUE)
  )
}

write_csv_utf8(rbindlist(sample_qc, fill = TRUE), file.path(OUT_QC, "01_sample_QC.csv"))
write_csv_utf8(rbindlist(model_qc, fill = TRUE), file.path(OUT_QC, "01_expected_PEF_model_QC.csv"))

cat("\nPublic dataset-construction stage completed.\n")
cat("Exposure algorithm: pooled PEF ~ age + sex + height; no wave term;\n")
cat("deficit = max(0, -residual z); two-assessment mean; cohort-specific 1-SD scaling.\n")
