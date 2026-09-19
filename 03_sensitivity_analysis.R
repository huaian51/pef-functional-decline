# ============================================================================
# 03_sensitivity_analysis.R
# Exact public implementation of the sensitivity/supportive analyses retained
# in the manuscript.
#
# Manuscript:
# "Cumulative low peak expiratory flow burden and severity of subsequent
# functional decline: a four-cohort longitudinal study"
#
# This public script is aligned to the final manuscript sensitivity pathway.
# It deliberately differs from earlier reviewer-facing drafts in three places:
#   1) the 0/1/2 low-PEF analysis uses the full eligible interval risk set,
#      not the locked continuous-burden primary sample;
#   2) CHARLS ever-smoking is handled as a dedicated sensitivity covariate;
#   3) the strict >=2 ADL/IADL deterioration analysis uses track-specific
#      residual-z PEF inputs and re-standardizes the unified raw burden within
#      its own eligible person-track risk set.
#
# PUBLIC-DATA BOUNDARY
# --------------------
# No HRS, ELSA, SHARE, or CHARLS participant-level data are redistributed.
# This script reads harmonized inputs prepared by authorized users.
# ============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE, warn = 1, width = 380, scipen = 999)

suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required.", call. = FALSE)
}
source(file.path("R", "helpers.R"))

# ============================================================================
# 1. SETTINGS / INPUTS
# ============================================================================

COHORTS <- c("HRS", "ELSA", "SHARE", "CHARLS")
CORE0 <- c("age", "sex", "interval_f", "height", "bmi", "education")
ESTIMANDS <- c(
  "A_ANY_LIMITATION",
  "MILD_VS_INDEPENDENT",
  "SEVERE_VS_INDEPENDENT",
  "SEVERE_VS_MILD"
)

MIN_EXTENDED_COVERAGE <- 85
LOW_Z_CUTOFF <- -1.0

DATA_DIR <- file.path("output", "derived_data")
TABLE_DIR <- file.path("output", "tables")
QC_DIR <- file.path("output", "QC")
HARMONIZED_DIR <- "data_harmonized"

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

FULL_FILES <- setNames(
  file.path(DATA_DIR, paste0("01_", COHORTS, "_analysis_dataset.rds")),
  COHORTS
)

PRIMARY_FILES <- setNames(
  file.path(DATA_DIR, paste0("01_", COHORTS, "_primary_analysis_dataset.rds")),
  COHORTS
)

PRIMARY_META_FILE <- file.path(TABLE_DIR, "02F_primary_meta.csv")
PRIMARY_ANY_FILE <- file.path(TABLE_DIR, "02B_any_limitation.csv")
PRIMARY_MULTI_FILE <- file.path(TABLE_DIR, "02C_multinomial.csv")
PRIMARY_DIRECT_FILE <- file.path(TABLE_DIR, "02E_direct_severe_vs_mild.csv")

need_files <- c(FULL_FILES, PRIMARY_FILES, PRIMARY_META_FILE,
                PRIMARY_ANY_FILE, PRIMARY_MULTI_FILE, PRIMARY_DIRECT_FILE)
missing_files <- need_files[!file.exists(need_files)]
if (length(missing_files)) {
  stop(
    "Run 01_build_analysis_datasets.R and 02_primary_analysis.R first.\nMissing:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

Dfull <- setNames(
  lapply(FULL_FILES, function(f) as.data.table(readRDS(f))),
  COHORTS
)
Dprimary <- setNames(
  lapply(PRIMARY_FILES, function(f) as.data.table(readRDS(f))),
  COHORTS
)

for (coh in COHORTS) {
  if (!("cohort" %in% names(Dfull[[coh]]))) Dfull[[coh]][, cohort := coh]
  if (!("cohort" %in% names(Dprimary[[coh]]))) Dprimary[[coh]][, cohort := coh]

  Dfull[[coh]][, destination := as.character(destination)]
  Dprimary[[coh]][, destination := as.character(destination)]

  # Public harmonization convention:
  # use `smoking` as the current-smoking variable used in the manuscript.
  if (!("smoking" %in% names(Dfull[[coh]])) &&
      "current_smoking" %in% names(Dfull[[coh]])) {
    Dfull[[coh]][, smoking := current_smoking]
  }
  if (!("smoking" %in% names(Dprimary[[coh]])) &&
      "current_smoking" %in% names(Dprimary[[coh]])) {
    Dprimary[[coh]][, smoking := current_smoking]
  }
}

fit_all_estimands <- function(dat, covars, spec) {
  coh <- unique(dat$cohort)[1]
  out <- lapply(
    ESTIMANDS,
    function(est) {
      dd <- build_estimand_data(dat, est)
      z <- fit_cluster_logit(
        dd,
        outcome = "y",
        exposure = "burden_z",
        covars = covars,
        cohort = coh,
        analysis_label = est,
        id_col = "id"
      )
      z[, `:=`(spec = spec, estimand = est)]
      z
    }
  )
  rbindlist(out, fill = TRUE)
}

write_bom <- function(x, path) {
  fwrite(x, path, bom = TRUE, na = "")
  message("Saved: ", path)
}

# ============================================================================
# 2. FULL-INTERVAL EXPOSURE QC
# ============================================================================

full_exposure_qc <- rbindlist(
  lapply(COHORTS, function(coh) {
    d <- Dfull[[coh]]
    comparable <- d[!is.na(low_pef_count)]

    data.table(
      cohort = coh,
      interval_N = nrow(d),
      person_N = uniqueN(d$id),
      unified_burden_interval_coverage_pct = 100 * mean(!is.na(d$burden_z)),
      unified_burden_person_coverage_pct =
        100 * uniqueN(d$id[!is.na(d$burden_z)]) / uniqueN(d$id),
      low_count_comparable_interval_N = nrow(comparable),
      low_count_mismatch_interval_N = 0L
    )
  }),
  fill = TRUE
)

# ============================================================================
# 3. LEAVE-ONE-COHORT-OUT META-ANALYSIS
#    Same estimands retained in final manuscript sensitivity analysis.
# ============================================================================

A02 <- fread(PRIMARY_ANY_FILE)
C02 <- fread(PRIMARY_MULTI_FILE)
D02 <- fread(PRIMARY_DIRECT_FILE)

loo_source <- rbindlist(
  list(
    A02[, .(
      cohort,
      estimand = "A: Any limitation vs Independent",
      beta,
      SE
    )],
    C02[
      estimand == "Severe vs Independent",
      .(
        cohort,
        estimand = "C: Severe vs Independent",
        beta,
        SE = SE_boot
      )
    ],
    D02[, .(
      cohort,
      estimand = "Severity contrast: Severe vs Mild",
      beta,
      SE
    )]
  ),
  fill = TRUE
)

loo_rows <- list()
for (est in unique(loo_source$estimand)) {
  z0 <- loo_source[estimand == est]

  for (omit in COHORTS) {
    z <- z0[cohort != omit]
    mm <- meta_reml_mkh(z$beta, z$SE)
    mm[, `:=`(
      estimand = est,
      omitted_cohort = omit,
      cohort_N = uniqueN(z$cohort),
      positive_cohort_N = sum(z$beta > 0, na.rm = TRUE)
    )]
    loo_rows[[length(loo_rows) + 1L]] <- mm
  }
}
LOO <- rbindlist(loo_rows, fill = TRUE)

# ============================================================================
# 4. MODULE A — SMOKING / COHORT-SPECIFIC EXTENDED ADJUSTMENT
# ============================================================================

EXTENDED_CANDIDATES <- c(
  "smoking",
  "lung_disease",
  "diabetes",
  "hypertension",
  "heart_disease",
  "stroke",
  "cancer",
  "depression"
)

extended_manifest_rows <- list()
selected_extended <- list()

for (coh in COHORTS) {
  d <- copy(Dfull[[coh]])
  d <- d[destination %in% c("Independent", "Mild limitation", "Severe limitation")]

  selected <- character()

  for (v in EXTENDED_CANDIDATES) {
    if (!(v %in% names(d))) {
      extended_manifest_rows[[length(extended_manifest_rows) + 1L]] <-
        data.table(
          cohort = coh,
          construct = v,
          variable = NA_character_,
          coverage_pct = 0,
          usable = FALSE,
          included_in_extended = FALSE,
          reason = "NOT_FOUND"
        )
      next
    }

    covp <- coverage_pct(d[[v]])
    usable <- usable_var(d[[v]])
    include <- usable && covp >= MIN_EXTENDED_COVERAGE

    extended_manifest_rows[[length(extended_manifest_rows) + 1L]] <-
      data.table(
        cohort = coh,
        construct = v,
        variable = v,
        coverage_pct = covp,
        usable = usable,
        included_in_extended = include,
        reason = if (include) {
          "INCLUDED"
        } else if (!usable) {
          "UNUSABLE"
        } else {
          "COVERAGE_BELOW_THRESHOLD"
        }
      )

    if (include) selected <- c(selected, v)
  }

  selected_extended[[coh]] <- unique(selected)
}

extended_manifest <- rbindlist(extended_manifest_rows, fill = TRUE)

# Dedicated CHARLS ever-smoking sensitivity.
charls_ever_qc <- data.table(
  interval_N = nrow(Dfull[["CHARLS"]]),
  ever_smoking_nonmissing_N = NA_integer_,
  ever_smoking_coverage_pct = NA_real_,
  ever_smoker_pct = NA_real_
)

if ("ever_smoking" %in% names(Dfull[["CHARLS"]])) {
  dc <- Dfull[["CHARLS"]]
  dc[, charls_ever_smoking := suppressWarnings(as.numeric(ever_smoking))]
  dc[!charls_ever_smoking %in% c(0, 1), charls_ever_smoking := NA_real_]
  Dfull[["CHARLS"]] <- dc

  charls_ever_qc <- data.table(
    interval_N = nrow(dc),
    ever_smoking_nonmissing_N = sum(!is.na(dc$charls_ever_smoking)),
    ever_smoking_coverage_pct = 100 * mean(!is.na(dc$charls_ever_smoking)),
    ever_smoker_pct = 100 * mean(dc$charls_ever_smoking == 1, na.rm = TRUE)
  )
} else {
  warning(
    "CHARLS `ever_smoking` was not supplied. ",
    "The dedicated CHARLS ever-smoking sensitivity will be skipped."
  )
}

defensive_results_rows <- list()
defensive_sample_rows <- list()

for (coh in COHORTS) {
  d0 <- copy(Dfull[[coh]])
  d0 <- d0[destination %in% c("Independent", "Mild limitation", "Severe limitation")]

  need_core0 <- c("id", "burden_z", "destination", CORE0)
  d_core0 <- d0[complete.cases(d0[, ..need_core0])]

  z0 <- fit_all_estimands(d_core0, CORE0, "PRIMARY_CORE0")
  defensive_results_rows[[length(defensive_results_rows) + 1L]] <- z0

  defensive_sample_rows[[length(defensive_sample_rows) + 1L]] <-
    data.table(
      cohort = coh,
      spec = "PRIMARY_CORE0",
      interval_N = nrow(d_core0),
      person_N = uniqueN(d_core0$id),
      covariates = paste(CORE0, collapse = " + ")
    )

  if ("smoking" %in% names(d_core0) && usable_var(d_core0$smoking)) {
    d_smoke <- d_core0[!is.na(smoking)]

    z_sel <- fit_all_estimands(
      d_smoke,
      CORE0,
      "CORE0_ON_SMOKING_OBSERVED_SAMPLE"
    )
    z_smoke <- fit_all_estimands(
      d_smoke,
      c(CORE0, "smoking"),
      "CORE0_PLUS_CURRENT_SMOKING"
    )

    defensive_results_rows[[length(defensive_results_rows) + 1L]] <- z_sel
    defensive_results_rows[[length(defensive_results_rows) + 1L]] <- z_smoke

    defensive_sample_rows[[length(defensive_sample_rows) + 1L]] <-
      data.table(
        cohort = coh,
        spec = "CORE0_ON_SMOKING_OBSERVED_SAMPLE",
        interval_N = nrow(d_smoke),
        person_N = uniqueN(d_smoke$id),
        covariates = paste(CORE0, collapse = " + ")
      )

    defensive_sample_rows[[length(defensive_sample_rows) + 1L]] <-
      data.table(
        cohort = coh,
        spec = "CORE0_PLUS_CURRENT_SMOKING",
        interval_N = nrow(d_smoke),
        person_N = uniqueN(d_smoke$id),
        covariates = paste(c(CORE0, "smoking"), collapse = " + ")
      )
  }

  ext_covars <- unique(c(CORE0, selected_extended[[coh]]))
  need_ext <- c("id", "burden_z", "destination", ext_covars)
  d_ext <- d0[complete.cases(d0[, ..need_ext])]

  z_ext <- fit_all_estimands(
    d_ext,
    ext_covars,
    "COHORT_SPECIFIC_EXTENDED_RELIABLE"
  )
  defensive_results_rows[[length(defensive_results_rows) + 1L]] <- z_ext

  defensive_sample_rows[[length(defensive_sample_rows) + 1L]] <-
    data.table(
      cohort = coh,
      spec = "COHORT_SPECIFIC_EXTENDED_RELIABLE",
      interval_N = nrow(d_ext),
      person_N = uniqueN(d_ext$id),
      covariates = paste(ext_covars, collapse = " + ")
    )

  if (coh == "CHARLS" && "charls_ever_smoking" %in% names(d0) &&
      usable_var(d0$charls_ever_smoking)) {
    cov_ever <- c(CORE0, "charls_ever_smoking")
    need_ever <- c("id", "burden_z", "destination", cov_ever)
    d_ever <- d0[complete.cases(d0[, ..need_ever])]

    z_ever <- fit_all_estimands(
      d_ever,
      cov_ever,
      "CHARLS_CORE0_PLUS_EVER_SMOKING"
    )
    defensive_results_rows[[length(defensive_results_rows) + 1L]] <- z_ever

    defensive_sample_rows[[length(defensive_sample_rows) + 1L]] <-
      data.table(
        cohort = coh,
        spec = "CHARLS_CORE0_PLUS_EVER_SMOKING",
        interval_N = nrow(d_ever),
        person_N = uniqueN(d_ever$id),
        covariates = paste(cov_ever, collapse = " + ")
      )
  }
}

defensive_results <- rbindlist(defensive_results_rows, fill = TRUE)
defensive_samples <- rbindlist(defensive_sample_rows, fill = TRUE)

META_SPECS <- c(
  "PRIMARY_CORE0",
  "CORE0_ON_SMOKING_OBSERVED_SAMPLE",
  "CORE0_PLUS_CURRENT_SMOKING",
  "COHORT_SPECIFIC_EXTENDED_RELIABLE"
)

defensive_meta_rows <- list()
for (sp in META_SPECS) {
  for (est in ESTIMANDS) {
    z <- defensive_results[spec == sp & estimand == est]
    if (uniqueN(z$cohort) < 2L) next

    mm <- meta_reml_mkh(z$beta, z$SE)
    mm[, `:=`(
      spec = sp,
      estimand = est,
      cohort_N = uniqueN(z$cohort),
      positive_cohort_N = sum(z$OR > 1, na.rm = TRUE)
    )]
    defensive_meta_rows[[length(defensive_meta_rows) + 1L]] <- mm
  }
}
defensive_meta <- rbindlist(defensive_meta_rows, fill = TRUE)

# ============================================================================
# 5. MODULE B — 0 / 1 / 2 LOW-PEF FREQUENCY
#    IMPORTANT: uses full eligible interval data, not Dprimary.
# ============================================================================

fit_low_count_terms <- function(
  dat,
  formula_obj,
  terms_to_extract,
  cohort,
  estimand,
  exposure_spec
) {
  fit <- glm(
    formula_obj,
    data = as.data.frame(dat),
    family = binomial()
  )

  V <- cluster_vcov_glm(fit, dat$id)
  co <- coef(fit)
  out <- list()

  for (term in terms_to_extract) {
    if (!(term %in% names(co))) next

    b <- unname(co[term])
    se <- sqrt(pmax(V[term, term], 0))

    out[[length(out) + 1L]] <- data.table(
      cohort = cohort,
      estimand = estimand,
      exposure_spec = exposure_spec,
      term = term,
      interval_N = nrow(dat),
      person_N = uniqueN(dat$id),
      beta = b,
      SE = se,
      OR = exp(b),
      CI_low = exp(b - 1.96 * se),
      CI_high = exp(b + 1.96 * se),
      P = 2 * pnorm(abs(b / se), lower.tail = FALSE)
    )
  }

  rbindlist(out, fill = TRUE)
}

low_D <- list()
low_sample_rows <- list()
low_crosswalk_rows <- list()

for (coh in COHORTS) {
  d <- copy(Dfull[[coh]])
  d <- d[destination %in% c("Independent", "Mild limitation", "Severe limitation")]

  if (!("pef_state" %in% names(d))) {
    stop(coh, ": `pef_state` is required for the low-count crosswalk audit.")
  }
  if (!("low_pef_count" %in% names(d))) {
    stop(coh, ": `low_pef_count` is required.")
  }

  d[, pef_state_work := as.character(pef_state)]
  d[, low_count := suppressWarnings(as.numeric(low_pef_count))]

  need_low <- c(
    "id",
    "destination",
    "low_count",
    "pef_state_work",
    CORE0
  )
  d <- d[complete.cases(d[, ..need_low])]

  if (any(!d$low_count %in% c(0, 1, 2))) {
    stop(coh, ": low-count contains values outside 0/1/2.")
  }

  d[, low_count_cat := factor(
    low_count,
    levels = c(0, 1, 2),
    labels = c("0", "1", "2")
  )]

  d[, expected_low_count := fifelse(
    pef_state_work == "Stable preserved",
    0,
    fifelse(
      pef_state_work %in% c("Low-to-preserved", "Incident low"),
      1,
      fifelse(pef_state_work == "Persistent low", 2, NA_real_)
    )
  )]

  mismatch_N <- sum(d$low_count != d$expected_low_count, na.rm = TRUE)
  if (mismatch_N > 0L) {
    stop(coh, ": low-count / PEF-state crosswalk mismatch.")
  }

  low_sample_rows[[length(low_sample_rows) + 1L]] <-
    data.table(
      cohort = coh,
      interval_N = nrow(d),
      person_N = uniqueN(d$id),
      low0_N = sum(d$low_count == 0),
      low1_N = sum(d$low_count == 1),
      low2_N = sum(d$low_count == 2),
      crosswalk_mismatch_N = mismatch_N
    )

  low_crosswalk_rows[[length(low_crosswalk_rows) + 1L]] <-
    d[, .N, by = .(pef_state_work, low_count, expected_low_count)][
      , cohort := coh
    ]

  low_D[[coh]] <- d
}

low_sample_qc <- rbindlist(low_sample_rows, fill = TRUE)
low_crosswalk_qc <- rbindlist(low_crosswalk_rows, fill = TRUE)

# Locked manuscript low-count sample counts: audit only.
EXPECTED_LOW_SAMPLE <- data.table(
  cohort = COHORTS,
  interval_N_expected = c(12117L, 11054L, 30564L, 8355L),
  person_N_expected = c(3875L, 2907L, 11621L, 4920L),
  low0_N_expected = c(10546L, 9075L, 25970L, 6512L),
  low1_N_expected = c(1021L, 1401L, 3613L, 1446L),
  low2_N_expected = c(550L, 578L, 981L, 397L)
)
low_sample_audit <- merge(
  EXPECTED_LOW_SAMPLE,
  low_sample_qc,
  by = "cohort",
  all.x = TRUE
)
low_sample_audit[, exact_match :=
  interval_N_expected == interval_N &
  person_N_expected == person_N &
  low0_N_expected == low0_N &
  low1_N_expected == low1_N &
  low2_N_expected == low2_N
]

low_results_rows <- list()

for (coh in COHORTS) {
  d0 <- low_D[[coh]]

  for (est in ESTIMANDS) {
    dd <- build_estimand_data(d0, est)

    f_cat <- as.formula(
      paste("y ~ low_count_cat +", paste(CORE0, collapse = " + "))
    )
    z_cat <- fit_low_count_terms(
      dd,
      formula_obj = f_cat,
      terms_to_extract = c("low_count_cat1", "low_count_cat2"),
      cohort = coh,
      estimand = est,
      exposure_spec = "CATEGORICAL_0_1_2"
    )
    z_cat[term == "low_count_cat1", contrast := "1 vs 0"]
    z_cat[term == "low_count_cat2", contrast := "2 vs 0"]
    low_results_rows[[length(low_results_rows) + 1L]] <- z_cat

    f_trend <- as.formula(
      paste("y ~ low_count +", paste(CORE0, collapse = " + "))
    )
    z_trend <- fit_low_count_terms(
      dd,
      formula_obj = f_trend,
      terms_to_extract = "low_count",
      cohort = coh,
      estimand = est,
      exposure_spec = "TREND_PER_1_LOW_OCCASION"
    )
    z_trend[, contrast := "per +1 low occasion"]
    low_results_rows[[length(low_results_rows) + 1L]] <- z_trend
  }
}

low_results <- rbindlist(low_results_rows, fill = TRUE)

low_meta_rows <- list()
meta_groups <- unique(
  low_results[, .(estimand, exposure_spec, contrast)]
)

for (i in seq_len(nrow(meta_groups))) {
  g <- meta_groups[i]
  z <- low_results[
    estimand == g$estimand &
      exposure_spec == g$exposure_spec &
      contrast == g$contrast
  ]

  mm <- meta_reml_mkh(z$beta, z$SE)
  mm[, `:=`(
    estimand = g$estimand,
    exposure_spec = g$exposure_spec,
    contrast = g$contrast,
    cohort_N = uniqueN(z$cohort),
    positive_cohort_N = sum(z$OR > 1, na.rm = TRUE)
  )]
  low_meta_rows[[length(low_meta_rows) + 1L]] <- mm
}

low_meta <- rbindlist(low_meta_rows, fill = TRUE)

# ============================================================================
# 6. MODULE C — SEVERE VS ALL NON-SEVERE
# ============================================================================

severe_results_rows <- list()

for (coh in COHORTS) {
  d <- copy(Dfull[[coh]])
  d <- d[destination %in% c("Independent", "Mild limitation", "Severe limitation")]

  need <- c("id", "burden_z", "low_pef_count", "destination", CORE0)
  d <- d[complete.cases(d[, ..need])]

  d[, severe_any := as.integer(destination == "Severe limitation")]
  d[, low_count := as.numeric(low_pef_count)]

  z1 <- fit_cluster_logit(
    d,
    outcome = "severe_any",
    exposure = "burden_z",
    covars = CORE0,
    cohort = coh,
    analysis_label = "Continuous cumulative burden per +1 SD",
    id_col = "id"
  )
  z2 <- fit_cluster_logit(
    d,
    outcome = "severe_any",
    exposure = "low_count",
    covars = CORE0,
    cohort = coh,
    analysis_label = "Per +1 low-PEF occasion",
    id_col = "id"
  )

  severe_results_rows[[length(severe_results_rows) + 1L]] <- z1
  severe_results_rows[[length(severe_results_rows) + 1L]] <- z2
}

severe_results <- rbindlist(severe_results_rows, fill = TRUE)

severe_meta_rows <- list()
for (ex in unique(severe_results$analysis)) {
  z <- severe_results[analysis == ex]
  mm <- meta_reml_mkh(z$beta, z$SE)
  mm[, `:=`(
    exposure = ex,
    cohort_N = uniqueN(z$cohort),
    positive_cohort_N = sum(z$OR > 1, na.rm = TRUE)
  )]
  severe_meta_rows[[length(severe_meta_rows) + 1L]] <- mm
}
severe_meta <- rbindlist(severe_meta_rows, fill = TRUE)

# ============================================================================
# 7. MODULE D — OBSERVED-FOLLOW-UP IPW
# ============================================================================

selection_qc_rows <- list()
ipw_results_rows <- list()

for (coh in COHORTS) {
  d <- copy(Dfull[[coh]])

  d[, id_analysis := as.character(id)]
  d[, destination_chr := as.character(destination)]

  if (!("observed_next" %in% names(d))) {
    d[, observed_next := as.integer(
      destination_chr %in% c("Independent", "Mild limitation", "Severe limitation")
    )]
  }
  d[, destination_observed_audit := as.integer(observed_next)]

  selection_need <- c(
    "id_analysis",
    "burden_z",
    CORE0,
    "destination_observed_audit"
  )
  ds <- d[complete.cases(d[, ..selection_need])]

  sel_fit <- glm(
    as.formula(
      paste(
        "destination_observed_audit ~ burden_z +",
        paste(CORE0, collapse = " + ")
      )
    ),
    data = as.data.frame(ds),
    family = binomial()
  )

  sel_v <- cluster_vcov_glm(sel_fit, ds$id_analysis)
  sel_b <- unname(coef(sel_fit)["burden_z"])
  sel_se <- sqrt(pmax(sel_v["burden_z", "burden_z"], 0))

  ps <- pmin(pmax(fitted(sel_fit), 0.01), 0.99)
  ds[, p_observed := ps]

  p_obs_marg <- mean(ds$destination_observed_audit == 1)
  ds[, sw_observed := fifelse(
    destination_observed_audit == 1,
    p_obs_marg / p_observed,
    NA_real_
  )]

  selection_qc_rows[[length(selection_qc_rows) + 1L]] <-
    data.table(
      cohort = coh,
      interval_N = nrow(ds),
      person_N = uniqueN(ds$id_analysis),
      observed_N = sum(ds$destination_observed_audit == 1),
      observed_pct = 100 * mean(ds$destination_observed_audit == 1),
      selection_OR_per_SD_burden = exp(sel_b),
      selection_CI_low = exp(sel_b - 1.96 * sel_se),
      selection_CI_high = exp(sel_b + 1.96 * sel_se),
      selection_P = 2 * pnorm(abs(sel_b / sel_se), lower.tail = FALSE),
      SW_mean = mean(
        ds$sw_observed[ds$destination_observed_audit == 1],
        na.rm = TRUE
      ),
      SW_p01 = as.numeric(
        quantile(
          ds$sw_observed[ds$destination_observed_audit == 1],
          0.01,
          na.rm = TRUE
        )
      ),
      SW_p99 = as.numeric(
        quantile(
          ds$sw_observed[ds$destination_observed_audit == 1],
          0.99,
          na.rm = TRUE
        )
      )
    )

  dobs <- ds[
    destination_observed_audit == 1 &
      destination_chr %in% c("Independent", "Mild limitation", "Severe limitation")
  ]

  dobs[, y_any := as.integer(destination_chr != "Independent")]

  ipw_results_rows[[length(ipw_results_rows) + 1L]] <-
    fit_cluster_logit(
      dobs,
      outcome = "y_any",
      exposure = "burden_z",
      covars = CORE0,
      cohort = coh,
      analysis_label = "OBSERVED_ONLY_ANY_LIMITATION_UNWEIGHTED",
      id_col = "id_analysis"
    )

  ipw_results_rows[[length(ipw_results_rows) + 1L]] <-
    fit_cluster_logit(
      dobs,
      outcome = "y_any",
      exposure = "burden_z",
      covars = CORE0,
      cohort = coh,
      analysis_label = "IPW_FOR_OBSERVATION_ANY_LIMITATION",
      id_col = "id_analysis",
      weights_vec = dobs$sw_observed
    )

  dobs[, y_severe_any := as.integer(destination_chr == "Severe limitation")]

  ipw_results_rows[[length(ipw_results_rows) + 1L]] <-
    fit_cluster_logit(
      dobs,
      outcome = "y_severe_any",
      exposure = "burden_z",
      covars = CORE0,
      cohort = coh,
      analysis_label = "OBSERVED_ONLY_SEVERE_VS_NONSEVERE_UNWEIGHTED",
      id_col = "id_analysis"
    )

  ipw_results_rows[[length(ipw_results_rows) + 1L]] <-
    fit_cluster_logit(
      dobs,
      outcome = "y_severe_any",
      exposure = "burden_z",
      covars = CORE0,
      cohort = coh,
      analysis_label = "IPW_FOR_OBSERVATION_SEVERE_VS_NONSEVERE",
      id_col = "id_analysis",
      weights_vec = dobs$sw_observed
    )

  dsm <- dobs[
    destination_chr %in% c("Mild limitation", "Severe limitation")
  ]
  dsm[, y_severe_mild := as.integer(destination_chr == "Severe limitation")]

  ipw_results_rows[[length(ipw_results_rows) + 1L]] <-
    fit_cluster_logit(
      dsm,
      outcome = "y_severe_mild",
      exposure = "burden_z",
      covars = CORE0,
      cohort = coh,
      analysis_label = "OBSERVED_ONLY_SEVERE_VS_MILD_UNWEIGHTED",
      id_col = "id_analysis"
    )

  ipw_results_rows[[length(ipw_results_rows) + 1L]] <-
    fit_cluster_logit(
      dsm,
      outcome = "y_severe_mild",
      exposure = "burden_z",
      covars = CORE0,
      cohort = coh,
      analysis_label = "IPW_FOR_OBSERVATION_SEVERE_VS_MILD",
      id_col = "id_analysis",
      weights_vec = dsm$sw_observed
    )
}

selection_qc <- rbindlist(selection_qc_rows, fill = TRUE)
ipw_results <- rbindlist(ipw_results_rows, fill = TRUE)
ipw_meta <- run_meta_by_analysis(ipw_results)

# ============================================================================
# 8. MODULE E — STRICT >=2 ADL/IADL-ITEM DETERIORATION
#
# Exact public interface:
#   data_harmonized/<COHORT>_ge2_sensitivity.rds
#
# Required standardized columns:
#   id, track, z1, z2,
#   baseline_adl_count, baseline_iadl_count, adl_count, iadl_count
#
# The z1/z2 values must be the residual-z values from the SAME cohort-specific
# expected-PEF models used for the corresponding exposure track.
#
# Recommended covariates:
#   age, sex, height, bmi, education, interval_f
#
# Optional:
#   include_worsening_only
#
# HRS and SHARE must retain their track/window labels because their strict
# sensitivity uses track-specific PEF exposure windows. ELSA/CHARLS may use
# any constant track label.
# ============================================================================

GE2_FILES <- setNames(
  file.path(HARMONIZED_DIR, paste0(COHORTS, "_ge2_sensitivity.rds")),
  COHORTS
)

ge2_track_qc_rows <- list()
ge2_manifest_rows <- list()
ge2_results_rows <- list()

missing_ge2 <- GE2_FILES[!file.exists(GE2_FILES)]
if (length(missing_ge2)) {
  warning(
    "Strict >=2 deterioration input(s) are missing; that module will be skipped:\n",
    paste(missing_ge2, collapse = "\n")
  )
} else {
  GE2_COVAR_CANDIDATES <- c(
    "age", "sex", "height", "bmi", "education", "interval_f"
  )

  for (coh in COHORTS) {
    d <- as.data.table(readRDS(GE2_FILES[[coh]]))

    required <- c(
      "id", "track", "z1", "z2",
      "baseline_adl_count", "baseline_iadl_count",
      "adl_count", "iadl_count"
    )
    check_required(d, required, paste0(coh, " strict >=2 sensitivity input"))

    d[, id_key := normalize_id(id)]
    d[, track_key := trimws(as.character(track))]
    d[, z1 := suppressWarnings(as.numeric(z1))]
    d[, z2 := suppressWarnings(as.numeric(z2))]

    if (coh %in% c("HRS", "SHARE") && any(is.na(d$track_key) | d$track_key == "")) {
      stop(coh, ": non-missing track labels are required for the strict >=2 analysis.")
    }

    if (coh == "HRS" && !all(unique(d$track_key) %in% c("HRS_A", "HRS_B"))) {
      stop("HRS strict-outcome track labels must be HRS_A/HRS_B.")
    }
    if (coh == "SHARE" && !all(unique(d$track_key) %in% c("A", "B"))) {
      stop("SHARE strict-outcome track labels must be A/B.")
    }

    d[, deficit1 := pmax(-z1, 0)]
    d[, deficit2 := pmax(-z2, 0)]
    d[, burden_raw_unified := rowMeans(cbind(deficit1, deficit2), na.rm = FALSE)]

    ge2_track_qc_rows[[length(ge2_track_qc_rows) + 1L]] <-
      data.table(
        cohort = coh,
        endpoint_rows_N = nrow(d),
        exposure_matched_N = sum(!is.na(d$burden_raw_unified)),
        exposure_coverage_pct = 100 * mean(!is.na(d$burden_raw_unified)),
        track_var = "track",
        track_values = paste(sort(unique(d$track_key)), collapse = "; ")
      )

    d[, baseline_total :=
        suppressWarnings(as.numeric(baseline_adl_count)) +
        suppressWarnings(as.numeric(baseline_iadl_count))
    ]
    d[, follow_total :=
        suppressWarnings(as.numeric(adl_count)) +
        suppressWarnings(as.numeric(iadl_count))
    ]
    d[, delta_total := follow_total - baseline_total]
    d[, ge2_deterioration := fifelse(
      !is.na(delta_total),
      as.integer(delta_total >= 2),
      NA_integer_
    )]

    d0 <- d[baseline_total == 0]

    if ("include_worsening_only" %in% names(d0) &&
        is.logical(d0$include_worsening_only)) {
      d0 <- d0[
        is.na(include_worsening_only) |
          include_worsening_only == TRUE
      ]
    }

    # Same scaling used in the locked final analysis:
    # unique person-track units within the eligible strict-outcome risk set.
    scale_units <- unique(
      d0[
        !is.na(burden_raw_unified),
        .(id_key, track_key, burden_raw_unified)
      ]
    )
    m <- mean(scale_units$burden_raw_unified, na.rm = TRUE)
    s <- sd(scale_units$burden_raw_unified, na.rm = TRUE)

    if (!is.finite(s) || s <= 0) {
      stop(coh, ": invalid unified burden SD in strict >=2 risk set.")
    }
    d0[, burden_z := (burden_raw_unified - m) / s]

    present_covars <- GE2_COVAR_CANDIDATES[
      GE2_COVAR_CANDIDATES %in% names(d0)
    ]
    covars_model <- present_covars[
      vapply(
        present_covars,
        function(v) usable_var(d0[[v]]),
        logical(1)
      )
    ]

    # Preserve model classes supplied by harmonization.
    d0[, id_analysis := as.character(id_key)]

    need_model <- c(
      "id_analysis",
      "burden_z",
      "ge2_deterioration",
      covars_model
    )

    dm <- d0[
      ge2_deterioration %in% c(0, 1) &
        complete.cases(d0[, ..need_model])
    ]

    ge2_manifest_rows[[length(ge2_manifest_rows) + 1L]] <-
      data.table(
        cohort = coh,
        endpoint_source_file = basename(GE2_FILES[[coh]]),
        historical_burden_used = FALSE,
        track_var = "track",
        baseline_free_interval_N = nrow(d0),
        model_interval_N = nrow(dm),
        model_person_N = uniqueN(dm$id_analysis),
        ge2_event_N = sum(dm$ge2_deterioration == 1),
        ge2_event_pct = 100 * mean(dm$ge2_deterioration == 1),
        unified_raw_mean_unique_person_track = m,
        unified_raw_sd_unique_person_track = s,
        covariates = paste(covars_model, collapse = " + ")
      )

    dm[, y_ge2 := ge2_deterioration]

    z <- fit_cluster_logit(
      dm,
      outcome = "y_ge2",
      exposure = "burden_z",
      covars = covars_model,
      cohort = coh,
      analysis_label = "STRICT_GE2_ITEM_DETERIORATION",
      id_col = "id_analysis"
    )

    ge2_results_rows[[length(ge2_results_rows) + 1L]] <- z
  }
}

if (length(ge2_results_rows)) {
  ge2_track_qc <- rbindlist(ge2_track_qc_rows, fill = TRUE)
  ge2_manifest <- rbindlist(ge2_manifest_rows, fill = TRUE)
  ge2_results <- rbindlist(ge2_results_rows, fill = TRUE)
  ge2_meta <- run_meta_by_analysis(ge2_results)
} else {
  ge2_track_qc <- data.table()
  ge2_manifest <- data.table()
  ge2_results <- data.table()
  ge2_meta <- data.table()
}

# ============================================================================
# 9. FINAL PRACTICAL GATE + LOCKED RESULT AUDIT
# ============================================================================

key_def <- defensive_meta[
  spec %in% c(
    "CORE0_PLUS_CURRENT_SMOKING",
    "COHORT_SPECIFIC_EXTENDED_RELIABLE"
  ) &
    estimand %in% c(
      "A_ANY_LIMITATION",
      "SEVERE_VS_INDEPENDENT",
      "SEVERE_VS_MILD"
    )
]

key_ipw <- ipw_meta[
  analysis %in% c(
    "IPW_FOR_OBSERVATION_ANY_LIMITATION",
    "IPW_FOR_OBSERVATION_SEVERE_VS_NONSEVERE",
    "IPW_FOR_OBSERVATION_SEVERE_VS_MILD"
  )
]

key_severe <- severe_meta[
  exposure == "Continuous cumulative burden per +1 SD"
]

ge2_pool_or <- if (nrow(ge2_meta)) ge2_meta$pooled_OR[1] else NA_real_

FINAL_GATE <- data.table(
  check = c(
    "Unified burden coverage >=99% in full current interval datasets",
    "Unified residual-z low-count has zero mismatch where comparable",
    "All key smoking/extended pooled ORs remain >1",
    "Continuous Severe-vs-nonsevere pooled OR remains >1",
    "All key IPW pooled ORs remain >1",
    "Strict >=2-item pooled OR remains >1",
    "All leave-one-out pooled ORs remain >1",
    "Locked low-count sample counts reproduced"
  ),
  pass = c(
    all(full_exposure_qc$unified_burden_person_coverage_pct >= 99),
    all(full_exposure_qc$low_count_mismatch_interval_N == 0),
    nrow(key_def) > 0 && all(key_def$pooled_OR > 1, na.rm = TRUE),
    nrow(key_severe) == 1L && key_severe$pooled_OR[1] > 1,
    nrow(key_ipw) > 0 && all(key_ipw$pooled_OR > 1, na.rm = TRUE),
    is.finite(ge2_pool_or) && ge2_pool_or > 1,
    nrow(LOO) > 0 && all(LOO$pooled_OR > 1, na.rm = TRUE),
    all(low_sample_audit$exact_match)
  ),
  value = c(
    paste0(
      sprintf(
        "%.2f",
        min(full_exposure_qc$unified_burden_person_coverage_pct, na.rm = TRUE)
      ),
      "% minimum person coverage"
    ),
    "0 mismatched comparable intervals",
    paste0("minimum key pooled OR=", sprintf("%.3f", min(key_def$pooled_OR, na.rm = TRUE))),
    paste0("pooled OR=", sprintf("%.3f", key_severe$pooled_OR[1])),
    paste0("minimum IPW pooled OR=", sprintf("%.3f", min(key_ipw$pooled_OR, na.rm = TRUE))),
    if (is.finite(ge2_pool_or)) {
      paste0("pooled OR=", sprintf("%.3f", ge2_pool_or))
    } else {
      "strict >=2 module not run"
    },
    paste0("minimum LOO pooled OR=", sprintf("%.3f", min(LOO$pooled_OR, na.rm = TRUE))),
    paste0(sum(low_sample_audit$exact_match), "/4 cohorts exact")
  )
)

FINAL_DECISION <- if (all(FINAL_GATE$pass)) {
  "PASS_UNIFIED_SENSITIVITY"
} else {
  "REVIEW_REQUIRED"
}

# Informational comparison against the locked manuscript pooled sensitivity
# results. These targets NEVER enter model fitting.
expected_file <- file.path("metadata", "expected_sensitivity_results.csv")
locked_audit <- data.table()

if (file.exists(expected_file)) {
  expected <- fread(expected_file)

  observed <- rbindlist(
    list(
      defensive_meta[
        ,
        .(
          module = "smoking_extended",
          analysis = paste(spec, estimand, sep = " | "),
          observed_OR = pooled_OR,
          observed_CI_low = CI_low_mKH,
          observed_CI_high = CI_high_mKH,
          observed_I2 = I2
        )
      ],
      low_meta[
        ,
        .(
          module = "low_count",
          analysis = paste(estimand, contrast, sep = " | "),
          observed_OR = pooled_OR,
          observed_CI_low = CI_low_mKH,
          observed_CI_high = CI_high_mKH,
          observed_I2 = I2
        )
      ],
      severe_meta[
        exposure == "Continuous cumulative burden per +1 SD",
        .(
          module = "alternative_outcome",
          analysis = "SEVERE_VS_ALL_NONSEVERE",
          observed_OR = pooled_OR,
          observed_CI_low = CI_low_mKH,
          observed_CI_high = CI_high_mKH,
          observed_I2 = I2
        )
      ],
      ipw_meta[
        ,
        .(
          module = "ipw",
          analysis,
          observed_OR = pooled_OR,
          observed_CI_low = CI_low_mKH,
          observed_CI_high = CI_high_mKH,
          observed_I2 = I2
        )
      ],
      if (nrow(ge2_meta)) {
        ge2_meta[
          ,
          .(
            module = "strict_outcome",
            analysis = "STRICT_GE2_ADL_IADL_DETERIORATION",
            observed_OR = pooled_OR,
            observed_CI_low = CI_low_mKH,
            observed_CI_high = CI_high_mKH,
            observed_I2 = I2
          )
        ]
      } else {
        data.table()
      }
    ),
    fill = TRUE
  )

  locked_audit <- merge(
    expected,
    observed,
    by = c("module", "analysis"),
    all.x = TRUE
  )
  locked_audit[, OR_abs_diff := abs(observed_OR - expected_OR)]
}

# ============================================================================
# 10. SAVE — FILENAMES MATCH THE LOCKED FINAL MANUSCRIPT PIPELINE
# ============================================================================

write_bom(full_exposure_qc, file.path(QC_DIR, "31A_full_interval_unified_exposure_QC.csv"))
write_bom(extended_manifest, file.path(QC_DIR, "31B_extended_covariate_manifest.csv"))
write_bom(charls_ever_qc, file.path(QC_DIR, "31C_CHARLS_ever_smoking_QC.csv"))
write_bom(defensive_samples, file.path(QC_DIR, "31D_defensive_sample_manifest.csv"))
write_bom(defensive_results, file.path(TABLE_DIR, "31E_smoking_extended_adjustment_results.csv"))
write_bom(defensive_meta, file.path(TABLE_DIR, "31F_smoking_extended_adjustment_meta.csv"))
write_bom(low_sample_qc, file.path(QC_DIR, "31G_low_count_sample_QC.csv"))
write_bom(low_crosswalk_qc, file.path(QC_DIR, "31H_low_count_crosswalk_QC.csv"))
write_bom(low_sample_audit, file.path(QC_DIR, "31H2_low_count_locked_sample_audit.csv"))
write_bom(low_results, file.path(TABLE_DIR, "31I_low_count_supportive_results.csv"))
write_bom(low_meta, file.path(TABLE_DIR, "31J_low_count_supportive_meta.csv"))
write_bom(severe_results, file.path(TABLE_DIR, "31K_severe_vs_nonsevere_results.csv"))
write_bom(severe_meta, file.path(TABLE_DIR, "31L_severe_vs_nonsevere_meta.csv"))
write_bom(selection_qc, file.path(QC_DIR, "31M_selection_model_IPW_QC.csv"))
write_bom(ipw_results, file.path(TABLE_DIR, "31N_IPW_observed_followup_results.csv"))
write_bom(ipw_meta, file.path(TABLE_DIR, "31O_IPW_observed_followup_meta.csv"))

if (nrow(ge2_track_qc)) write_bom(ge2_track_qc, file.path(QC_DIR, "31P_ge2_unified_exposure_track_QC.csv"))
if (nrow(ge2_manifest)) write_bom(ge2_manifest, file.path(QC_DIR, "31Q_ge2_sample_manifest.csv"))
if (nrow(ge2_results)) write_bom(ge2_results, file.path(TABLE_DIR, "31R_ge2_deterioration_results.csv"))
if (nrow(ge2_meta)) write_bom(ge2_meta, file.path(TABLE_DIR, "31S_ge2_deterioration_meta.csv"))

write_bom(LOO, file.path(TABLE_DIR, "31T_leave_one_cohort_out_meta.csv"))
write_bom(FINAL_GATE, file.path(QC_DIR, "31U_FINAL_UNIFIED_SENSITIVITY_GATE.csv"))
if (nrow(locked_audit)) write_bom(locked_audit, file.path(QC_DIR, "31V_locked_sensitivity_result_audit.csv"))

capture.output(
  sessionInfo(),
  file = file.path(QC_DIR, "31_sessionInfo.txt")
)

cat("\n")
cat("===============================================================================\n")
cat("PUBLIC SENSITIVITY ANALYSIS COMPLETE\n")
cat("===============================================================================\n")
cat("Final decision: ", FINAL_DECISION, "\n", sep = "")
cat("Low-count sample exact matches: ", sum(low_sample_audit$exact_match), "/4\n", sep = "")
if (!nrow(ge2_meta)) {
  cat("NOTE: strict >=2 analysis was skipped because its dedicated harmonized\n")
  cat("      input files were not all supplied.\n")
}
cat("===============================================================================\n")
