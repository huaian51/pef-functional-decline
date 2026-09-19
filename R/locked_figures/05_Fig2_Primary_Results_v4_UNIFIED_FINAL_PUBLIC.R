# ============================================================================== 
# 05_Fig2_Primary_Results_FINAL_v3.R
#
# FIGURE 2 -- PRIMARY UNIFIED CONTINUOUS-BURDEN RESULTS (FINAL V4)
# ============================================================================== 

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  warn = 1,
  width = 300,
  scipen = 999
)

required_pkgs <- c("data.table", "ggplot2")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1. PATHS
# ------------------------------------------------------------------------------
PROJECT_ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)
UNIFIED_ROOT <- file.path(PROJECT_ROOT, "output", "compat", "unified_burden_check")
TABLE_DIR <- file.path(UNIFIED_ROOT, "tables")
QC_DIR <- file.path(UNIFIED_ROOT, "QC")
FIG_DIR <- file.path(PROJECT_ROOT, "output", "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

F_A    <- file.path(TABLE_DIR, "30C_any_limitation_unified_burden.csv")
F_C    <- file.path(TABLE_DIR, "30D_multinomial_unified_burden.csv")
F_D    <- file.path(TABLE_DIR, "30F_direct_severe_vs_mild_unified_burden.csv")
F_META <- file.path(TABLE_DIR, "30G_meta_unified_burden.csv")
F_GATE <- file.path(QC_DIR, "30J_unified_burden_stability_gate.csv")

INPUT_FILES <- c(F_A, F_C, F_D, F_META, F_GATE)
missing_files <- INPUT_FILES[!file.exists(INPUT_FILES)]
if (length(missing_files) > 0L) {
  stop(
    paste0(
      "Missing Figure 2 input file(s):\n",
      paste(missing_files, collapse = "\n"),
      "\n\nRun Script 30 (30_UNIFIED_CONTINUOUS_PEF_BURDEN_REANALYSIS.R) first."
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 2. HELPERS / ORDERS
# ------------------------------------------------------------------------------
COHORT_ORDER <- c("HRS", "ELSA", "SHARE", "CHARLS", "Pooled")

PANEL_MAP <- data.table(
  estimand = c(
    "Any limitation vs Independent",
    "Mild limitation vs Independent",
    "Severe limitation vs Independent",
    "Severe vs Mild"
  ),
  meta_estimand = c(
    "A: Any limitation vs Independent",
    "C: Mild vs Independent",
    "C: Severe vs Independent",
    "Severity contrast: Severe vs Mild"
  ),
  panel_letter = c("A", "B", "C", "D"),
  panel_title = c(
    "Any limitation vs Independent",
    "Mild limitation vs Independent",
    "Severe limitation vs Independent",
    "Severe vs Mild"
  )
)

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_or_ci <- function(est, low, high, digits = 2) {
  paste0(fmt_num(est, digits), " (", fmt_num(low, digits), "-", fmt_num(high, digits), ")")
}

# ------------------------------------------------------------------------------
# 3. READ INPUTS
# ------------------------------------------------------------------------------
A <- fread(F_A, encoding = "UTF-8")
C <- fread(F_C, encoding = "UTF-8")
D <- fread(F_D, encoding = "UTF-8")
META <- fread(F_META, encoding = "UTF-8")
GATE <- fread(F_GATE, encoding = "UTF-8")

if (!all(c("check", "pass", "value") %in% names(GATE))) {
  stop("30J unified-burden gate has unexpected columns.", call. = FALSE)
}

gate_pass <- if (is.logical(GATE$pass)) {
  GATE$pass
} else {
  tolower(trimws(as.character(GATE$pass))) %in% c("true", "t", "1", "yes")
}

if (!all(gate_pass)) {
  cat("\nScript 30 stability gate:\n")
  print(GATE)
  stop(
    "Script 30 unified-burden stability gate did not fully pass. Figure 2 generation stopped.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 3B. SCRIPT 30 INPUT QC
# ------------------------------------------------------------------------------
EXPECTED_COHORTS <- c("HRS", "ELSA", "SHARE", "CHARLS")

need_A <- c("cohort", "effect", "CI_low", "CI_high", "P", "interval_N", "person_N")
need_C <- c(
  "cohort", "estimand", "effect",
  "CI_low_percentile", "CI_high_percentile",
  "P_boot_normal", "interval_N", "person_N"
)
need_D <- c("cohort", "effect", "CI_low", "CI_high", "P", "interval_N", "person_N")
need_META <- c(
  "estimand", "pooled_effect", "CI_low_mKH", "CI_high_mKH",
  "P_mKH", "I2", "k"
)

check_cols <- function(dt, need, label) {
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0L) {
    stop(
      label, " missing required column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
}

check_cols(A, need_A, "30C")
check_cols(C, need_C, "30D")
check_cols(D, need_D, "30F")
check_cols(META, need_META, "30G")

for (nm in c("A", "C", "D")) {
  z <- get(nm)
  got <- sort(unique(as.character(z$cohort)))
  expected <- sort(EXPECTED_COHORTS)
  if (!identical(got, expected)) {
    stop(
      nm, " does not contain exactly the four expected cohorts. Found: ",
      paste(got, collapse = ", "),
      call. = FALSE
    )
  }
}

sample_qc <- unique(A[, .(cohort, interval_N, person_N)])

EXPECTED_SAMPLE <- data.table(
  cohort = EXPECTED_COHORTS,
  interval_N_expected = c(12020L, 11054L, 30561L, 8355L),
  person_N_expected = c(3843L, 2907L, 11620L, 4920L)
)

sample_qc <- merge(
  sample_qc,
  EXPECTED_SAMPLE,
  by = "cohort",
  all.x = TRUE
)

sample_qc[, exact_match :=
  interval_N == interval_N_expected &
  person_N == person_N_expected
]

if (!all(sample_qc$exact_match)) {
  cat("\nFigure 2 Script-30 sample QC:\n")
  print(sample_qc)
  stop("Locked primary-sample QC failed. Figure 2 generation stopped.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 4. COHORT-SPECIFIC ESTIMATES
# ------------------------------------------------------------------------------
A_plot <- A[, .(
  cohort,
  estimand = "Any limitation vs Independent",
  effect,
  CI_low,
  CI_high,
  P,
  interval_N,
  person_N,
  source = "cohort"
)]

C_plot <- C[
  estimand %in% c("Mild vs Independent", "Severe vs Independent"),
  .(
    cohort,
    estimand = fifelse(
      estimand == "Mild vs Independent",
      "Mild limitation vs Independent",
      "Severe limitation vs Independent"
    ),
    effect,
    CI_low = CI_low_percentile,
    CI_high = CI_high_percentile,
    P = P_boot_normal,
    interval_N,
    person_N,
    source = "cohort"
  )
]

D_plot <- D[, .(
  cohort,
  estimand = "Severe vs Mild",
  effect,
  CI_low,
  CI_high,
  P,
  interval_N,
  person_N,
  source = "cohort"
)]

COHORT_PLOT <- rbindlist(list(A_plot, C_plot, D_plot), fill = TRUE)

# ------------------------------------------------------------------------------
# 5. POOLED META-ANALYSIS ESTIMATES
# ------------------------------------------------------------------------------
meta_use <- merge(
  PANEL_MAP,
  META[, .(
    meta_estimand = estimand,
    pooled_effect,
    CI_low_mKH,
    CI_high_mKH,
    P_mKH,
    I2,
    k
  )],
  by = "meta_estimand",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(meta_use$pooled_effect)) {
  miss <- meta_use[is.na(pooled_effect), meta_estimand]
  stop(
    paste0("Could not match pooled meta-analysis row(s): ", paste(miss, collapse = ", ")),
    call. = FALSE
  )
}

POOLED_PLOT <- meta_use[, .(
  cohort = "Pooled",
  estimand,
  effect = pooled_effect,
  CI_low = CI_low_mKH,
  CI_high = CI_high_mKH,
  P = P_mKH,
  interval_N = NA_integer_,
  person_N = NA_integer_,
  I2,
  k,
  source = "pooled"
)]

PLOT_DT <- rbindlist(list(COHORT_PLOT, POOLED_PLOT), fill = TRUE)
PLOT_DT <- merge(PLOT_DT, PANEL_MAP, by = "estimand", all.x = TRUE, sort = FALSE)

I2_MAP <- meta_use[, .(
  estimand,
  strip_label = paste0(
    panel_letter, "  ", panel_title,
    "\nPooled random-effects I² = ", fmt_num(I2, 1), "%"
  )
)]

PLOT_DT <- merge(PLOT_DT, I2_MAP, by = "estimand", all.x = TRUE, sort = FALSE)

PLOT_DT[, cohort := factor(cohort, levels = rev(COHORT_ORDER))]
PLOT_DT[, panel := factor(strip_label, levels = I2_MAP$strip_label)]
PLOT_DT[, effect_label := fmt_or_ci(effect, CI_low, CI_high, digits = 2)]
PLOT_DT[, p_label := fifelse(is.na(P), "", fifelse(P < 0.001, "P<0.001", paste0("P=", fmt_num(P, 3))))]
PLOT_DT[, label_right := effect_label]

# ------------------------------------------------------------------------------
# 6. AXIS / LABEL RANGE
# ------------------------------------------------------------------------------
xmin_data <- min(PLOT_DT$CI_low, na.rm = TRUE)
xmax_data <- max(PLOT_DT$CI_high, na.rm = TRUE)

x_left <- min(0.90, floor((xmin_data - 0.03) * 20) / 20)
x_core_right <- max(1.45, ceiling((xmax_data + 0.03) * 20) / 20)
x_text <- x_core_right + 0.04
x_right <- x_core_right + 0.34
x_breaks <- seq(floor(x_left * 10) / 10, ceiling(x_core_right * 10) / 10, by = 0.1)

# ------------------------------------------------------------------------------
# 7. EXPORT PLOT DATA / CAPTION
# ------------------------------------------------------------------------------
plot_export <- copy(PLOT_DT)[
  order(panel, cohort),
  .(
    panel_letter,
    panel_title,
    cohort = as.character(cohort),
    estimand,
    effect,
    CI_low,
    CI_high,
    effect_label,
    P,
    p_label,
    interval_N,
    person_N,
    source
  )
]

fwrite(plot_export, file.path(FIG_DIR, "Fig2_Primary_Results_plot_data_UNIFIED_FINAL.csv"))

caption_txt <- paste(
  "Figure 2. Primary unified continuous-burden results across the four cohorts.",
  "Each panel shows odds ratios (95% confidence intervals) per +1 SD increase in the unified cumulative low-PEF burden.",
  "At each PEF assessment, low-PEF deficit was defined as max(0, -residual z); the two deficits were averaged and then standardized within cohort among unique participants in the locked primary sample.",
  "Panel A shows the primary binary end point (any limitation versus remaining independent).",
  "Panels B and C show multinomial severity-model estimates for mild limitation and severe limitation, each versus remaining independent; these cohort-specific confidence intervals are participant-cluster bootstrap intervals.",
  "Panel D shows the formal direct severe-versus-mild contrast based on the severe/mild sub-sample.",
  "Pooled estimates are random-effects meta-analyses using REML with modified Hartung-Knapp inference.",
  "All cohort-specific models use the common adjustment set: age, sex, analysis interval, height, BMI, and education.",
  sep = " "
)
writeLines(caption_txt, con = file.path(FIG_DIR, "Fig2_caption_UNIFIED_FINAL.txt"), useBytes = TRUE)

# ------------------------------------------------------------------------------
# 8. BUILD FIGURE
# ------------------------------------------------------------------------------
base_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle = element_text(size = 10.5, colour = "#5A6B7C"),
    strip.text = element_text(face = "bold", size = 11, hjust = 0),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 10.5, colour = "#1A1A1A"),
    axis.text.x = element_text(size = 9.8, colour = "#1A1A1A"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#D9DEE5", linewidth = 0.45),
    strip.background = element_rect(fill = "#F2F5F9", colour = "#D0D7E2", linewidth = 0.7),
    panel.border = element_rect(fill = NA, colour = "#D0D7E2", linewidth = 0.7),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position = "none",
    plot.margin = margin(10, 110, 10, 12)
  )

p <- ggplot(PLOT_DT, aes(x = effect, y = cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, colour = "#7F7F7F") +
  geom_segment(
    aes(x = CI_low, xend = CI_high, y = cohort, yend = cohort),
    linewidth = 0.8,
    colour = "#2E5F87",
    lineend = "round"
  ) +
  geom_point(
    data = PLOT_DT[source == "cohort"],
    shape = 21,
    size = 3.2,
    stroke = 1.0,
    fill = "white",
    colour = "#2E5F87"
  ) +
  geom_point(
    data = PLOT_DT[source == "pooled"],
    shape = 23,
    size = 4.2,
    stroke = 0.9,
    fill = "#184E77",
    colour = "#184E77"
  ) +
  geom_text(
    aes(x = x_text, label = label_right),
    hjust = 0,
    size = 3.10,
    colour = "#1A1A1A"
  ) +
  facet_wrap(~ panel, ncol = 2, scales = "fixed") +
  scale_x_continuous(
    limits = c(x_left, x_right),
    breaks = x_breaks,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Figure 2. Primary unified continuous-burden results",
    subtitle = paste0(
      "Cohort-specific estimates are shown together with the pooled random-effects estimate in each panel. ",
      "All effects are odds ratios per +1 SD unified cumulative low-PEF burden."
    ),
    x = "Odds ratio (95% CI)"
  ) +
  base_theme

# ------------------------------------------------------------------------------
# 9. EXPORT FIGURE
# ------------------------------------------------------------------------------
OUT_PNG  <- file.path(FIG_DIR, "Fig2_Primary_Results_FINAL_v3.png")
OUT_PDF  <- file.path(FIG_DIR, "Fig2_Primary_Results_FINAL_v3.pdf")
OUT_TIFF <- file.path(FIG_DIR, "Fig2_Primary_Results_FINAL_v3.tiff")

plot_width <- 14.4
plot_height <- 9.8

ggsave(OUT_PNG, plot = p, width = plot_width, height = plot_height, units = "in", dpi = 320, bg = "white")
ggsave(OUT_PDF, plot = p, width = plot_width, height = plot_height, units = "in", device = cairo_pdf, bg = "white")
ggsave(OUT_TIFF, plot = p, width = plot_width, height = plot_height, units = "in", dpi = 600, compression = "lzw", bg = "white")

# ------------------------------------------------------------------------------
# 10. CONSOLE SUMMARY
# ------------------------------------------------------------------------------
cat("\n======================================================================\n")
cat("FIGURE 2 UNIFIED-BURDEN FINAL COMPLETED\n")
cat("======================================================================\n\n")

cat("[1] SCRIPT 30 LOCKED-SAMPLE QC\n\n")
print(sample_qc)

cat("\n[2] SCRIPT 30 STABILITY GATE\n\n")
print(GATE)

cat("\n[3] PANEL-LEVEL META I²\n\n")
print(meta_use[, .(
  panel = paste0(panel_letter, ". ", panel_title),
  pooled_k = k,
  I2 = round(I2, 1),
  pooled_OR_95CI = fmt_or_ci(pooled_effect, CI_low_mKH, CI_high_mKH, digits = 2)
)])

cat("\n[4] OUTPUT FILES\n\n")
cat("  ", OUT_PNG, "\n", sep = "")
cat("  ", OUT_PDF, "\n", sep = "")
cat("  ", OUT_TIFF, "\n", sep = "")
cat("  ", file.path(FIG_DIR, "Fig2_Primary_Results_plot_data_UNIFIED_FINAL.csv"), "\n", sep = "")
cat("  ", file.path(FIG_DIR, "Fig2_caption_UNIFIED_FINAL.txt"), "\n", sep = "")

cat("\nExposure source: Script 30 unified burden.\n")
cat("FINAL STATUS: PASS\n")
