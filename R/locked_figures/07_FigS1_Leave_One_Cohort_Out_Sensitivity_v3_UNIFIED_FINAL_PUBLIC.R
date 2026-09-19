# ==============================================================================
# Supplementary Figure S1 -- unified-burden final v3
# Leave-one-cohort-out robustness of the pooled primary results
# ------------------------------------------------------------------------------
# This script reads the unified-burden meta-analysis outputs from Scripts 30
# and 31 and creates Supplementary Figure S1 showing how the pooled estimates
# change when each cohort is omitted in turn.
#
# Required upstream files:
#   30G_meta_unified_burden.csv
#   31T_leave_one_cohort_out_meta.csv
#   31U_FINAL_UNIFIED_SENSITIVITY_GATE.csv
#
# Suggested manuscript use:
#   Supplementary Figure S1. Leave-one-cohort-out pooled estimates for the
#   primary unified continuous-burden analyses.
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

# ==============================================================================
# 1. PATHS TO EDIT IF NEEDED
# ==============================================================================
PROJECT_ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)

U30_ROOT <- file.path(PROJECT_ROOT, "output", "compat", "unified_burden_check"
)

U31_ROOT <- file.path(PROJECT_ROOT, "output", "compat", "unified_burden_final"
)

U30_TABLE <- file.path(U30_ROOT, "tables")
U31_TABLE <- file.path(U31_ROOT, "tables")
U31_QC <- file.path(U31_ROOT, "QC")

FIG_DIR <- file.path(PROJECT_ROOT, "output", "figures")

PRIMARY_META_FILE <- file.path(
  U30_TABLE,
  "30G_meta_unified_burden.csv"
)

LOO_META_FILE <- file.path(
  U31_TABLE,
  "31T_leave_one_cohort_out_meta.csv"
)

GATE_FILE <- file.path(
  U31_QC,
  "31U_FINAL_UNIFIED_SENSITIVITY_GATE.csv"
)

OUT_PNG <- file.path(FIG_DIR, "FigS1_Leave_One_Cohort_Out_Sensitivity_FINAL.png")
OUT_PDF <- file.path(FIG_DIR, "FigS1_Leave_One_Cohort_Out_Sensitivity_FINAL.pdf")

INPUT_FILES <- c(
  PRIMARY_META_FILE,
  LOO_META_FILE,
  GATE_FILE
)

missing_files <- INPUT_FILES[!file.exists(INPUT_FILES)]

if (length(missing_files) > 0L) {
  stop(
    "Missing Supplementary Figure S1 input file(s):\n",
    paste(missing_files, collapse = "\n"),
    "\n\nRun Scripts 30 and 31 first.",
    call. = FALSE
  )
}

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. READ DATA
# ==============================================================================
primary_meta <- as.data.table(fread(PRIMARY_META_FILE))
loo_meta <- as.data.table(fread(LOO_META_FILE))
GATE <- as.data.table(fread(GATE_FILE))

# Script 30 names the pooled effect "pooled_effect"; Script 31 leave-one-out
# output names the same quantity "pooled_OR".  Standardize only the column name
# for plotting; no estimate is recomputed here.
if (!"pooled_effect" %in% names(primary_meta)) {
  stop("30G is missing pooled_effect.", call. = FALSE)
}

if (!"pooled_OR" %in% names(loo_meta)) {
  stop("31T is missing pooled_OR.", call. = FALSE)
}

loo_meta[, pooled_effect := pooled_OR]

# Script 31 final practical gate must pass before the sensitivity figure is made.
if (!all(c("check", "pass", "value") %in% names(GATE))) {
  stop(
    "31U unified-sensitivity gate has unexpected columns.",
    call. = FALSE
  )
}

gate_pass <- if (is.logical(GATE$pass)) {
  GATE$pass
} else {
  tolower(trimws(as.character(GATE$pass))) %in%
    c("true", "t", "1", "yes")
}

if (!all(gate_pass)) {
  cat("\nScript 31 final sensitivity gate:\n")
  print(GATE)
  stop(
    "Script 31 unified sensitivity gate did not fully pass. FigS1 generation stopped.",
    call. = FALSE
  )
}

# Keep the three key pooled estimands that are most central to the final paper.
keep_estimands <- c(
  "A: Any limitation vs Independent",
  "C: Severe vs Independent",
  "Severity contrast: Severe vs Mild"
)

primary_meta <- primary_meta[estimand %in% keep_estimands]
loo_meta <- loo_meta[estimand %in% keep_estimands]

if (nrow(primary_meta) == 0L || nrow(loo_meta) == 0L) {
  stop("Expected estimands were not found in the primary or LOO meta files.")
}

EXPECTED_COHORTS <- c("HRS", "ELSA", "SHARE", "CHARLS")

if (nrow(primary_meta) != length(keep_estimands)) {
  stop(
    "30G does not contain exactly one row for each of the three FigS1 estimands.",
    call. = FALSE
  )
}

loo_qc <- loo_meta[
  ,
  .(
    row_N = .N,
    omitted_N = uniqueN(omitted_cohort),
    all_positive = all(pooled_effect > 1, na.rm = TRUE)
  ),
  by = estimand
]

if (
  nrow(loo_qc) != length(keep_estimands) ||
  any(loo_qc$row_N != 4L) ||
  any(loo_qc$omitted_N != 4L) ||
  !setequal(unique(loo_meta$omitted_cohort), EXPECTED_COHORTS)
) {
  cat("\nLeave-one-cohort-out structural QC:\n")
  print(loo_qc)
  stop(
    "31T does not contain the expected 3 estimands x 4 omitted cohorts.",
    call. = FALSE
  )
}

# ==============================================================================
# 3. BUILD PLOT DATASET
# ==============================================================================
primary_plot <- copy(primary_meta)
primary_plot[, `:=`(
  analysis = "Full 4-cohort pooled estimate",
  row_type = "Primary pooled"
)]

loo_plot <- copy(loo_meta)
loo_plot[, `:=`(
  analysis = paste0("Omit ", omitted_cohort),
  row_type = "Leave-one-out"
)]

plot_dt <- rbindlist(
  list(
    primary_plot[, .(
      estimand,
      analysis,
      row_type,
      pooled_effect,
      CI_low_mKH,
      CI_high_mKH,
      P_mKH,
      I2
    )],
    loo_plot[, .(
      estimand,
      analysis,
      row_type,
      pooled_effect,
      CI_low_mKH,
      CI_high_mKH,
      P_mKH,
      I2
    )]
  ),
  fill = TRUE
)

panel_map <- data.table(
  estimand = keep_estimands,
  panel = c(
    "A  Any limitation vs Independent",
    "B  Severe limitation vs Independent",
    "C  Severe vs Mild"
  )
)
plot_dt <- merge(plot_dt, panel_map, by = "estimand", all.x = TRUE)

# Order rows within each panel.
row_levels <- c(
  "Full 4-cohort pooled estimate",
  "Omit HRS",
  "Omit ELSA",
  "Omit SHARE",
  "Omit CHARLS"
)
plot_dt[, analysis := factor(analysis, levels = rev(row_levels))]
plot_dt[, panel := factor(panel, levels = c(
  "A  Any limitation vs Independent",
  "B  Severe limitation vs Independent",
  "C  Severe vs Mild"
))]

# Formatted labels shown at the right side of each panel.
plot_dt[, p_text := fifelse(
  is.na(P_mKH),
  "",
  fifelse(P_mKH < 0.001, "P<0.001", paste0("P=", formatC(P_mKH, format = "f", digits = 3)))
)]

plot_dt[, estimate_text := paste0(
  sprintf("%.2f", pooled_effect),
  " (",
  sprintf("%.2f", CI_low_mKH),
  "-",
  sprintf("%.2f", CI_high_mKH),
  ")"
)]

plot_dt[, i2_text := paste0("I",
                            intToUtf8(178),
                            "=",
                            sprintf("%.1f", I2),
                            "%")]
plot_dt[, label := paste0(estimate_text, "; ", i2_text)]

# Dynamic x-range with room for right-side labels.
x_min <- min(c(0.95, plot_dt$CI_low_mKH), na.rm = TRUE)
x_max <- max(plot_dt$CI_high_mKH, na.rm = TRUE)
text_x <- x_max + 0.18
x_upper <- text_x + 0.10

# ==============================================================================
# 4. DRAW FIGURE
# ==============================================================================
p <- ggplot(plot_dt, aes(x = pooled_effect, y = analysis)) +
  facet_wrap(~panel, ncol = 1, scales = "fixed") +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, colour = "grey45") +
  geom_errorbarh(
    aes(xmin = CI_low_mKH, xmax = CI_high_mKH),
    height = 0.16,
    linewidth = 0.9,
    colour = "#2F6DA3"
  ) +
  geom_point(
    data = plot_dt[row_type == "Leave-one-out"],
    shape = 21,
    size = 3.2,
    stroke = 1.0,
    fill = "white",
    colour = "#2F6DA3"
  ) +
  geom_point(
    data = plot_dt[row_type == "Primary pooled"],
    shape = 23,
    size = 4.0,
    stroke = 1.0,
    fill = "#1D4F78",
    colour = "#1D4F78"
  ) +
  geom_text(
    aes(x = text_x, label = label),
    hjust = 0,
    size = 3.4,
    colour = "#333333"
  ) +
  scale_x_continuous(
    limits = c(x_min, x_upper),
    expand = expansion(mult = c(0.01, 0.01)),
    breaks = pretty(c(x_min, x_max), n = 6)
  ) +
  labs(
    title = "Supplementary Figure S1. Leave-one-cohort-out sensitivity analyses",
    subtitle = "Full four-cohort and leave-one-cohort-out pooled ORs per +1 SD unified cumulative low-PEF burden.",
    x = "Odds ratio (95% CI)",
    y = NULL,
    caption = paste(
      "Diamond = full four-cohort pooled estimate; hollow circles = leave-one-cohort-out estimates.",
      "I² is shown for each pooled model."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, colour = "#1F2D3D"),
    plot.subtitle = element_text(size = 11, colour = "#5E6E7E", margin = margin(b = 8)),
    plot.caption = element_text(size = 10, colour = "#5E6E7E", hjust = 0),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10, colour = "#2C3E50"),
    axis.text.y = element_text(size = 11, colour = "#1F2D3D", face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "#D5DBE3", linewidth = 0.45),
    strip.text = element_text(face = "bold", size = 12, colour = "#1F2D3D"),
    strip.background = element_rect(fill = "#EEF2F6", colour = "#C9D2DC", linewidth = 0.8),
    panel.spacing = unit(1.1, "lines"),
    plot.margin = margin(14, 150, 12, 12)
  ) +
  coord_cartesian(clip = "off")

# ==============================================================================
# 5. SAVE
# ==============================================================================
ggsave(
  filename = OUT_PNG,
  plot = p,
  width = 12.2,
  height = 11.2,
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = OUT_PDF,
  plot = p,
  width = 12.2,
  height = 11.2,
  device = cairo_pdf,
  bg = "white"
)

cat("============================================================\n")
cat("Supplementary Figure S1 unified-burden final completed.\n")
cat("============================================================\n")

cat("\n[1] SCRIPT 31 FINAL SENSITIVITY GATE\n\n")
print(GATE)

cat("\n[2] LEAVE-ONE-COHORT-OUT STRUCTURAL QC\n\n")
print(loo_qc)

cat("\n[3] FULL AND LOO PLOT DATA\n\n")
print(
  plot_dt[
    ,
    .(
      estimand,
      analysis = as.character(analysis),
      pooled_effect,
      CI_low_mKH,
      CI_high_mKH,
      I2
    )
  ]
)

cat("\n[4] OUTPUT FILES\n\n")
cat("PNG: ", OUT_PNG, "\n", sep = "")
cat("PDF: ", OUT_PDF, "\n", sep = "")
cat("\nExposure source: Script 30 unified primary estimates + Script 31 unified LOO reanalysis.\n")
cat("FINAL STATUS: PASS\n")
cat("============================================================\n")
