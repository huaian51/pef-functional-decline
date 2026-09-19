# ============================================================================== 
# Supplementary Figure S2
# Pooled unified-burden sensitivity analyses for the final four-cohort paper
# ------------------------------------------------------------------------------
# This script summarizes the main pooled robustness analyses around the primary
# unified continuous cumulative low-PEF burden exposure.
#
# Figure structure:
#   Panel A: Any limitation vs Independent
#   Panel B: Severe limitation vs Independent
#   Panel C: Severe vs Mild
#
# Rows shown within panels:
#   - Primary common-core model
#   - Primary on smoking-observed sample
#   - + smoking adjustment
#   - Cohort-specific extended model
#   - IPW for observed follow-up (available for Panel A and Panel C)
#
# Required upstream files:
#   30G_meta_unified_burden.csv
#   31F_smoking_extended_adjustment_meta.csv
#   31O_IPW_observed_followup_meta.csv
#   31U_FINAL_UNIFIED_SENSITIVITY_GATE.csv
#
# Suggested manuscript use:
#   Supplementary Figure S2. Robustness of the pooled primary associations to
#   smoking adjustment, cohort-specific extended covariate sets, and follow-up
#   observation weighting.
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

DEFENSIVE_META_FILE <- file.path(
  U31_TABLE,
  "31F_smoking_extended_adjustment_meta.csv"
)

IPW_META_FILE <- file.path(
  U31_TABLE,
  "31O_IPW_observed_followup_meta.csv"
)

GATE_FILE <- file.path(
  U31_QC,
  "31U_FINAL_UNIFIED_SENSITIVITY_GATE.csv"
)

OUT_PNG <- file.path(FIG_DIR, "FigS2_Pooled_Sensitivity_Analyses_v2_FINAL.png")
OUT_PDF <- file.path(FIG_DIR, "FigS2_Pooled_Sensitivity_Analyses_v2_FINAL.pdf")
OUT_CSV <- file.path(FIG_DIR, "FigS2_Pooled_Sensitivity_Analyses_plot_data_UNIFIED_FINAL.csv")

need_files <- c(
  PRIMARY_META_FILE,
  DEFENSIVE_META_FILE,
  IPW_META_FILE,
  GATE_FILE
)
missing_files <- need_files[!file.exists(need_files)]
if (length(missing_files)) {
  stop(
    "Required file(s) missing:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nRun Scripts 30 and 31 first.",
    call. = FALSE
  )
}
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================== 
# 2. HELPERS
# ============================================================================== 
`%||%` <- function(x, y) {
  if (!is.null(x) && length(x)) x else y
}

pick_col <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)]
  if (!length(hit)) return(NULL)
  hit[1]
}

std_or_meta <- function(dt) {
  dt <- as.data.table(copy(dt))

  or_col   <- pick_col(dt, c("pooled_OR", "pooled_effect", "OR", "effect"))
  low_col  <- pick_col(dt, c("CI_low_mKH", "CI_low", "lower", "lower_CI"))
  high_col <- pick_col(dt, c("CI_high_mKH", "CI_high", "upper", "upper_CI"))
  p_col    <- pick_col(dt, c("P_mKH", "P", "p", "p_value"))
  i2_col   <- pick_col(dt, c("I2", "I2_pct", "I2_percent"))

  if (is.null(or_col) || is.null(low_col) || is.null(high_col)) {
    stop("Meta file is missing pooled OR / CI columns.")
  }

  dt[, pooled_or := as.numeric(get(or_col))]
  dt[, ci_low    := as.numeric(get(low_col))]
  dt[, ci_high   := as.numeric(get(high_col))]
  dt[, p_value   := if (!is.null(p_col)) as.numeric(get(p_col)) else NA_real_]
  dt[, i2_value  := if (!is.null(i2_col)) as.numeric(get(i2_col)) else NA_real_]

  dt
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "",
    ifelse(p < 0.001, "P<0.001", paste0("P=", formatC(p, format = "f", digits = 3)))
  )
}

fmt_i2 <- function(i2) {
  ifelse(
    is.na(i2),
    "",
    paste0("I", intToUtf8(178), "=", formatC(i2, format = "f", digits = 1), "%")
  )
}

# ============================================================================== 
# 3. READ + STANDARDIZE META FILES
# ============================================================================== 
primary_meta <- std_or_meta(fread(PRIMARY_META_FILE))
defensive_meta <- std_or_meta(fread(DEFENSIVE_META_FILE))
ipw_meta <- std_or_meta(fread(IPW_META_FILE))
GATE <- as.data.table(fread(GATE_FILE))

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
    "Script 31 unified sensitivity gate did not fully pass. FigS2 generation stopped.",
    call. = FALSE
  )
}

# ============================================================================== 
# 4. BUILD FIGURE DATASET
# ============================================================================== 
# ---- 4A. Primary pooled results from Script 30
primary_map <- data.table(
  estimand = c(
    "A: Any limitation vs Independent",
    "C: Severe vs Independent",
    "Severity contrast: Severe vs Mild"
  ),
  panel = c(
    "A  Any limitation vs Independent",
    "B  Severe limitation vs Independent",
    "C  Severe vs Mild"
  )
)

primary_plot <- merge(primary_meta, primary_map, by = "estimand", all.y = FALSE)
if (nrow(primary_plot) == 0L) {
  stop("Could not find the expected primary pooled estimands in 30G_meta_unified_burden.csv")
}
primary_plot <- primary_plot[, .(
  panel,
  analysis = "Primary common-core model",
  pooled_or,
  ci_low,
  ci_high,
  p_value,
  i2_value,
  source_type = "Primary"
)]

# ---- 4B. Smoking / extended-adjustment meta results from Script 31
spec_map <- data.table(
  spec = c(
    "CORE0_ON_SMOKING_OBSERVED_SAMPLE",
    "CORE0_PLUS_CURRENT_SMOKING",
    "COHORT_SPECIFIC_EXTENDED_RELIABLE"
  ),
  analysis = c(
    "Primary on smoking-observed sample",
    "+ smoking adjustment",
    "Cohort-specific extended model"
  )
)

estimand_map_def <- data.table(
  estimand = c(
    "A_ANY_LIMITATION",
    "SEVERE_VS_INDEPENDENT",
    "SEVERE_VS_MILD"
  ),
  panel = c(
    "A  Any limitation vs Independent",
    "B  Severe limitation vs Independent",
    "C  Severe vs Mild"
  )
)

defensive_plot <- merge(defensive_meta, spec_map, by = "spec", all = FALSE)
defensive_plot <- merge(defensive_plot, estimand_map_def, by = "estimand", all = FALSE)
if (nrow(defensive_plot) == 0L) {
  stop("Could not find the expected defensive pooled analyses in 31F_smoking_extended_adjustment_meta.csv")
}
defensive_plot <- defensive_plot[, .(
  panel,
  analysis,
  pooled_or,
  ci_low,
  ci_high,
  p_value,
  i2_value,
  source_type = "Sensitivity"
)]

# ---- 4C. IPW pooled results from Script 31
ipw_map <- data.table(
  analysis = c(
    "IPW_FOR_OBSERVATION_ANY_LIMITATION",
    "IPW_FOR_OBSERVATION_SEVERE_VS_MILD"
  ),
  panel = c(
    "A  Any limitation vs Independent",
    "C  Severe vs Mild"
  ),
  analysis_label = c(
    "IPW for observed follow-up",
    "IPW for observed follow-up"
  )
)

ipw_plot <- merge(ipw_meta, ipw_map, by = "analysis", all = FALSE)
if (nrow(ipw_plot) == 0L) {
  warning("No eligible IPW pooled rows found; the figure will be drawn without IPW rows.")
  ipw_plot <- data.table()
} else {
  ipw_plot <- ipw_plot[, .(
    panel,
    analysis = analysis_label,
    pooled_or,
    ci_low,
    ci_high,
    p_value,
    i2_value,
    source_type = "Sensitivity"
  )]
}

# ---- 4D. Structural QC before combining
# Three primary estimands are shown.
if (nrow(primary_plot) != 3L) {
  stop(
    "FigS2 expected exactly 3 primary pooled rows from 30G.",
    call. = FALSE
  )
}

# Three sensitivity specifications x three estimands = nine rows.
expected_defensive_rows <- 9L
if (nrow(defensive_plot) != expected_defensive_rows) {
  cat("\nDefensive sensitivity rows found:\n")
  print(defensive_plot)
  stop(
    "FigS2 expected exactly 9 smoking/extended rows from 31F.",
    call. = FALSE
  )
}

# IPW is deliberately shown only for estimands directly matching the panels:
# Any limitation vs Independent and Severe vs Mild.
# The Script 31 severe-vs-nonsevere IPW analysis is a different estimand from
# Panel B (Severe vs Independent) and therefore must not be inserted here.
if (nrow(ipw_plot) != 2L) {
  cat("\nEligible IPW rows found:\n")
  print(ipw_plot)
  stop(
    "FigS2 expected exactly 2 directly matched IPW rows from 31O.",
    call. = FALSE
  )
}

# ---- 4E. Combine
plot_dt <- rbindlist(
  list(primary_plot, defensive_plot, ipw_plot),
  fill = TRUE,
  use.names = TRUE
)

if (nrow(plot_dt) == 0L) {
  stop("No rows available for plotting.")
}

# Remove duplicates if the same panel-analysis combination somehow appears twice.
setorder(plot_dt, panel, analysis)
plot_dt <- unique(plot_dt, by = c("panel", "analysis"))

row_levels <- c(
  "Primary common-core model",
  "Primary on smoking-observed sample",
  "+ smoking adjustment",
  "Cohort-specific extended model",
  "IPW for observed follow-up"
)

plot_dt[, analysis := factor(analysis, levels = rev(row_levels))]
plot_dt[, panel := factor(panel, levels = c(
  "A  Any limitation vs Independent",
  "B  Severe limitation vs Independent",
  "C  Severe vs Mild"
))]
plot_dt[, source_type := factor(source_type, levels = c("Sensitivity", "Primary"))]

plot_dt[, estimate_text := paste0(
  sprintf("%.2f", pooled_or),
  " (",
  sprintf("%.2f", ci_low),
  "-",
  sprintf("%.2f", ci_high),
  ")"
)]

plot_dt[, stat_text := ifelse(
  nzchar(fmt_i2(i2_value)),
  paste0(estimate_text, "; ", fmt_i2(i2_value)),
  estimate_text
)]

# ------------------------------------------------------------------------------
# Axis range and right-margin label placement.
# ------------------------------------------------------------------------------
x_min <- min(c(0.95, plot_dt$ci_low), na.rm = TRUE)
x_max <- max(plot_dt$ci_high, na.rm = TRUE)
text_x <- x_max + 0.20
x_upper <- text_x + 0.17

# Save plot dataset for transparency.
fwrite(plot_dt, OUT_CSV)

# ============================================================================== 
# 5. DRAW FIGURE
# ============================================================================== 
p <- ggplot(plot_dt, aes(x = pooled_or, y = analysis)) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey45"
  ) +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    height = 0.16,
    linewidth = 0.95,
    colour = "#2F6DA3"
  ) +
  geom_point(
    data = plot_dt[source_type == "Sensitivity"],
    shape = 21,
    size = 3.2,
    stroke = 1.0,
    fill = "white",
    colour = "#2F6DA3"
  ) +
  geom_point(
    data = plot_dt[source_type == "Primary"],
    shape = 23,
    size = 4.0,
    stroke = 1.0,
    fill = "#1D4F78",
    colour = "#1D4F78"
  ) +
  geom_text(
    aes(x = text_x, label = stat_text),
    hjust = 0,
    size = 3.35,
    colour = "#333333"
  ) +
  scale_x_continuous(
    limits = c(x_min, x_upper),
    expand = expansion(mult = c(0.01, 0.01)),
    breaks = pretty(c(x_min, x_max), n = 6)
  ) +
  labs(
    title = "Supplementary Figure S2. Pooled sensitivity analyses",
    subtitle = "Unified primary pooled estimates compared with smoking, extended-adjustment, and follow-up-weighting sensitivity models.",
    x = "Odds ratio (95% CI)",
    y = NULL,
    caption = paste(
      "Diamond = unified primary pooled model; hollow circles = sensitivity estimates.",
      "IPW is shown only for directly matched estimands; labels show OR (95% CI) and I²."
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
    plot.margin = margin(14, 180, 12, 12)
  ) +
  coord_cartesian(clip = "off")

# ============================================================================== 
# 6. SAVE
# ============================================================================== 
ggsave(
  filename = OUT_PNG,
  plot = p,
  width = 13.4,
  height = 11.8,
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = OUT_PDF,
  plot = p,
  width = 13.4,
  height = 11.8,
  bg = "white"
)

cat("============================================================\n")
cat("Supplementary Figure S2 unified-burden final completed.\n")
cat("============================================================\n")

cat("\n[1] SCRIPT 31 FINAL SENSITIVITY GATE\n\n")
print(GATE)

cat("\n[2] FIGS2 PLOT DATA\n\n")
print(
  plot_dt[
    ,
    .(
      panel = as.character(panel),
      analysis = as.character(analysis),
      pooled_or,
      ci_low,
      ci_high,
      p_value,
      i2_value,
      source_type = as.character(source_type)
    )
  ]
)

cat("\n[3] OUTPUT FILES\n\n")
cat("PNG : ", OUT_PNG, "\n", sep = "")
cat("PDF : ", OUT_PDF, "\n", sep = "")
cat("DATA: ", OUT_CSV, "\n", sep = "")

cat("\nSources:\n")
cat("  Primary pooled : Script 30 / 30G\n")
cat("  Smoking/extended: Script 31 / 31F\n")
cat("  IPW             : Script 31 / 31O\n")
cat("FINAL STATUS: PASS\n")
cat("============================================================\n")
