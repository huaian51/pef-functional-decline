# ==============================================================================
# FIGURE 1 UNIFIED FINAL -- STUDY TIMELINE, ANALYTICAL FRAMEWORK, AND FOUR-COHORT SYNTHESIS
#
# Manuscript:
# Cumulative low peak expiratory flow burden and severity of subsequent
# functional decline: a four-cohort longitudinal study
#
# OUTPUT DIRECTORY:
# output/figures
#
# OUTPUTS:
#   Fig1_Study_Design_FINAL.pdf
#   Fig1_Study_Design_FINAL.tiff
#   Fig1_Study_Design_FINAL.png
#   Fig1_primary_sample_summary_FINAL.csv
#   Fig1_caption_FINAL.txt
#
# PANEL A:
#   Study timeline adapted from the earlier project timeline.
#   Circle  = first PEF assessment
#   Diamond = second PEF assessment / index
#   Square  = subsequent functional-state assessment
#
# PANEL B:
#   Exposure construction + functional-transition framework
#
# PANEL C:
#   Four-cohort primary samples -> harmonized cohort-specific estimation
#   -> REML random-effects meta-analysis with modified Hartung-Knapp inference
#
# IMPORTANT:
# - The primary exposure is the Script 30 unified continuous low-PEF burden.
# - At each PEF assessment: deficit = max(0, -residual_z), where residual_z
#   comes from the pooled two-wave cohort-specific model PEF ~ age + sex + height.
# - Unified raw burden = mean(deficit_1, deficit_2), then standardized among
#   unique persons in the locked cohort-specific primary sample.
# - The 0/1/2 exposure is supportive only and must NOT be called a trajectory.
# - All primary transition intervals begin in functional independence.
# ==============================================================================


# ==============================================================================
# 0. CLEAN SESSION / PACKAGES
# ==============================================================================

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  warn = 1,
  width = 260,
  scipen = 999
)

pkgs <- c(
  "ggplot2",
  "data.table"
)

missing_pkgs <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing package(s): ",
    paste(
      missing_pkgs,
      collapse = ", "
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(grid)
})


# ==============================================================================
# 1. PATHS / SETTINGS
# ==============================================================================

PROJECT_ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)

FIG_DIR <- file.path(
  PROJECT_ROOT,
  "output",
  "figures"
)

DATA_DIR <- file.path(PROJECT_ROOT, "output", "compat", "unified_burden_check",
  "data"
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

COHORTS <- c(
  "HRS",
  "ELSA",
  "SHARE",
  "CHARLS"
)

FILES <- setNames(
  file.path(
    DATA_DIR,
    paste0(
      "30_",
      COHORTS,
      "_primary_analysis_UNIFIED_burden.rds"
    )
  ),
  COHORTS
)

CORE0 <- c(
  "age",
  "sex",
  "interval_f",
  "height",
  "bmi",
  "education"
)

# If the manuscript explicitly discusses excluded 2020 waves, set TRUE.
# Otherwise FALSE gives the cleanest main-text figure.
SHOW_EXCLUDED <- FALSE


# ==============================================================================
# 2. FROZEN PRIMARY-SAMPLE TARGETS
# ==============================================================================

EXPECTED <- data.table(
  cohort = COHORTS,
  interval_N = c(
    12020L,
    11054L,
    30561L,
    8355L
  ),
  person_N = c(
    3843L,
    2907L,
    11620L,
    4920L
  ),
  independent_N = c(
    10751L,
    9994L,
    28290L,
    6528L
  ),
  mild_N = c(
    799L,
    641L,
    1482L,
    971L
  ),
  severe_N = c(
    470L,
    419L,
    789L,
    856L
  )
)


# ==============================================================================
# 3. READ LOCKED PRIMARY ANALYSIS SAMPLES
# ==============================================================================

missing_files <- FILES[
  !file.exists(
    FILES
  )
]

if (length(missing_files) > 0L) {
  stop(
    paste0(
      "Missing reviewer-facing dataset(s):\n",
      paste(
        missing_files,
        collapse = "\n"
      ),
      "\nRun Script 30 (30_UNIFIED_CONTINUOUS_PEF_BURDEN_REANALYSIS.R) first."
    ),
    call. = FALSE
  )
}

sample_rows <- list()

for (coh in COHORTS) {

  d <- as.data.table(
    readRDS(
      FILES[[coh]]
    )
  )

  need <- c(
    "id",
    "burden_z",
    "burden_raw_unified",
    "burden_z_unified",
    "destination",
    CORE0
  )

  miss <- setdiff(
    need,
    names(d)
  )

  if (length(miss) > 0L) {
    stop(
      paste0(
        "[",
        coh,
        "] missing variable(s): ",
        paste(
          miss,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  # Script 30 deliberately writes the unified standardized exposure back to
  # burden_z.  Confirm that the plotting input really follows that contract.
  same_unified <- isTRUE(
    all.equal(
      as.numeric(d$burden_z),
      as.numeric(d$burden_z_unified),
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  )

  if (!same_unified) {
    stop(
      "[", coh,
      "] burden_z is not identical to burden_z_unified in the Script 30 dataset.",
      call. = FALSE
    )
  }

  d[
    ,
    destination := as.character(
      destination
    )
  ]

  d_primary <- d[
    destination %in%
      c(
        "Independent",
        "Mild limitation",
        "Severe limitation"
      ) &
    complete.cases(
      d[
        ,
        c(
          "id",
          "burden_z",
          "destination",
          CORE0
        ),
        with = FALSE
      ]
    )
  ]

  sample_rows[[
    length(sample_rows) + 1L
  ]] <- data.table(
    cohort = coh,
    interval_N = nrow(d_primary),
    person_N = uniqueN(
      d_primary$id
    ),
    independent_N = sum(
      d_primary$destination ==
        "Independent"
    ),
    mild_N = sum(
      d_primary$destination ==
        "Mild limitation"
    ),
    severe_N = sum(
      d_primary$destination ==
        "Severe limitation"
    )
  )
}

sample_summary <- rbindlist(
  sample_rows
)

sample_summary[
  ,
  cohort := factor(
    cohort,
    levels = COHORTS
  )
]

setorder(
  sample_summary,
  cohort
)

sample_summary[
  ,
  cohort := as.character(
    cohort
  )
]


# ==============================================================================
# 4. QC -- STOP IF PRIMARY SAMPLE DOES NOT MATCH FROZEN ANALYSIS
# ==============================================================================

qc <- merge(
  sample_summary,
  EXPECTED,
  by = "cohort",
  suffixes = c(
    "_observed",
    "_expected"
  ),
  all.x = TRUE
)

qc[
  ,
  exact_match :=
    interval_N_observed ==
      interval_N_expected &
    person_N_observed ==
      person_N_expected &
    independent_N_observed ==
      independent_N_expected &
    mild_N_observed ==
      mild_N_expected &
    severe_N_observed ==
      severe_N_expected
]

cat("\n============================================================\n")
cat("FIGURE 1 FINAL PRIMARY-SAMPLE QC\n")
cat("============================================================\n")
print(qc)

if (!all(qc$exact_match)) {
  stop(
    "Primary-sample QC failed. Figure generation stopped.",
    call. = FALSE
  )
}

fwrite(
  sample_summary,
  file.path(
    FIG_DIR,
    "Fig1_primary_sample_summary_FINAL.csv"
  ),
  bom = TRUE
)


# ==============================================================================
# 5. TIMELINE DATA
# ==============================================================================
#
# This preserves the earlier manuscript's study-window chronology while
# simplifying the visual language:
#
#   first PEF assessment -> second PEF assessment / index -> outcome waves
#
# If the final Methods use a different calendar-wave mapping, edit ONLY this
# block; the rest of the plotting code does not need to change.
# ==============================================================================

timeline_rows <- data.table(
  window = c(
    "HRS-A",
    "HRS-B",
    "ELSA",
    "SHARE-A",
    "SHARE-B",
    "CHARLS"
  ),

  pef1_year = c(
    2006,
    2010,
    2004.5,
    2006.5,
    2011,
    2011
  ),

  index_year = c(
    2010,
    2014,
    2008.5,
    2011,
    2015,
    2013
  )
)

outcome_rows <- rbindlist(
  list(
    data.table(
      window = "HRS-A",
      year = c(
        2012,
        2014
      )
    ),
    data.table(
      window = "HRS-B",
      year = c(
        2016,
        2018
      )
    ),
    data.table(
      window = "ELSA",
      year = c(
        2010.5,
        2012.5,
        2014.5,
        2016.5,
        2018.5
      )
    ),
    data.table(
      window = "SHARE-A",
      year = c(
        2013,
        2015,
        2017.5,
        2019.5
      )
    ),
    data.table(
      window = "SHARE-B",
      year = c(
        2017.5,
        2019.5
      )
    ),
    data.table(
      window = "CHARLS",
      year = c(
        2015,
        2018
      )
    )
  )
)

excluded_rows <- data.table(
  window = c(
    "HRS-B",
    "CHARLS"
  ),
  year = c(
    2020,
    2020
  )
)

window_levels <- rev(
  timeline_rows$window
)

timeline_rows[
  ,
  y := match(
    window,
    window_levels
  )
]

outcome_rows[
  ,
  y := match(
    window,
    window_levels
  )
]

excluded_rows[
  ,
  y := match(
    window,
    window_levels
  )
]

# End of follow-up used only to draw the pale horizontal backbone.
followup_end <- outcome_rows[
  ,
  .(
    followup_end = max(
      year
    )
  ),
  by = window
]

if (SHOW_EXCLUDED) {
  followup_end <- merge(
    followup_end,
    excluded_rows[
      ,
      .(
        window,
        excluded_year = year
      )
    ],
    by = "window",
    all.x = TRUE
  )

  followup_end[
    !is.na(
      excluded_year
    ),
    followup_end := pmax(
      followup_end,
      excluded_year
    )
  ]

  followup_end[
    ,
    excluded_year := NULL
  ]
}

timeline_rows <- merge(
  timeline_rows,
  followup_end,
  by = "window",
  all.x = TRUE,
  sort = FALSE
)


# ==============================================================================
# 6. VISUAL SETTINGS
# ==============================================================================

BASE_FAMILY <- "Arial"

COL_TEXT    <- "#17232F"
COL_LIGHT   <- "#6D7885"
COL_GRID    <- "#E7EBEF"
COL_LINE    <- "#C9CED4"
COL_PEF     <- "#2C6DB2"
COL_OUTCOME <- "#D76A15"
COL_PRIMARY <- "#2C5E7D"
COL_SUPPORT <- "#8A6A3C"
COL_BOX     <- "#EEF3F6"
COL_BOX2    <- "#F6F1E9"
COL_SEVERE  <- "#DCE6EC"
COL_BORDER  <- "#31404C"
COL_META    <- "#E7F0F5"

theme_clean <- theme_void(
  base_family = BASE_FAMILY
) +
  theme(
    plot.title = element_text(
      family = BASE_FAMILY,
      face = "bold",
      size = 13.5,
      colour = COL_TEXT,
      hjust = 0,
      margin = margin(
        b = 8
      )
    ),
    plot.margin = margin(
      7,
      9,
      7,
      9
    )
  )


# ==============================================================================
# 7. PANEL A -- STUDY TIMELINE
# ==============================================================================

pA <- ggplot()

# Pale study backbone.
pA <- pA +
  geom_segment(
    data = timeline_rows,
    aes(
      x = pef1_year,
      xend = followup_end,
      y = y,
      yend = y
    ),
    linewidth = 1.0,
    colour = COL_LINE
  )

# Exposure-history segment.
pA <- pA +
  geom_segment(
    data = timeline_rows,
    aes(
      x = pef1_year,
      xend = index_year,
      y = y,
      yend = y
    ),
    linewidth = 1.4,
    colour = COL_PEF,
    alpha = 0.72
  )

# Post-index follow-up segment.
post_segments <- merge(
  timeline_rows[
    ,
    .(
      window,
      y,
      index_year
    )
  ],
  outcome_rows[
    ,
    .(
      first_outcome = min(
        year
      ),
      last_outcome = max(
        year
      )
    ),
    by = window
  ],
  by = "window",
  all.x = TRUE
)

pA <- pA +
  geom_segment(
    data = post_segments,
    aes(
      x = first_outcome,
      xend = last_outcome,
      y = y,
      yend = y
    ),
    linewidth = 1.35,
    colour = COL_OUTCOME,
    alpha = 0.72
  )

# First PEF assessment.
pA <- pA +
  geom_point(
    data = timeline_rows,
    aes(
      x = pef1_year,
      y = y
    ),
    shape = 21,
    size = 4.1,
    stroke = 1.25,
    fill = "white",
    colour = COL_PEF
  )

# Second PEF assessment / index.
pA <- pA +
  geom_segment(
    data = timeline_rows,
    aes(
      x = index_year,
      xend = index_year,
      y = y - 0.28,
      yend = y + 0.28
    ),
    linewidth = 1.1,
    colour = "#111111"
  ) +
  geom_point(
    data = timeline_rows,
    aes(
      x = index_year,
      y = y
    ),
    shape = 23,
    size = 4.4,
    stroke = 1.15,
    fill = "white",
    colour = "#111111"
  )

# Functional-state assessment waves.
pA <- pA +
  geom_point(
    data = outcome_rows,
    aes(
      x = year,
      y = y
    ),
    shape = 22,
    size = 3.5,
    stroke = 1.15,
    fill = "white",
    colour = COL_OUTCOME
  )

# Optional excluded-wave markers.
if (SHOW_EXCLUDED) {

  pA <- pA +
    geom_point(
      data = excluded_rows,
      aes(
        x = year,
        y = y
      ),
      shape = 4,
      size = 4.7,
      stroke = 1.15,
      colour = "#777777"
    )
}

# Legend is built from invisible points for full control.
legend_rows <- data.table(
  x = c(
    2004.1,
    2007.3,
    2012.4
  ),
  y = rep(
    0.25,
    3
  ),
  label = c(
    "First PEF assessment",
    "Second PEF assessment / index",
    "Functional-state assessment"
  ),
  type = c(
    "pef",
    "index",
    "outcome"
  )
)

pA <- pA +
  geom_point(
    data = legend_rows[
      type ==
        "pef"
    ],
    aes(
      x = x,
      y = y
    ),
    shape = 21,
    size = 3.8,
    stroke = 1.15,
    fill = "white",
    colour = COL_PEF
  ) +
  geom_point(
    data = legend_rows[
      type ==
        "index"
    ],
    aes(
      x = x,
      y = y
    ),
    shape = 23,
    size = 4.0,
    stroke = 1.1,
    fill = "white",
    colour = "#111111"
  ) +
  geom_point(
    data = legend_rows[
      type ==
        "outcome"
    ],
    aes(
      x = x,
      y = y
    ),
    shape = 22,
    size = 3.4,
    stroke = 1.1,
    fill = "white",
    colour = COL_OUTCOME
  ) +
  geom_text(
    data = legend_rows,
    aes(
      x = x + 0.45,
      y = y,
      label = label
    ),
    family = BASE_FAMILY,
    hjust = 0,
    vjust = 0.45,
    size = 3.2,
    colour = COL_TEXT
  )

if (SHOW_EXCLUDED) {

  pA <- pA +
    annotate(
      "point",
      x = 2017.7,
      y = 0.25,
      shape = 4,
      size = 4.3,
      stroke = 1.1,
      colour = "#777777"
    ) +
    annotate(
      "text",
      x = 2018.15,
      y = 0.25,
      label = "Excluded wave",
      family = BASE_FAMILY,
      hjust = 0,
      size = 3.2,
      colour = COL_TEXT
    )
}

pA <- pA +
  scale_y_continuous(
    breaks = timeline_rows$y,
    labels = timeline_rows$window,
    limits = c(
      0,
      max(
        timeline_rows$y
      ) + 0.45
    )
  ) +
  scale_x_continuous(
    breaks = seq(
      2004,
      2020,
      by = 4
    ),
    minor_breaks = seq(
      2004,
      2020,
      by = 2
    ),
    limits = c(
      2002.8,
      2020.7
    ),
    expand = c(
      0,
      0
    )
  ) +
  labs(
    title = "A  Study timeline",
    x = "Calendar year",
    y = NULL
  ) +
  theme_minimal(
    base_family = BASE_FAMILY,
    base_size = 11.5
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13.5,
      hjust = 0,
      colour = COL_TEXT,
      margin = margin(
        b = 5
      )
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 10.7,
      colour = COL_TEXT,
      margin = margin(
        r = 8
      )
    ),
    axis.text.x = element_text(
      size = 9.6,
      colour = COL_TEXT
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 10.7,
      colour = COL_TEXT,
      margin = margin(
        t = 7
      )
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(
      colour = COL_GRID,
      linewidth = 0.5
    ),
    panel.grid.minor.x = element_line(
      colour = "#F1F3F5",
      linewidth = 0.35
    ),
    axis.line.x = element_line(
      colour = "#111111",
      linewidth = 0.7
    ),
    axis.ticks.x = element_line(
      colour = "#111111",
      linewidth = 0.6
    ),
    plot.margin = margin(
      6,
      8,
      2,
      8
    )
  )


# ==============================================================================
# 8. PANEL B -- EXPOSURE + FUNCTIONAL-TRANSITION FRAMEWORK
# ==============================================================================

blank_panel <- function(
  xmax,
  ymax,
  title
) {

  ggplot() +
    coord_cartesian(
      xlim = c(
        0,
        xmax
      ),
      ylim = c(
        0,
        ymax
      ),
      clip = "off",
      expand = FALSE
    ) +
    theme_clean +
    labs(
      title = title
    )
}


add_box <- function(
  p,
  xmin,
  xmax,
  ymin,
  ymax,
  label,
  fill = "white",
  colour = COL_BORDER,
  linewidth = 0.55,
  text_size = 3.25,
  text_colour = COL_TEXT,
  fontface = "plain",
  lineheight = 1.05
) {

  p +
    annotate(
      "rect",
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill,
      colour = colour,
      linewidth = linewidth
    ) +
    annotate(
      "text",
      x = (
        xmin +
          xmax
      ) / 2,
      y = (
        ymin +
          ymax
      ) / 2,
      label = label,
      family = BASE_FAMILY,
      size = text_size,
      colour = text_colour,
      fontface = fontface,
      lineheight = lineheight
    )
}


add_arrow <- function(
  p,
  x,
  xend,
  y,
  yend,
  colour = COL_BORDER,
  linewidth = 0.62
) {

  p +
    annotate(
      "segment",
      x = x,
      xend = xend,
      y = y,
      yend = yend,
      colour = colour,
      linewidth = linewidth,
      arrow = grid::arrow(
        length = unit(
          0.15,
          "cm"
        ),
        type = "closed"
      )
    )
}


pB <- blank_panel(
  xmax = 12,
  ymax = 10,
  title = "B  Exposure and functional-transition framework"
)

# Left: repeated PEF -> burden.
pB <- add_box(
  pB,
  xmin = 0.35,
  xmax = 2.65,
  ymin = 7.25,
  ymax = 8.55,
  label = "PEF 1",
  fill = "#F4F6F8",
  fontface = "bold"
)

pB <- add_box(
  pB,
  xmin = 0.35,
  xmax = 2.65,
  ymin = 5.25,
  ymax = 6.55,
  label = "PEF 2",
  fill = "#F4F6F8",
  fontface = "bold"
)

# Merge node.
pB <- pB +
  annotate(
    "point",
    x = 3.45,
    y = 6.9,
    shape = 21,
    size = 3.4,
    stroke = 0.9,
    fill = "white",
    colour = COL_PRIMARY
  )

pB <- add_arrow(
  pB,
  x = 2.65,
  xend = 3.28,
  y = 7.9,
  yend = 7.03
)

pB <- add_arrow(
  pB,
  x = 2.65,
  xend = 3.28,
  y = 5.9,
  yend = 6.77
)

pB <- add_box(
  pB,
  xmin = 4.05,
  xmax = 7.45,
  ymin = 6.55,
  ymax = 8.55,
  label = paste0(
    "PRIMARY EXPOSURE\n",
    "Mean of 2 low-PEF deficits\n",
    "max(0, -residual z)\n",
    "standardized within cohort"
  ),
  fill = "#E7F0F5",
  colour = COL_PRIMARY,
  linewidth = 0.75,
  text_colour = COL_PRIMARY,
  fontface = "bold",
  text_size = 2.85
)

pB <- add_box(
  pB,
  xmin = 4.05,
  xmax = 7.45,
  ymin = 3.85,
  ymax = 5.55,
  label = paste0(
    "SUPPORTIVE EXPOSURE\n",
    "Low-PEF observations: 0 / 1 / 2"
  ),
  fill = COL_BOX2,
  colour = COL_SUPPORT,
  linewidth = 0.65,
  text_colour = COL_SUPPORT,
  fontface = "bold",
  text_size = 3.05
)

pB <- add_arrow(
  pB,
  x = 3.62,
  xend = 4.0,
  y = 6.95,
  yend = 7.55
)

pB <- add_arrow(
  pB,
  x = 3.62,
  xend = 4.0,
  y = 6.85,
  yend = 4.7
)

# Right: independent origin -> destination states.
pB <- add_box(
  pB,
  xmin = 8.05,
  xmax = 10.55,
  ymin = 6.0,
  ymax = 7.8,
  label = paste0(
    "ORIGIN\n",
    "Independent"
  ),
  fill = "#E7F0F5",
  colour = COL_PRIMARY,
  linewidth = 0.75,
  text_colour = COL_PRIMARY,
  fontface = "bold",
  text_size = 3.2
)

# Link the exposure definition to the subsequent functional transition.
pB <- add_arrow(
  pB,
  x = 7.45,
  xend = 7.98,
  y = 7.2,
  yend = 7.2,
  colour = COL_PRIMARY,
  linewidth = 0.72
)


pB <- add_box(
  pB,
  xmin = 10.80,
  xmax = 12.0,
  ymin = 8.0,
  ymax = 9.0,
  label = "Independent",
  fill = "white",
  text_size = 2.65,
  fontface = "bold"
)

pB <- add_box(
  pB,
  xmin = 10.80,
  xmax = 12.0,
  ymin = 6.2,
  ymax = 7.2,
  label = "Mild",
  fill = COL_BOX,
  text_size = 2.8,
  fontface = "bold"
)

pB <- add_box(
  pB,
  xmin = 10.80,
  xmax = 12.0,
  ymin = 4.4,
  ymax = 5.4,
  label = "Severe",
  fill = COL_SEVERE,
  text_size = 2.8,
  fontface = "bold"
)

pB <- add_arrow(
  pB,
  x = 10.55,
  xend = 10.75,
  y = 7.1,
  yend = 8.45
)

pB <- add_arrow(
  pB,
  x = 10.55,
  xend = 10.75,
  y = 6.9,
  yend = 6.7
)

pB <- add_arrow(
  pB,
  x = 10.55,
  xend = 10.75,
  y = 6.7,
  yend = 4.9
)

pB <- pB +
  annotate(
    "text",
    x = 6.0,
    y = 2.2,
    label = paste0(
      "Primary: any limitation vs Independent\n",
      "Severity: Independent / Mild / Severe\n",
      "Formal contrast: Severe vs Mild"
    ),
    family = BASE_FAMILY,
    size = 3.05,
    colour = COL_TEXT,
    lineheight = 1.07
  ) +
  annotate(
    "text",
    x = 6.0,
    y = 1.0,
    label = "All modeled transition intervals begin in functional independence.",
    family = BASE_FAMILY,
    size = 2.85,
    colour = COL_LIGHT
  )


# ==============================================================================
# 9. PANEL C -- FOUR-COHORT PRIMARY ANALYSIS AND SYNTHESIS
# ==============================================================================

pC <- blank_panel(
  xmax = 12,
  ymax = 10,
  title = "C  Four-cohort analysis and synthesis"
)

get_cohort_label <- function(coh) {

  z <- sample_summary[
    cohort ==
      coh
  ]

  paste0(
    "Intervals: ",
    format(
      z$interval_N,
      big.mark = ","
    ),
    "\n",
    "Participants: ",
    format(
      z$person_N,
      big.mark = ","
    )
  )
}

card_pos <- list(
  HRS = c(
    0.25,
    2.75
  ),
  ELSA = c(
    3.05,
    5.55
  ),
  SHARE = c(
    5.85,
    8.35
  ),
  CHARLS = c(
    8.65,
    11.15
  )
)

for (coh in COHORTS) {

  xr <- card_pos[[coh]]

  pC <- add_box(
    pC,
    xmin = xr[1],
    xmax = xr[2],
    ymin = 6.15,
    ymax = 8.75,
    label = get_cohort_label(
      coh
    ),
    fill = "white",
    colour = COL_BORDER,
    linewidth = 0.55,
    text_size = 3.25,
    fontface = "plain",
    lineheight = 1.08
  )

  pC <- pC +
    annotate(
      "text",
      x = mean(xr),
      y = 8.25,
      label = coh,
      family = BASE_FAMILY,
      size = 4.05,
      fontface = "bold",
      colour = COL_PRIMARY
    )
}

# Common model box.
pC <- add_box(
  pC,
  xmin = 1.05,
  xmax = 10.35,
  ymin = 3.45,
  ymax = 5.15,
  label = paste0(
    "Common adjustment set\n",
    "age + sex + analysis interval + height + BMI + education"
  ),
  fill = "#F3F6F8",
  colour = COL_BORDER,
  linewidth = 0.55,
  text_size = 3.45,
  fontface = "plain"
)

# Arrows from cohort cards to common analysis.
target_x <- c(
  2.6,
  4.5,
  6.7,
  8.7
)

for (i in seq_along(COHORTS)) {

  coh <- COHORTS[i]
  xr <- card_pos[[coh]]

  pC <- add_arrow(
    pC,
    x = mean(xr),
    xend = target_x[i],
    y = 6.15,
    yend = 5.15,
    colour = COL_LIGHT,
    linewidth = 0.5
  )
}

# Meta-analysis synthesis.
pC <- add_box(
  pC,
  xmin = 2.35,
  xmax = 9.05,
  ymin = 0.9,
  ymax = 2.55,
  label = paste0(
    "Pooled random-effects synthesis\n",
    "REML + modified Hartung-Knapp"
  ),
  fill = COL_META,
  colour = COL_PRIMARY,
  linewidth = 0.75,
  text_colour = COL_PRIMARY,
  fontface = "bold",
  text_size = 3.60
)

pC <- add_arrow(
  pC,
  x = 5.7,
  xend = 5.7,
  y = 3.45,
  yend = 2.6,
  colour = COL_PRIMARY,
  linewidth = 0.65
)

pC <- pC +
  annotate(
    "text",
    x = 5.7,
    y = 9.4,
    label = "Locked primary common-core complete-case samples",
    family = BASE_FAMILY,
    size = 2.95,
    colour = COL_LIGHT
  )


# ==============================================================================
# 10. ASSEMBLY
# ==============================================================================
#
# Layout:
#   Row 1: Panel A across full width
#   Row 2: Panel B (58%) + Panel C (42%)
# ==============================================================================

draw_fig1_final <- function() {

  grid.newpage()

  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = 2,
        heights = unit(
          c(
            0.54,
            0.46
          ),
          "npc"
        ),
        widths = unit(
          c(
            0.53,
            0.47
          ),
          "npc"
        )
      )
    )
  )

  print(
    pA,
    vp = viewport(
      layout.pos.row = 1,
      layout.pos.col = 1:2
    )
  )

  print(
    pB,
    vp = viewport(
      layout.pos.row = 2,
      layout.pos.col = 1
    )
  )

  print(
    pC,
    vp = viewport(
      layout.pos.row = 2,
      layout.pos.col = 2
    )
  )

  upViewport()
}


# ==============================================================================
# 11. EXPORT
# ==============================================================================

PDF_FILE <- file.path(
  FIG_DIR,
  "Fig1_Study_Design_FINAL.pdf"
)

TIFF_FILE <- file.path(
  FIG_DIR,
  "Fig1_Study_Design_FINAL.tiff"
)

PNG_FILE <- file.path(
  FIG_DIR,
  "Fig1_Study_Design_FINAL.png"
)

if (capabilities("cairo")) {

  cairo_pdf(
    filename = PDF_FILE,
    width = 13.5,
    height = 9.2,
    family = BASE_FAMILY
  )

} else {

  pdf(
    file = PDF_FILE,
    width = 13.5,
    height = 9.2,
    family = BASE_FAMILY,
    useDingbats = FALSE
  )
}

draw_fig1_final()
dev.off()

tiff(
  filename = TIFF_FILE,
  width = 13.5,
  height = 9.2,
  units = "in",
  res = 600,
  compression = "lzw",
  bg = "white"
)

draw_fig1_final()
dev.off()

png(
  filename = PNG_FILE,
  width = 13.5,
  height = 9.2,
  units = "in",
  res = 300,
  bg = "white"
)

draw_fig1_final()
dev.off()


# ==============================================================================
# 12. CAPTION
# ==============================================================================

caption <- paste0(
  "Figure 1. Study timeline, exposure construction, functional-transition ",
  "framework, and four-cohort analytical synthesis. ",
  "(A) Each analytical window used two peak expiratory flow (PEF) assessments, ",
  "with the second PEF assessment defining the index for subsequent functional ",
  "follow-up. Functional state was assessed prospectively after the index. ",
  "(B) For each of the two PEF assessments, a low-PEF deficit was defined as ",
  "max(0, -residual z), where residual z was obtained from the cohort-specific ",
  "two-wave pooled PEF model adjusted for age, sex, and height. The two deficits ",
  "were averaged and then standardized within cohort among unique participants ",
  "in the locked primary sample to define the primary continuous burden exposure. ",
  "The number of low-PEF observations across the two assessments (0, 1, or 2) ",
  "was used as a supportive exposure ",
  "representation. Analyses were restricted to intervals beginning in functional ",
  "independence and evaluated subsequent Independent, Mild limitation, or Severe ",
  "limitation states, including a formal Severe-versus-Mild contrast. ",
  "(C) Cohort-specific estimates from HRS, ELSA, SHARE, and CHARLS used a common ",
  "adjustment set of age, sex, analysis interval, height, body mass index, and ",
  "education and were synthesized using random-effects meta-analysis with REML ",
  "estimation and modified Hartung-Knapp inference."
)

writeLines(
  caption,
  con = file.path(
    FIG_DIR,
    "Fig1_caption_FINAL.txt"
  ),
  useBytes = TRUE
)


# ==============================================================================
# 13. CONSOLE DASHBOARD
# ==============================================================================

cat("\n\n")
cat("================================================================================\n")
cat("FIGURE 1 FINAL COMPLETE\n")
cat("================================================================================\n")

cat("\n[A] Timeline windows:\n")
print(
  timeline_rows[
    ,
    .(
      window,
      pef1_year,
      index_year,
      followup_end
    )
  ]
)

cat("\n[C] Locked primary samples:\n")
print(
  sample_summary
)

cat("\nSaved files:\n")
cat("  PDF : ", PDF_FILE, "\n", sep = "")
cat("  TIFF: ", TIFF_FILE, "\n", sep = "")
cat("  PNG : ", PNG_FILE, "\n", sep = "")

cat("\nExposure source:\n")
cat("  Script 30 unified burden: mean[max(0, -residual_z)] across two PEF assessments,\n")
cat("  standardized within cohort among unique persons in the locked primary sample.\n")

cat("\nScientific reading order:\n")
cat("  A. When repeated PEF and subsequent functional follow-up occurred.\n")
cat("  B. How cumulative burden and functional-transition severity were defined.\n")
cat("  C. How four cohort-specific estimates were synthesized.\n")

cat("\nFINAL STATUS: PASS\n")
cat("Figure 1 generated from Script 30 locked primary samples with unified burden.\n")
cat("================================================================================\n")
