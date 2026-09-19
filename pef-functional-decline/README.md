# Cumulative low PEF burden and subsequent functional decline

Public analysis code accompanying the manuscript:

> **Cumulative low peak expiratory flow burden and severity of subsequent functional decline: a four-cohort longitudinal study**

**Authors:** Wei Jiang; Xiaohua Dai (corresponding author)

This repository contains a clean, public-facing implementation of the final analytical pathway reported in the manuscript. It is designed for transparency, review, and reuse by researchers who have legitimate access to the source cohort data.

## What this repository does

The code implements the final cross-cohort analysis in four longitudinal ageing studies:

- Health and Retirement Study (HRS)
- English Longitudinal Study of Ageing (ELSA)
- Survey of Health, Ageing and Retirement in Europe (SHARE)
- China Health and Retirement Longitudinal Study (CHARLS)

The main workflow:

1. constructs a harmonized repeated low-PEF exposure;
2. creates the primary complete-case analysis datasets;
3. estimates loss of functional independence and destination severity;
4. pools cohort-specific estimates using REML random-effects meta-analysis with modified Hartung-Knapp inference;
5. performs the key supportive and sensitivity analyses reported in the final manuscript;
6. creates machine-readable summary tables and forest plots.

## Unified exposure algorithm

Within each cohort, the two designated PEF exposure assessments are pooled and the same expected-PEF model is fitted:

```text
PEF ~ age + sex + height
```

No wave term is included.

For each PEF assessment:

```text
residual_z = model residual / pooled residual SD
deficit_t  = max(0, -residual_z_t)
```

The raw repeated burden is:

```text
burden_raw = mean(deficit_1, deficit_2)
```

The burden is then standardized within cohort among unique participants in the locked primary complete-case sample. Main associations are reported per 1-SD higher burden.

A supportive threshold exposure counts the number of PEF assessments with residual z <= -1 SD (0, 1, or 2 low-PEF observations). This is a frequency measure, **not a trajectory**.

## Functional outcomes

All primary functional-transition intervals begin in **functional independence**. The subsequent destination is classified as:

- **Independent:** 0 ADL/IADL limitations
- **Mild limitation:** 1 limitation
- **Severe limitation:** >=2 limitations

The primary estimands are:

- Any limitation vs Independent
- Mild vs Independent
- Severe vs Independent
- Severe vs Mild among limitation destinations

The Severe-vs-Mild comparison does **not** represent progression from a Mild state to a Severe state. It compares the severity of the subsequent destination among intervals that end in Mild or Severe limitation.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── 00_run_all.R
├── 01_build_analysis_datasets.R
├── 02_primary_analysis.R
├── 03_sensitivity_analysis.R
├── 04_make_tables_figures.R
├── R/
│   └── helpers.R
├── metadata/
│   ├── expected_primary_results.csv
│   ├── expected_sensitivity_results.csv
│   ├── expected_sample_counts.csv
│   ├── study_windows.csv
│   └── variable_dictionary.csv
├── data_harmonized/
│   ├── README.md
│   └── templates/
├── examples/
│   └── 00_make_synthetic_inputs.R
└── output/
```

## Data are not redistributed

The source cohort data are **not** included in this repository. Researchers must obtain HRS, ELSA, SHARE, and CHARLS data from the corresponding repositories and comply with each study's access and data-use terms.

The public analysis layer starts from harmonized cohort extracts with a documented schema. This boundary is intentional: the repository shares the analytical method without redistributing restricted participant-level data.

See [`data_harmonized/README.md`](data_harmonized/README.md) and [`metadata/variable_dictionary.csv`](metadata/variable_dictionary.csv).

## Input files

For each cohort, place two RDS files in `data_harmonized/`:

```text
HRS_pef_exposure.rds
HRS_functional_intervals.rds
ELSA_pef_exposure.rds
ELSA_functional_intervals.rds
SHARE_pef_exposure.rds
SHARE_functional_intervals.rds
CHARLS_pef_exposure.rds
CHARLS_functional_intervals.rds
```

The first contains the two designated PEF assessments and the covariates needed to fit `PEF ~ age + sex + height`.

The second contains the harmonized functional-transition intervals, primary common-core covariates, and—where available—variables used in sensitivity analyses.

## Required software

The analyses were developed for R 4.5.x. Core packages:

```r
install.packages(c("data.table", "MASS", "nnet", "ggplot2"))
```

`ggplot2` is required only for the plotting stage.

## Execution

Run from the repository root:

```r
source("01_build_analysis_datasets.R")
source("02_primary_analysis.R")
source("03_sensitivity_analysis.R")
source("04_make_tables_figures.R")
```

or simply:

```r
source("00_run_all.R")
```

Outputs are written to:

```text
output/derived_data/
output/tables/
output/QC/
output/figures/
```

## Primary inference

### Binary outcomes

Cohort-specific logistic regression with participant-cluster robust standard errors.

### Multinomial severity outcome

Multinomial logistic regression with Independent as the reference category. Inference uses 1,000 participant-level bootstrap resamples. Participants are sampled with replacement and all intervals contributed by a sampled participant are retained through cluster multiplicity weights.

### Meta-analysis

Cohort-specific log odds ratios are pooled using:

- REML estimation of between-cohort variance;
- modified Hartung-Knapp uncertainty;
- I² as a descriptive measure of between-cohort heterogeneity.

Because only four cohorts contribute, heterogeneity estimates should not be overinterpreted.

## Sensitivity analyses

The public code includes:

- restriction to the smoking-observed sample;
- current-smoking adjustment where available;
- broader cohort-specific adjustment using variables with >=85% coverage and usable variation;
- 0/1/2 low-PEF frequency analysis;
- Severe limitation vs all non-severe destinations;
- stabilized inverse-probability weighting for observed next-wave functional status;
- strict >=2 ADL/IADL-item deterioration where paired counts are available;
- leave-one-cohort-out meta-analysis.

The IPW sensitivity addresses measured predictors of follow-up observation under the fitted model; it should not be interpreted as eliminating all selection bias.

## Locked reference results

`metadata/expected_primary_results.csv` and `metadata/expected_sensitivity_results.csv` record the final manuscript estimates as **audit targets only**. They are never used in model fitting.

The primary pooled odds ratios in the locked analysis were approximately:

| Estimand | Pooled OR |
|---|---:|
| Any limitation vs Independent | 1.159 |
| Mild vs Independent | 1.090 |
| Severe vs Independent | 1.270 |
| Severe vs Mild | 1.172 |

Exact agreement can depend on using the same cohort releases, harmonized inputs, variable coding, and bootstrap seed.

## Synthetic smoke test

To inspect the workflow without any cohort data:

```r
source("examples/00_make_synthetic_inputs.R")
source("01_build_analysis_datasets.R")
source("02_primary_analysis.R")
source("03_sensitivity_analysis.R")
```

The synthetic files are generated solely to demonstrate the main input schema and core code execution. They are not derived from any study participant and have no scientific meaning. The locked manuscript figures in Script 04 intentionally enforce the real manuscript sample/result gates and therefore are not part of the synthetic smoke test.

## Reproducibility scope

This is a **public release of the final analysis pathway**, not an archive of every exploratory script produced during study development.

Not included:

- exploratory analyses not retained in the manuscript;
- discarded model specifications;
- debugging scripts;
- private local paths;
- restricted cohort source data;
- intermediate data files containing participant-level information.

This repository should therefore be described as:

> **Public analysis code for the primary statistical models and key sensitivity analyses.**

## Interpretation

This is an observational study. The code estimates associations and does not establish that lower PEF causes later functional decline.

PEF is treated as a repeated physiological performance marker with broad pulmonary, muscular, neuromuscular, anthropometric, and health-related determinants. It should not be interpreted as a disease-specific mechanism or as equivalent to a multidimensional intrinsic-capacity measure.

## License

The code in this repository is released under the MIT License. The license applies to the code only and does not grant rights to any HRS, ELSA, SHARE, or CHARLS data.

## Contact

**Wei Jiang**  
Anhui University of Chinese Medicine, Hefei, Anhui Province, China

**Xiaohua Dai, corresponding author**  
The First Affiliated Hospital of Anhui University of Chinese Medicine, Hefei, Anhui Province, China


## Citation and permanent archive

The source repository may be hosted on GitHub. Versioned releases can be archived in Zenodo to provide a permanent DOI for the exact software version associated with the manuscript.

This repository includes a `CITATION.cff` file so that citation metadata can be recognized by GitHub and software-archiving services.


## Version

This repository is prepared as **v1.0.0**, corresponding to the analytical code accompanying the submitted manuscript.
