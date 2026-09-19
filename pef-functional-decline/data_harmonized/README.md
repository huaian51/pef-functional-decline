# Harmonized input boundary

This directory is intentionally empty in the public repository.

The analysis code does **not** redistribute original or derived participant-level HRS, ELSA, SHARE, or CHARLS data. Users with legitimate data access should create two harmonized RDS files per cohort.

## 1. PEF exposure file

Example filename:

```text
HRS_pef_exposure.rds
```

One row per participant. Required variables:

```text
id
pef_1 age_1 sex_1 height_1
pef_2 age_2 sex_2 height_2
```

The public script pools both PEF assessments and fits:

```text
PEF ~ age + sex + height
```

No wave term is included.

PEF and height must use internally consistent cohort-specific units across the two exposure assessments. Because the model is fitted within cohort, the public code does not require the same raw PEF or height unit across different cohorts.

## 2. Functional interval file

Example filename:

```text
HRS_functional_intervals.rds
```

One row per candidate functional-transition interval.

Required:

```text
id
interval_f
destination
age
sex
height
bmi
education
```

`destination` should be one of:

```text
Independent
Mild limitation
Severe limitation
NA
```

Rows with an unobserved next-wave destination should remain in this file when possible so that the observed-follow-up IPW sensitivity can be reproduced.

Optional sensitivity variables are documented in `../metadata/variable_dictionary.csv`.

## Important coding rule for education

Preserve the cohort-harmonized model class of `education`. Do not automatically convert education to a factor in every cohort. The submitted analysis retained cohort-specific harmonized coding while using the same common-core covariate domain.

## Templates

Empty CSV templates are provided in `templates/`. They are documentation aids only; the analysis scripts read RDS files.


## 3. Strict >=2 ADL/IADL deterioration sensitivity file

The exact manuscript sensitivity for an increase of at least two ADL/IADL limitations uses a dedicated track-aware input because HRS and SHARE contain more than one eligible PEF exposure window.

Example filename:

```text
HRS_ge2_sensitivity.rds
```

One row per candidate strict-deterioration interval. Required variables:

```text
id
track
z1
z2
baseline_adl_count
baseline_iadl_count
adl_count
iadl_count
```

`z1` and `z2` are the residual-z PEF values from the same cohort-specific expected-PEF models used for the corresponding exposure track. Script 03 reconstructs the strict-outcome exposure as:

```text
deficit_1 = max(0, -z1)
deficit_2 = max(0, -z2)
burden_raw = mean(deficit_1, deficit_2)
```

For this sensitivity only, `burden_raw` is standardized within cohort among unique participant-track exposure units in the eligible baseline-independent (`ADL + IADL = 0`) risk set. This matches the final manuscript analysis.

Recommended model covariates, using the harmonized model classes from the study pipeline:

```text
age
sex
height
bmi
education
interval_f
```

`include_worsening_only` may also be supplied as a logical variable if it was part of the harmonized strict-outcome endpoint construction.

For HRS and SHARE, `track` must identify the relevant exposure window for each interval. ELSA and CHARLS may use a constant track label.

The strict-outcome files are needed to reproduce the full final sensitivity set and Supplementary Figure S3.
