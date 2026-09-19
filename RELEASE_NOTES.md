# Release notes

## v1.0.0 — 2026-09-19

First public release accompanying the manuscript:

*Cumulative low peak expiratory flow burden and severity of subsequent functional decline: a four-cohort longitudinal study*


### Exact manuscript-alignment update

The public `03_sensitivity_analysis.R` and `04_make_tables_figures.R` were rebuilt against the locked final manuscript pipeline.

- The 0/1/2 low-PEF analysis now uses the full eligible interval risk set rather than the continuous-burden primary sample.
- CHARLS ever-smoking is a dedicated sensitivity analysis and is not folded into the generic extended-covariate selection.
- The strict >=2 ADL/IADL deterioration analysis uses track-specific residual-z PEF inputs and re-standardizes the unified burden within its own eligible unique person-track risk set.
- Script 04 uses public-path copies of the locked final Figure 1-3 and Figure S1-S3 plotting scripts and reconstructs the final Table 1 and Table 2 definitions.

### Included

- unified repeated low-PEF exposure construction;
- primary complete-case analysis;
- participant-cluster robust binary models;
- multinomial models with participant-level bootstrap inference;
- Severe-versus-Mild formal contrast;
- REML random-effects meta-analysis with modified Hartung-Knapp inference;
- smoking and extended-adjustment sensitivity analyses;
- 0/1/2 low-PEF frequency analysis;
- alternative severe outcome analysis;
- stabilized inverse-probability weighting;
- strict >=2 ADL/IADL deterioration analysis;
- leave-one-cohort-out analysis;
- input schemas and synthetic examples;
- locked manuscript estimates used only as reproducibility audit targets.

### Not included

- restricted participant-level HRS, ELSA, SHARE, or CHARLS data;
- exploratory or discarded analyses;
- local development paths;
- private intermediate datasets.

The public analysis layer begins from harmonized cohort inputs documented in `data_harmonized/README.md` and `metadata/variable_dictionary.csv`.
