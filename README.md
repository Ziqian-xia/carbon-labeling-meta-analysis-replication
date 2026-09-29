# Carbon-labeling meta-analysis replication package

Data and R code for “A meta-analysis of the real-world impact of carbon labeling on consumer choices.”

## Current dataset and manuscript targets

The canonical input is `data/metadata_change.xlsx`: **52 effect-size estimates from 22 articles**.
The 2026-09-28 extraction update corrected six directions (IDs 35, 38, 46–49),
Nowak (ID 4), and participant-level Rahmani comparisons (IDs 54–55).
The current numerical targets correspond to the manuscript/SI reviewed on **2026-09-29**.

The main multilevel estimate is d = 0.187 [0.070, 0.305]; PET is
0.074 [−0.107, 0.254]. Historical submitted results (d = 0.133; PET ≈ 0)
remain available at commit `54b059e`; they are not the current targets.
Do not use a same-named k = 55 workbook from another project directory.

Risk-of-bias ratings are linked by `Key` + `studyid` to the canonical effects and
variances in `code/00_setup.R`. The data workbook also contains the updated values.
Changing row order therefore cannot silently attach ratings to a different comparison.
The separate regression workbook must match the same canonical effect IDs and values.

## Run the numerical reproduction check

From the package root:

```sh
Rscript --vanilla run_all.R
```

Dependencies: R, readxl, metafor, dplyr, weightr, zcurve.
Verified environment: R 4.5.2; metafor 4.8.0; weightr 2.0.2; zcurve 2.4.6.
The check regenerates:

- `output/reproduction_check.csv` and `output/reproduction_summary.txt`
- `output/SI_S7.csv`–`output/SI_S10.csv`
- `output/Table_S5.csv` and `output/moderator_tests.csv`
- `output/sessionInfo.txt` and `output/input_checksums.csv`

It checks the primary and PET estimates, heterogeneity/prediction intervals,
selection-model and z-curve results, small-study slopes, setting moderator tests,
subgroup predicted means, and the five S5 omnibus p-values against manuscript targets.
It exits with an error if a checked target fails. z-curve uses `set.seed(1)`.

**Passing verifies these numerical targets, not every manuscript result or methodological assumption.**
It does not run Bayesian MCMC, certify source-paper extraction/eligibility, resolve
TES assumptions, or validate the illustrative power analysis in Figure S5.
The workbook notes retain unresolved source/definition questions for author review.

## Figures and supplementary analyses

- `code/analysis.R`: multilevel moderators and PET/PEESE.
- `code/supplementary_information.R`: SI figures and tables. Figure S2 now accumulates
  by publication year (ties by effect ID), with colors representing years. Figure S4
  uses model-predicted subgroup means and their CIs, not reference-category contrasts.
- `code/reproduce_figures_S2_S4.R`: standalone regeneration of those two corrected figures.
- `code/15_si_tables_S7_S10.R`: regenerated S7–S10 tables.
- `code/13_reviewer_response_diagnostics.R`: additional diagnostics.

S5 retains the manuscript's existing univariate mixed-effects specifications;
this update fixes its inputs and reproduction, not its design assumptions.
The historical S3 ordinary ANOVA and Figure S5 power routine remain in the source
for review; neither is certified by the numerical reproduction check.

## Bayesian analyses

`code/table.R` and `code/16_bayesian_meta_regression.R` use legacy RoBMA APIs.
`code/14_robma_bias_decomposition.R` uses RoBMA 4.0.0 and JAGS >= 4.3.1.
The PET/PEESE-only ensemble differs from selection-model ensembles.
Do not silently substitute a different API or ensemble to claim exact reproduction.
The current manuscript's Bayesian output still requires a pinned, verified environment
and archived fit/summary/convergence records; it is not covered by `run_all.R`.
