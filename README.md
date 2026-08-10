# Carbon-labeling meta-analysis replication package

This repository contains data and R code for reproducing the analyses in
“A meta-analysis of the real-world impact of carbon labeling on consumer choices.”

The canonical dataset is `data/metadata_change.xlsx` (k = 52 effect sizes, 22 articles).
Reproduces the manuscript headline numbers exactly: multilevel d = 0.133 [0.029, 0.236];
Q(51) = 260.77; PET-corrected g = 0.000 [−0.142, 0.142].

> NOTE: a same-named file `plot_v1108/data/metadata_change.xlsx` elsewhere is a DIFFERENT
> (k = 55) cut and must not be used. Always use the copy in this package.

## Main analysis scripts
- `code/analysis.R`   — FDR moderator tests, PET/PEESE, heterogeneity
- `code/table.R`      — RoBMA (NOTE: `priors_bias` = PET + PEESE only; no selection models), summary tables
- `code/plot.R`, `code/supplementary_information.R` — figures / SI

## Additional diagnostic analyses
| Script | Produces | Notes |
|---|---|---|
| `code/13_reviewer_response_diagnostics.R` | prediction interval, moderated PET, 3PSM, caliper, z-curve | deterministic + z-curve seed |
| `code/14_robma_bias_decomposition.R` | RoBMA bias-inclusion BF under 3 bias ensembles (PET/PEESE-only vs default-with-selection vs selection-only) | RoBMA 4.0.0 API; run with `parallel = FALSE` |
| `code/15_si_tables_S7_S10.R` | Supplementary Tables S7–S10 (CSV in `output/`) | z-curve `set.seed(1)` |
| `output/SI_S7.csv` to `output/SI_S10.csv` | table contents | regenerate by running script 15 |
| `output/reproduction_check.csv`, `output/reproduction_summary.txt` | automated consistency checks against the submitted manuscript | regenerate by running `run_all.R` |

## One-command check
Run this from the package root:

```r
Rscript --vanilla run_all.R
```

This loads `data/metadata_change.xlsx` and `data/data_regression_rob.xlsx`, regenerates
Supplementary Tables S7-S10, and writes:

- `output/reproduction_check.csv`
- `output/reproduction_summary.txt`

All automated checks should pass.

## Reproducibility / environment
- R 4.5.2; metafor 4.8.0; weightr 2.0.2; zcurve 2.4.6.
- z-curve uses bootstrap CIs → **always `set.seed(1)`** (point estimates are stable; CIs vary without a seed).
- RoBMA: requires JAGS ≥ 4.3.1. The **submitted RoBMA used RoBMA 2.x**; CRAN is now 4.0.0 (breaking API change:
  `d/v/priors_bias` → `yi/vi/measure/prior_bias`). To reproduce the submitted bias BF exactly, pin RoBMA 2.x;
  script 14 is written for 4.0.0. On Apple Silicon: `brew install jags`, then reinstall `rjags` and `RoBMA`
  from source with `PKG_CONFIG_PATH=/opt/homebrew/lib/pkgconfig`, and run with
  `DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib` and `parallel = FALSE`.

## Key findings reproduced by this package
- 3PSM (heterogeneity-robust selection model): no detectable selection; likelihood-ratio tests are non-significant across cut-points (smallest p = 0.196).
- z-curve: ERR = 0.77, EDR = 0.73, implied false-discovery rate ≈ 2% → significant findings mostly genuine.
- Prediction interval (overall): [−0.35, 0.63]; I² ≈ 96%.
- Small-study slope survives adjustment for design moderators and an SMD-artifact-robust √(1/n) predictor;
  precision uncorrelated with lab/field (p = 0.61) → asymmetry real but origin (selection vs heterogeneity) undetermined.
- RoBMA Fit A (PET/PEESE-only, submitted spec) reproduces bias inclusion BF ≈ 12.1-12.5,
  confirming the submitted bias evidence is PET/PEESE-based, not selection-based.
