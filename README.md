# SSFairnessAudit

`SSFairnessAudit` is an R package for auditing group fairness metrics when labels are fully observed or only partially available. The repository contains the package source, generated documentation, and a small set of research scripts used for simulation experiments.

## Main workflow

Use `Audit_Fairness()` as the high-level entry point or call `SSFairness()`
directly when you want more control over the semi-supervised estimator.

`SSFairness()` now supports two useful switches:

- `cross_fit_variance = TRUE` to use the cross-fitted imputation path for variance estimation
- `return_imputation_quality = TRUE` to return imputation diagnostics plus the labeled and unlabeled imputations

## Repository layout

- `R/`: package source code
- `man/`: generated `.Rd` documentation
- `docs/`: pkgdown site output
- `scripts/`: simulation and exploratory analysis scripts

## Main user-facing functions

- `Audit_Fairness()`: wrapper for supervised and semi-supervised auditing
- `SupervisedFairness()`: fairness estimation with labeled outcomes
- `SSFairness()`: semi-supervised fairness estimation and optional imputation diagnostics
- `DataGeneration()`: synthetic data generator for simulations
- `ImputeQuality()`: compatibility wrapper for imputation-quality assessment
- `Select_Model()`: candidate-model selection helper

## Notes

- `docs/` is committed so the package website can be served from GitHub Pages.
- `scripts/` is intentionally excluded from package builds because it supports experiments rather than the package API.
