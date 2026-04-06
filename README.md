# SSFairnessAudit

`SSFairnessAudit` is an R package for auditing group fairness metrics when labels are fully observed or only partially available. The repository contains the package source, generated documentation, and a small set of research scripts used for simulation experiments.

## Main workflow

Use `Audit_Fairness()` as the high-level entry point or call `SSFairness()`
directly when you want more control over the semi-supervised estimator.

`SSFairness()` now supports useful controls:

- `cross_fit_variance = TRUE` to use the cross-fitted imputation path for variance estimation
- `return_imputation_quality = TRUE` to return imputation diagnostics plus the labeled and unlabeled imputations
- `folds = ...` to reuse the same labeled-data folds across candidate models when comparing them with `Select_Model()`

The semi-supervised basis options now include polynomial, natural spline,
interaction, beta-calibration, and kernel branches. The natural spline path is
available through `basis = "Spline(S)"`, `basis = "Spline(S) + X"`, and
`basis = "Spline Interaction"`. The additive spline branch uses a shared
smooth in `S` plus additive covariate effects; the spline interaction branch
adds spline-by-covariate interactions so the shape in `S` can vary with `X`.

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
- `Select_Model()`: candidate-model selection helper

When comparing candidate semi-supervised models with cross-fitted imputation
quality, reuse the same `folds` object across all `SSFairness()` calls so the
comparison is based on the same labeled-data splits.

## Notes

- `docs/` is committed so the package website can be served from GitHub Pages.
- `scripts/` is intentionally excluded from package builds because it supports experiments rather than the package API.
