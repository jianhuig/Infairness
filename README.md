# SSFairnessAudit

`SSFairnessAudit` is an R package for auditing group fairness metrics when labels are fully observed or only partially available. The repository contains the package source, generated documentation, and a small set of research scripts used for simulation experiments.

## Repository layout

- `R/`: package source code
- `man/`: generated `.Rd` documentation
- `docs/`: pkgdown site output
- `scripts/`: simulation and exploratory analysis scripts

## Main user-facing functions

- `Audit_Fairness()`: wrapper for supervised and semi-supervised auditing
- `SupervisedFairness()`: fairness estimation with labeled outcomes
- `SSFairness()`: semi-supervised fairness estimation
- `DataGeneration()`: synthetic data generator for simulations
- `ImputeQuality()`: imputation-quality assessment by group
- `Select_Model()`: candidate-model selection helper

## Notes

- `docs/` is committed so the package website can be served from GitHub Pages.
- `scripts/` is intentionally excluded from package builds because it supports experiments rather than the package API.
