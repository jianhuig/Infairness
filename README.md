# 📦 Infairness: Reliable Fairness Auditing with Semi-Supervised Inference

This R package provides tools for evaluating group fairness metrics using both supervised and semi-supervised estimators. It implements the methodology described in our paper:  
**“Reliable Fairness Auditing with Semi-Supervised Inference.”**

---

## 🚀 Installation

To install the development version directly from GitHub:

```r
devtools::install_github("jianhuig/Infairness")
```

---

## 📘 Example Usage

Here is a minimal example using simulated data to compare supervised and semi-supervised fairness estimates:

```r
library(Infairness)
library(dplyr)

# Generate semi-supervised dataset
set.seed(123)
dat <- DataGeneration(n_labeled = 1e3,
                      N_unlabeled = 1e4,
                      prot_att_prevalence = 0.5,
                      model = "scenario 1",
                      rho = 0.4)

# Independent training data for prediction model
indep <- DataGeneration(n_labeled = 3e3,
                        N_unlabeled = 0,
                        prot_att_prevalence = 0.5,
                        model = "scenario 1",
                        rho = 0.4)

# Fit group-specific logistic regression
model_0 <- glm(Y ~ ., family = binomial(), data = indep %>% filter(A == 0) %>% select(-Y_miss, -A, -starts_with("W")))
model_1 <- glm(Y ~ ., family = binomial(), data = indep %>% filter(A == 1) %>% select(-Y_miss, -A, -starts_with("W")))

# Predict risk scores
dat$S <- NA
dat$S[dat$A == 0] <- predict(model_0, newdata = dat %>% filter(A == 0), type = "response")
dat$S[dat$A == 1] <- predict(model_1, newdata = dat %>% filter(A == 1), type = "response")

# Supervised Estimation
sup <- Audit_Fairness(Y = dat$Y_miss,
                      S = dat$S,
                      A = dat$A,
                      threshold = 0.5,
                      method = "supervised")

# Semi-Supervised Estimation
ss <- Audit_Fairness(Y = dat$Y_miss,
                     S = dat$S,
                     A = dat$A,
                     threshold = 0.5,
                     method = "Infairness",
                     W = dat %>% select(contains("W")))
```

---

## 📊 Output: Point Estimates

### Supervised

| Metric | Group 0 | Group 1 | Delta |
|--------|---------|---------|--------|
| TPR    | 0.3629  | 0.5294  | -0.1665 |
| FPR    | 0.0485  | 0.0603  | -0.0118 |
| PPV    | 0.7031  | 0.7412  | -0.0381 |
| NPV    | 0.8252  | 0.8596  | -0.0344 |
| F1     | 0.4787  | 0.6176  | -0.1389 |
| ACC    | 0.8101  | 0.8388  | -0.0288 |
| BS     | 0.1325  | 0.1168  |  0.0157 |

### Semi-Supervised

| Metric | Group 0 | Group 1 | Delta |
|--------|---------|---------|--------|
| TPR    | 0.4040  | 0.5146  | -0.1105 |
| FPR    | 0.0594  | 0.0490  |  0.0104 |
| PPV    | 0.6986  | 0.7717  | -0.0731 |
| NPV    | 0.8224  | 0.8589  | -0.0365 |
| F1     | 0.5120  | 0.6174  | -0.1055 |
| ACC    | 0.8042  | 0.8448  | -0.0405 |
| BS     | 0.1355  | 0.1159  |  0.0196 |

---

## 📈 Output: Variance Estimates

### Supervised

| Metric | Group 0     | Group 1     | Delta      |
|--------|-------------|-------------|------------|
| TPR    | 1.86e-03    | 2.09e-03    | 3.96e-03   |
| FPR    | 1.18e-04    | 1.55e-04    | 2.73e-04   |
| PPV    | 3.26e-03    | 2.26e-03    | 5.52e-03   |
| NPV    | 3.19e-04    | 3.02e-04    | 6.21e-04   |
| F1     | 2.02e-03    | 1.60e-03    | 3.62e-03   |
| ACC    | 2.98e-04    | 2.79e-04    | 5.77e-04   |
| BS     | 9.94e-05    | 8.95e-05    | 1.89e-04   |

### Semi-Supervised

| Metric | Group 0     | Group 1     | Delta      |
|--------|-------------|-------------|------------|
| TPR    | 8.66e-04    | 1.04e-03    | 1.91e-03   |
| FPR    | 7.27e-05    | 1.04e-04    | 1.77e-04   |
| PPV    | 3.03e-03    | 2.10e-03    | 5.13e-03   |
| NPV    | 2.67e-04    | 2.64e-04    | 5.31e-04   |
| F1     | 1.18e-03    | 1.08e-03    | 2.27e-03   |
| ACC    | 2.52e-04    | 2.44e-04    | 4.96e-04   |
| BS     | 8.44e-05    | 7.62e-05    | 1.61e-04   |

---

## 📬 Contact

For questions, suggestions, or contributions, please contact: jianhui.gao@mail.utoronto.ca
