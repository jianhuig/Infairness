---
output:
  pdf_document: default
  html_document: default
---
# Semi-supervised Fairness Aduiting


# Installation

```{R, eval = FALSE}
devtools::install_github(repo = "https://github.com/jianhuig/SS-Fairness-Audit")
```

# Example
```{R}
library(SSFairnessAudit)

dat <- DataGeneration(n_labeled = 1e3,
                      N_unlabeled = 1e4,
                      prot_att_prevalence = 0.5,
                      model = "correct",
                      p = 10,
                      rho = 0.4)
                      

# Indepedent Training Dataset
indep <- DataGeneration(n_labeled = 3e3,
                      N_unlabeled = 0,
                      prot_att_prevalence = 0.5,
                      model = "correct",
                      p = 10,
                      rho = 0.4)
                      
# Train Logistic Regression Model
fit <- glm(Y ~ ., data = indep[, -3], family = binomial(link = "logit"))
dat$S <- predict(fit, newdata = dat, type = "response")
                      
# Supervised Estimation                      
sup <- Audit_Fairness(Y = dat$Y,
                      S = dat$S,
                      A = dat$A,
                      threshold = 0.5,
                      method = "supervised")
                      
sup$est
```
