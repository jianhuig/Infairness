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
  Metric     Group0    Group1       Delta
1    TPR 0.65950226 0.4003580  0.25914427
2    TNR 0.96246287 0.8650402  0.09742264
3    FPR 0.03753713 0.1349598 -0.09742264
4    FNR 0.34049774 0.5996420 -0.25914427
5    NPV 0.85549688 0.7683264  0.08717046
6    PPV 0.89348659 0.5633921  0.33009448
7    ACC 0.86455858 0.7241816  0.14037699
8     F1 0.75886756 0.4680851  0.29078245
9     BS 0.10163888 0.1894139 -0.08777504
```
