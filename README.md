# Semi-supervised Fairness Aduiting


# Installation

```{R, eval = FALSE}
devtools::install_github(repo = "https://github.com/jianhuig/SS-Fairness-Audit")
```

# Example
```{R}
library(SSFairnessAudit)
library(dplyr)
library(glmnet)

set.seed(123)

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
sup <- Audit_Fairness(Y = dat$Y_miss,
                      S = dat$S,
                      A = dat$A,
                      threshold = 0.5,
                      method = "supervised")

# Semi-supervised Estimation
ss <- Audit_Fairness(Y = dat$Y_miss,
                     S = dat$S,
                     A = dat$A,
                     threshold = 0.5,
                     method = "semi-supervised",
                     X = dat %>% select(contains("X")))
                      

# Point Estimate
                      
sup$est
  Metric     Group0     Group1       Delta
1    TPR 0.67515924 0.34246575  0.33269348
2    TNR 0.95965418 0.90857143  0.05108275
3    FPR 0.04034582 0.09142857 -0.05108275
4    FNR 0.32484076 0.65753425 -0.33269348
5    NPV 0.86718750 0.76811594  0.09907156
6    PPV 0.88333333 0.60975610  0.27357724
7    ACC 0.87103175 0.74193548  0.12909626
8     F1 0.76534296 0.43859649  0.32674647
9     BS 0.10372732 0.17957283 -0.07584551

ss$est
  Metric     Group0    Group1       Delta
1    TPR 0.69000174 0.4217571  0.26824467
2    TNR 0.96296384 0.8827309  0.08023296
3    FPR 0.03703616 0.1172691 -0.08023296
4    FNR 0.30999826 0.5782429 -0.26824467
5    NPV 0.87127872 0.7694673  0.10181142
6    PPV 0.89528839 0.6219136  0.27337476
7    ACC 0.87709954 0.7380640  0.13903553
8     F1 0.77935330 0.5026422  0.27671109
9     BS 0.10200908 0.1818544 -0.07984536

# Variance Estimate
sup$var
  Metric       Group0       Group1        Delta
1    TPR 1.396938e-03 0.0015423490 0.0029392869
2    TNR 1.115794e-04 0.0002373411 0.0003489205
3    FPR 1.115794e-04 0.0002373411 0.0003489205
4    FNR 1.396938e-03 0.0015423490 0.0029392869
5    NPV 2.999306e-04 0.0004302267 0.0007301572
6    PPV 8.587963e-04 0.0029018732 0.0037606695
7    ACC 2.228878e-04 0.0003860226 0.0006089104
8     F1 8.004906e-04 0.0016862446 0.0024867353
9     BS 3.793972e-05 0.0001209641 0.0001589038

ss$var
  Metric       Group0       Group1        Delta
1    TPR 6.676291e-04 0.0008522866 0.0015199157
2    TNR 7.805201e-05 0.0001260279 0.0002040800
3    FPR 7.805201e-05 0.0001260279 0.0002040800
4    FNR 6.676291e-04 0.0008522866 0.0015199157
5    NPV 2.206002e-04 0.0003915897 0.0006121900
6    PPV 7.004769e-04 0.0027704193 0.0034708962
7    ACC 1.677678e-04 0.0003485351 0.0005163029
8     F1 4.533742e-04 0.0011296395 0.0015830137
9     BS 2.103213e-05 0.0001111294 0.0001321615
```
