# r-mbest
Parallel implementation of Moment-Based Estimation for Hierarchical Models (R Package)

DEMO: fit three-levels multilevel GLM. 

```
## generate data
set.seed(12345)
n_level_1 <- 10
n_level_2 <- 10
n_level_3 <- 10
n_obs <- 30

x <- runif(n = n_level_1 * n_obs * n_level_2 * n_level_3, min = -1, max = 1)
g1 <- rep(seq_len(n_level_1), each = n_obs * n_level_2 * n_level_3)
g2 <- rep(seq_len(n_level_2 * n_level_1), each = n_obs * n_level_3)
g3 <- rep(seq_len(n_level_2 * n_level_1 * n_level_3), each = n_obs)

intercept_fix <- runif(1)
intercept_g1 <- rnorm(n_level_1)
intercept_g2 <- rnorm(n_level_1 * n_level_2, sd = 0.5)
intercept_g3 <- rnorm(n_level_1 * n_level_2 * n_level_3, sd = sqrt(0.25))

slope_fix <- runif(1)
slope_g1 <- rnorm(n_level_1)
slope_g2 <- rnorm(n_level_1 * n_level_2, sd = sqrt(0.5))
slope_g3 <- rnorm(n_level_1 * n_level_2 * n_level_3, sd = 0.25)

noise <- rnorm(length(x))

## The model is y ~ 1 + x + (1 + x || g1/g2/g3),
## i.e. nested groups and independent covariates.
y <- x * (slope_fix + slope_g1[g1] + slope_g2[g2] + slope_g3[g3]) + 
  intercept_fix + intercept_g1[g1] + intercept_g2[g2] + intercept_g3[g3] + noise

## Fit using multilevel mhglm
no <- length(x)
mhglm.fit <- mhglm.fit.multilevel(x = cbind(rep(1,no),x),
                                  z = cbind(rep(1,no),x),
                                  y = y,
                                  group =  cbind.data.frame(factor(g1),factor(g2),factor(g3)),
                                  control = list(standardize = FALSE,diagcov = TRUE))

## Empirical bayes inference
coef.eb <- ebayes.est.multilevel(mhglm.fit)

## Estimated fixed effects
mhglm.fit$fit$coefficient.mean
##                   x
## 0.8035564 1.1102478

## Estimated covariance matrix for each level
mhglm.fit$coef.cov.all
## [[1]]
##                    x
##   0.9258569 0.000000
## x 0.0000000 1.056306
## 
## [[2]]
##                     x
##   0.2130053 0.0000000
## x 0.0000000 0.4533855
## 
## [[3]]
##                      x
##   0.2372565 0.00000000
## x 0.0000000 0.06529185

## Estimated random effects
head(coef.eb[[1]])
##          [,1]        [,2]
## 1  0.33337467 -0.08160486
## 2 -1.37982807 -0.17602401
## 3  1.16104275 -1.16100055
## 4  0.05734561 -1.61044158
## 5 -0.20181388  0.64579372
## 6 -1.49864456 -0.34370688

head(coef.eb[[2]])
##           [,1]        [,2]
## 1:1 0.03285750  1.18729870
## 1:2 0.46912947 -0.01926404
## 1:3 0.09913568 -0.86173770
## 1:4 0.29749991 -0.13917048
## 1:5 0.04359381 -0.06830378
## 1:6 0.10773847 -0.10646400

head(coef.eb[[3]])
##              [,1]        [,2]
## 1:1:1 -0.57773686  0.11023873
## 1:1:2  0.32047161  0.06185304
## 1:1:3 -0.03093421 -0.26831903
## 1:1:4 -0.02303099  0.19919693
## 1:1:5  0.32545163 -0.15849075
## 1:1:6 -0.14832079 -0.02758480



```
