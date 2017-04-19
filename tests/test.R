library(devtools)
load_all('./R',export_all = TRUE)

# One level
set.seed(12345)
n_level_1 <- 10
n_obs <- 30

X <- runif(n = n_level_1 * n_obs , min = -1, max = 1)
g1 <- as.factor(rep(seq_len(n_level_1), each = n_obs ))

intercept_g0 <- rnorm(1)
intercept_g1 <- rnorm(n_level_1)

slope_g0 <- rnorm(1)
slope_g1 <- rnorm(n_level_1)

noise <- rnorm(length(X),sd = 0.5)
Y <- X * (slope_g0 + slope_g1[g1] ) + 
  intercept_g0 + intercept_g1[g1]  +  noise

no <- length(X)
x <- cbind(rep(1,no), X)
y <- Y
z <- list(x)
group <- list(g1)
weights = rep(1, no)
start <- etastart <- mustart <- NULL
offset = rep(0, no)
family = gaussian()
control = list(standardize = FALSE)
intercept = TRUE

fit <- mhglm.fit.multilevel(x,z,y,group,control = control)

# Two levels
set.seed(12345)
n_level_1 <- 10
n_level_2 <- 10
n_level_3 <- 10
n_obs <- 30

X <- runif(n = n_level_1 * n_obs * n_level_2 , min = -1, max = 1)
g1 <- as.factor(rep(seq_len(n_level_1), each = n_obs * n_level_2 ))
g2 <- as.factor(rep(seq_len(n_level_2 * n_level_1), each = n_obs ))

intercept_g0 <- rnorm(1)
intercept_g1 <- rnorm(n_level_1)
intercept_g2 <- rnorm(n_level_1 * n_level_2, sd = 0.5)

slope_g0 <- rnorm(1)
slope_g1 <- rnorm(n_level_1)
slope_g2 <- rnorm(n_level_1 * n_level_2, sd = sqrt(0.5))

noise <- rnorm(length(X),sd = 0.5)
Y <- X * (slope_g0 + slope_g1[g1] + slope_g2[g2] ) + 
  intercept_g0 + intercept_g1[g1] + intercept_g2[g2] +  noise

no <- length(X)
x <- cbind(rep(1,no), X)
y <- Y
z <- list(x, x)
group <- list(g1, g2)
weights = rep(1, no)
start <- etastart <- mustart <- NULL
offset = rep(0, no)
family = gaussian()
control = list(standardize = FALSE)
intercept = TRUE

fit <- mhglm.fit.multilevel(x,x,y,group = data.frame(g1,g2),control = control)
fit2 <- mhglm.fit.multilevel2(x,z,y,group,control = control)



## Three levels
set.seed(12345)
n_level_1 <- 10
n_level_2 <- 10
n_level_3 <- 10
n_obs <- 30

X <- runif(n = n_level_1 * n_obs * n_level_2 * n_level_3, min = -1, max = 1)
X2 <- runif(n = n_level_1 * n_obs * n_level_2 * n_level_3, min = -1, max = 1)
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

noise <- rnorm(length(X))
y.means <- (X + X2) * (slope_fix + slope_g1[g1] + slope_g2[g2] + slope_g3[g3]) + 
  intercept_fix + intercept_g1[g1] + intercept_g2[g2] + intercept_g3[g3] 
Y <- y.means + noise


no <- length(X)
x <- cbind(rep(1,no), X)
y <- Y
z <- list(x, x,x)
group <- list(g1, g2,g3)
weights = rep(1, no)
start <- etastart <- mustart <- NULL
offset = rep(0, no)
family = gaussian()
control = list(standardize = FALSE)
intercept = TRUE
data <- data.frame(X,X2, g1, g2, g3, Y)
 
fit <- mhglm.fit.multilevel(x,x,y,group = data.frame(factor(g1),factor(g2),factor(g3)),control = control)
fit2 <- mhglm.fit.multilevel(x,z,y,group,control = control)
ranef2 <- ebayes.est.multilevel(fit2)



## mhglm.R
data <- data.frame(X,X2, g1, g2, g3, Y)
object <- mhglm( Y ~ X + ( X|g1) + (X|g1:g2), data = data, control = list(standardize = FALSE))
object <- mhglm( Y ~ X + ( X|g1) + (X+X2|g1:g2:g3), data = data, control = list(standardize = FALSE,diagcov = TRUE))

object
pred <- predict(object,se.fit = FALSE)
pred <- predict(object,se.fit = TRUE)
fixef(object)
vcov(object)
VarCorr(object)
fitted(object)
weights(object)
residuals(object)
ranef(object,FALSE)
ranef(object,TRUE)
summary(object)

1==0
#control = list(standardize = FALSE)
#start = NULL
#model = TRUE
#method = "mhglm.fit.multilevel"
#x = FALSE 
#z = FALSE
#y = TRUE
#group = TRUE
#contrasts = NULL
#
#object <- mhglm( Y ~ X + ( X|g1) + (X|g2), data = data, control = list(standardize = FALSE))
#ranef.mhglm(object)
#
#fit <- mhglm( Y ~ X + X2 +  ( X + X2|g1) + (X|g1:g2) + (X|g1:g2:g3), data = data, control = list(standardize = FALSE))
#mf <- tmp(Y ~ X + ( X|g1) + (X|g1:g2), data = data, control = list(standardize = FALSE))
#
#fit <- mhglm.fit.multilevel(x, z= list(x,x),y = y,
#
##mf <- tmp(Y ~ X + ( X|g1) + (X|g1:g2), data = data, control = list(standardize = FALSE))
#



# 
