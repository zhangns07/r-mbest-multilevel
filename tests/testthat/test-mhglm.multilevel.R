
context("mhglm")

if(1==0){
  library(testthat)
  library(devtools)
  load_all('../R',export_all = TRUE)
}
library("lme4")

test_that("succeeds on sleepstudy with two levels", {
    set.seed(0)
    model.fit <- mhglm.fit.multilevel(x = as.matrix(sleepstudy[,'Days',drop = FALSE]),
				      z = as.matrix(sleepstudy[,'Days',drop = FALSE]),
				      y = as.matrix(sleepstudy[,'Reaction']),
				      group = cbind.data.frame(sleepstudy[,'Subject'],
							       factor( sample(c(1:2),180,replace = TRUE))),
				      control = list(standardize = FALSE),
				      family = gaussian())

    # fixef
    fixef0 <- c("Days" = 46.3)
    expect_that(round(model.fit$fit$coefficient.mean,1), equals(fixef0))

    # vcov
    vcov0 <- matrix(c(14.9), 1, 1)
    rownames(vcov0) <- colnames(vcov0) <- c("Days")
    expect_that(round(model.fit$fit$coefficient.mean.cov , 1), equals(vcov0))

    # VarCorr
    varcor1 <- matrix(c(269), 1, 1)
    rownames(varcor1) <- colnames(varcor1) <- c("Days")
    varcor2 <- matrix(c(4.8), 1, 1)
    rownames(varcor2) <- colnames(varcor2) <- c("Days")
    expect_that(round(model.fit$coef.cov.all[[1]] , 1), equals(varcor1))
    expect_that(round(model.fit$coef.cov.all[[2]] , 1), equals(varcor2))

    # ranef
    ranef.est  <- ebayes.est.multilevel(model.fit)
    ranef1 <- matrix(c( 10.9, -8.9, -6.2, 1.9, 3.1, 3.9, 4.9, 3, -5.8, 14.3, 
		       0.9, 6.9, 1.1, 8.5, 4.1, 3.9, 2.3, 5.6),
		     18,1)
    rownames(ranef1) <- as.character(c(308, 309, 310, 330, 331, 332, 333, 334,
                                       335, 337, 349, 350, 351, 352, 369, 370,
                                       371, 372))
    expect_that(round(ranef.est[[1]],1), equals(ranef1))

    ranef2 <- matrix(c( 0.2, 0, -0.1, -0.1, 0, -0.1, 0.1, -0.1, -0.2, 
		       0.2, 0.5, -0.4, 0.1, 0, 0.1, -0.1, -0.2, 0.1, 
		       0.1, 0.2, 0, 0, 0, 0.1, 0, 0, 0.2, 
		       -0.1, 0, 0.1, -0.1, 0.2, 0, 0.1, 0.1, 0),
		       36,1)
    rownames(ranef2) <- c("308:1", "308:2", "309:1", "309:2", "310:1", "310:2", "330:1", "330:2", "331:1", 
			  "331:2", "332:1", "332:2", "333:1", "333:2", "334:1", "334:2", "335:1", "335:2", 
			  "337:1", "337:2", "349:1", "349:2", "350:1", "350:2", "351:1", "351:2", "352:1",
			  "352:2", "369:1", "369:2", "370:1", "370:2", "371:1", "371:2", "372:1", "372:2")

    expect_that(round(ranef.est[[2]],1), equals(ranef2))

})


