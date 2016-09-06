# Copyright 2014 Patrick O. Perry
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# - xlevels is NULL
# - passing in group in newdata doesn't matter


mhglm.control <- function(standardize = TRUE, steps = 1, 
			  parallel = FALSE, diagcov = FALSE,
                          fit.method = "firthglm.fit",
                          fit.control = list(...), ...)
{
    if (!is.logical(standardize) || is.na(standardize))
        stop("value of 'standardize' must be TRUE or FALSE")
    if (!is.numeric(steps) || steps < 0)
        stop("number of steps must be >= 0")
    if (!is.logical(parallel) || is.na(parallel))
        stop("value of 'parallel' must be TRUE or FALSE")
    if (!is.logical(diagcov) || is.na(diagcov))
      stop("value of 'diagcov' must be TRUE or FALSE")

    if (!is.character(fit.method) && !is.function(fit.method))
        stop("invalid 'fit.method' argument")
    if (identical(fit.method, "firthglm.fit"))
        fit.control <- do.call("firthglm.control", fit.control)
    if (identical(fit.method, "glm.fit"))
        fit.control <- do.call("glm.control", fit.control)

    list(standardize = standardize, steps = steps, parallel = parallel, diagcov = diagcov,
         fit.method = fit.method, fit.control = fit.control)
}


mhglm <- function(formula, family = gaussian, data, weights, subset,
                  na.action, start = NULL, etastart, mustart, offset,
                  control = list(), model = TRUE, method = "mhglm.fit.multilevel",
                  x = FALSE, z = FALSE, y = TRUE, group = TRUE,
                  contrasts = NULL)
{
    # call
    call <- match.call()

    # family
    if (is.character(family))
        family <- get(family, mode="function", envir=parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }

    # data
    if (missing(data))
        data <- environment(formula)

    # model frame
    mf <- match.call(expand.dots=FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf$formula <- lme4::subbars(mf$formula)
    mf <- eval(mf, parent.frame())

    # method
    if (identical(method, "model.frame"))
        return(mf)
    if (!is.character(method) && !is.function(method))
        stop("invalid 'method' argument")
    if (identical(method, "mhglm.fit"))
        control <- do.call("mhglm.control", control)

    # terms
    mt <- attr(mf, "terms")

    # response
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }

    logging::loginfo("Creating design matrix", logger="mbest.mhglm")

    mt.fixed <- delete.response(terms(lme4::nobars(formula), data=data))
    X <- if (!is.empty.model(mt.fixed))
        model.matrix(mt.fixed, mf, contrasts)
    else matrix(, NROW(Y), 0L)

    logging::loginfo("Grouping factor and random effect design matrix",
                     logger="mbest.mhglm")

    bars <- lme4::findbars(formula)
    if (length(bars) >= 1L) {
      Z <- list()
      Group <- list()
      group.call <- list()
      mt.random <- list()

      for(i in seq_along(bars)){
	b <- bars[[i]]
	mf1 <- mf
	for (v in all.vars(b[[3L]])) {
	  mf1[[v]] <- factor(mf1[[v]])
	}
	group.call.tmp <- substitute(factor(fac), list(fac = b[[3L]]))
	group.tmp <- eval(group.call.tmp, mf1)

	if (all(is.na(group.tmp)))
	  stop("Invalid grouping factor specification, ", deparse(b[[3L]]))

	mt.random.tmp <- terms(eval(substitute(~trms, list(trms = b[[2L]]))),
			   data = data)
	z.tmp <- if (!is.empty.model(mt.random.tmp)){
	  model.matrix(mt.random.tmp, mf, contrasts)
	}else matrix(, NROW(Y), 0L)

	Z[[i]] <- z.tmp
	Group[[i]] <- group.tmp
	group.call[[i]] <- group.call.tmp
	mt.random[[i]] <- mt.random.tmp
      } 
    } else { # length(bars) == 0L
      Group <- factor(character(NROW(Y)))
      mt.random <- terms(~ -1, data=data)
      Z <- matrix(, NROW(Y), 0L)
    }

    logging::loginfo("Setting weights", logger="mbest.mhglm")

    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")

    logging::loginfo("Setting offset", logger="mbest.mhglm")

    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
    }

    # starting values
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")

    # group-specific estimates
    logging::loginfo("Fitting model", logger="mbest.mhglm")
    fit <- eval(call(if (is.function(method)) "method" else method,
                     x = X, z = Z, y = Y, group = Group,
                     weights = weights, start = start, etastart = etastart,
                     mustart = mustart, offset = offset, family = family,
                     control = control,
                     intercept = attr(mt.fixed, "intercept") > 0L))
    fit$contrasts.fixed <- attr(X, "contrasts")
    fit$contrasts.random <- attr(Z, "contrasts")

#    xlevels <- .getXlevels(mt.fixed, mf)
#    #xlevels.random <- .getXlevels(mt.random, mf)
#    xlevels.random <- lapply(mt.random, function(x){.getXlevels(x, mf)})
#
#    for (i in names(xlevels.random)) {
#        xlevels[[i]] <- xlevels.random[[i]]
#    }
#    fit$xlevels <- xlevels

    fit$group.levels <- lapply(Group, levels)
    fit$call <- call
    fit$control <- control
    fit$terms <- mt
    fit$terms.fixed <- mt.fixed
    fit$terms.random <- mt.random
    fit$group.call <- group.call
    if (model)
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x)
        fit$x <- X
    if (z)
        fit$z <- Z
    if (!y)
        fit$y <- NULL
    if (!group)
        fit$group <- NULL
    class(fit) <- "mhglm"

    fit
}


family.mhglm <- function(object, ...)
{
    object$family
}


terms.mhglm <- function(x, type=c("fixed", "random"), ...)
{
    type <- match.arg(type)
    switch(type, fixed = x$terms.fixed, random = x$terms.random)
}


model.frame.mhglm <- function(formula, ...)
{
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
    if (length(nargs) || is.null(formula$model)) {
        fcall <- formula$call
        fcall$method <- "model.frame"
        fcall[[1L]] <- quote(mhglm)
        fcall[names(nargs)] <- nargs
        env <- environment(formula$terms)
        if (is.null(env))
            env <- parent.frame()
        eval(fcall, env)
    }
    else formula$model
}


model.matrix.mhglm <- function(object, type=c("fixed", "random"), ...)
{
    type <- match.arg(type)

    mf <- model.frame(object)
    mt <- terms(object, type)

    if(type == 'fixed'){
      if (!is.empty.model(mt)){
	re <-  model.matrix(mt, mf, object$contrasts.fixed)
      } else {re <- matrix(, nobs(object), 0L)}
    } else if (type == 'random'){
      if(all( !unlist(lapply(mt,is.empty.model)))){
	re <-  lapply(mt,function(x){model.matrix(x, mf, object$contrasts.random)})
      } else {re <- matrix(, nobs(object), 0L)}
    } else {
      stop('Wrong type')
    }

    re
}

predict.mhglm <- function(object, newdata = NULL,
                          type = c("link", "response"),
                          se.fit = FALSE, na.action = na.pass, ...)
{
# TODO: when newdata is not NULL, things get weird
    type <- match.arg(type)

    tt.fixed <- terms(object, "fixed")
    tt.random <- terms(object, "random")

    if (missing(newdata) || is.null(newdata)) {
        x <- model.matrix(object, "fixed")
        z <- model.matrix(object, "random")
        m <- model.frame(object)
        offset <- object$offset
    } else {
        tt <- object$terms
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action,
                         xlev=object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)

        x <- model.matrix(tt.fixed, m, contrasts.arg=object$contrasts.fixed)
        z <- model.matrix(tt.random, m, contrasts.arg=object$contrasts.random)

        offset <- rep(0, nrow(x))
        if (!is.null(off.num <- attr(tt, "offset")))
            for (i in off.num)
                offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
        if (!is.null(object$call$offset))
            offset <- offset + eval(object$call$offset, newdata)
    }

    group <- list()
    for(r in seq_along(object$group.call)){
      mg <- m
      for (v in all.vars(object$group.call[[r]])) {
	mg[[v]] <- factor(mg[[v]])
      }
      group[[r]] <- eval(object$group.call[[r]],mg)
    }

    re <- ranef(object, condVar=se.fit)
    eta <- drop(x %*% fixef(object))

    for(r in seq_len(length(re))){
      group.ix <- match(group[[r]], object$group.levels[[r]])
      old <- !is.na(group.ix)
      eta[old] <- eta[old] + rowSums(z[[r]][old,,drop=FALSE]
				     * re[[r]][group.ix[old],,drop=FALSE])
    }

    if (!is.null(offset))
        eta <- eta + offset

    if (se.fit) {
        eta.se2 <- pmax(0, rowSums((x %*% vcov(object)) * x))
	var.corr <- VarCorr(object)

	for(r in seq_along(re)){
	  re.sigma <- attr(re[[r]], "postVar")
	  sigma0 <- var.corr[[r]]
	  group.ix <- match(group[[r]], object$group.levels[[r]])
	  old <- !is.na(group.ix)

	  for (i in seq_along(eta)) {
            zi <- drop(z[[r]][i,])
            sigma <- if (old[i])
                re.sigma[,,group.ix[i]]
            else sigma0

            eta.se2[i] <- eta.se2[i] + max(0, t(zi) %*% sigma %*% zi)
	  }
	}

	eta.se <- sqrt(eta.se2)
	names(eta.se) <- names(eta)
    } else {
        eta.se <- NULL
    }

    fit <- switch(type, link = eta, response =  family(object)$linkinv(eta))

    if (!se.fit) {
        pred <- fit
    } else {
        residual.scale <- as.vector(sqrt(object$dispersion))
        se.fit <- switch(type,
                         link = eta.se,
                         response = eta.se * abs(family(object)$mu.eta(eta)))
        pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
    }
    pred
}




fixef.mhglm <- function(object, ...)
{
    object$fit$coefficient.mean
}


vcov.mhglm <- function(object, ...)
{
    object$fit$coefficient.mean.cov
}


sigma.mhglm <- function(object, ...)
{
    sqrt(object$dispersion)
}


VarCorr.mhglm <- function(x, sigma=1, ...)
{
  coef.cov.all <- x$coef.cov.all
  varcor <- list()

  for(i in seq_along(coef.cov.all)){

    vc <- coef.cov.all[[i]]
    stddev <- sqrt(diag(vc))
    cor <- scale(vc, center=FALSE, scale=stddev) / stddev

    attr(cor, "scaled:scale") <- NULL
    attr(vc, "stddev") <- stddev
    attr(vc, "correlation") <- cor

    group.name <- deparse(x$group.call[[i]][[2L]])
    varcor[[group.name]] <- vc
    attr(varcor, "sc") <- sigma * sigma.mhglm(x)
    attr(varcor, "useSc") <- !(x$family$family %in% c("binomial", "poisson"))
    class(varcor) <- "VarCorr.mhglm"
  }
  varcor
}


coef.mhglm <- function(object, ...)
{
    stop("Use 'fixef' and 'ranef' methods for mhglm objects")
}


fitted.mhglm <- function(object, ...)
{
    predict(object, type="response")
}


weights.mhglm <- function (object, ...)
{
    res <- object$prior.weights
    if (is.null(object$na.action))
        res
    else naresid(object$na.action, res)
}


residuals.mhglm <- function(object, type = c("deviance", "pearson", "response"), ...)
{
    type <- match.arg(type)

    y <- model.response(model.frame(object), "any")
    mu <- fitted(object)
    wts <- object$prior.weights
    if (is.null(wts))
        wts <- rep(1, NROW(y))

    res <- switch(type, deviance = {
            d.res <- sqrt(pmax((object$family$dev.resids)(y, mu, wts), 0))
            ifelse(y > mu, d.res, -d.res)
        }, pearson = {
            (y - mu) * sqrt(wts)/sqrt(object$family$variance(mu))
        }, response = {
            y - mu
        })
    if (!is.null(object$na.action))
        res <- naresid(object$na.action, res)
    res
}


ranef.mhglm <- function
(object,condVar = FALSE){

  coef.mean  <- object$fit$coefficient.mean
  dispersion <- object$dispersion

  coef.cov.all <- object$coef.cov.all
  ngroupslevel <- length(coef.cov.all)

  # compute ranef
  r <- 1
  while (r <= ngroupslevel){
    object$fit <- ebayes.est.topdown(object$fit, condVar,
				     coef.mean,
				     coef.cov.all[[r]])
    r <- r+1
  }

  # print out ranef
  coefficients.eb<- as.list(rep(NULL,ngroupslevel))
  for(r in seq_len(ngroupslevel)){
    coefficients.eb[[r]] <- ebayes.est.print(object$fit,r, condVar)
  }
  class(coefficients.eb) <- 'ranef.mhglm'
  return(coefficients.eb)
}


summary.mhglm <- function(object, ...)
{
    # fixed effects
    vcov <- vcov(object)
    coefs <- cbind("Estimate" = fixef(object),
                   "Std. Error" = sqrt(diag(vcov)))
    coefs <- cbind(coefs, (cf3 <- coefs[,1]/coefs[,2]), deparse.level=0)
    colnames(coefs)[3] <- "z value"
    coefs <- cbind(coefs, "Pr(>|z|)" =
                          2*pnorm(abs(coefs[,3]), lower.tail = FALSE))

    # random effects
    varcor <- VarCorr(object)

    structure(list(call = object$call, family = object$family,
                   coefficients = coefs, dispersion = object$dispersion,
                   vcov = vcov, varcor = varcor),
              class = "summary.mhglm")
}


print.VarCorr.mhglm <- function(x, digits = max(3, getOption("digits") - 2),
                                var.print = FALSE, ...)
{
    dims <- sapply(x, ncol)
    pmax <- max(1L, max(dims))

    table <- matrix("", sum(dims + attr(x, "useSc")), pmax + 3L)

    rownames(table) <- rep("", nrow(table))
    colnames(table) <- c("Groups", "Name", "Variance", "Std.Dev", rep("", pmax - 1L))
    if (pmax > 1L) {
        colnames(table)[5L] <- "Corr"
    }

    off <- 0L
    for (i in seq_along(x)) {
        xx <- x[[i]]
        stddev <- attr(xx, "stddev")
        cor <- attr(xx, "correlation")
        p <- ncol(cor)

        tab <- matrix("", p, p + 3L)
        tab[1L,1L] <- names(x)[[i]]
        tab[,2L] <- names(stddev)
        tab[,3L] <- format(diag(xx), digits = digits, ...)
        tab[,4L] <- format(stddev, digits = digits, ...)

        if (p > 1L) {
            cor.str <- format(cor, digits = max(2L, digits - 2L), ...)
            cor.str[row(cor.str) <= col(cor.str)] <- ""
            tab[,5L:(p+3L)] <- cor.str[,-p]
        }
        table[off + seq_len(p), seq_len(p + 3L)] <- tab
        off <- off + p
    }

    if (attr(x, "useSc")) {
        sigma <- attr(x, "sc")
        table[off + 1L, 1L] <- "Residual"
        table[off + 1L, 3L] <- format(sigma^2, digits = digits, ...)
        table[off + 1L, 4L] <- format(sigma, digits = digits, ...)
    }

    if (!var.print) {
        table <- table[,-3L]
    }

    print.table(table)
}


print.summary.mhglm <- function(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

    cat("\nRandom effects:\n")
    print(x$varcor, var.print=TRUE)

    cat("\nFixed effects:\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)

    cat("\n")
}


print.mhglm <- function(x, digits = max(3L, getOption("digits") - 3L),
                        signif.stars = getOption("show.signif.stars"), ...)
{
    names.fixed <- names(x$fit$coefficient.mean)
    names.random <- lapply(x$coef.cov.all, colnames)
    #names.random <- colnames(x$coefficient.cov)

    # random effects
    random <- list()

    for(r in seq_along(x$coef.cov.all)){

      var <- pmax(0, diag(x$coef.cov.all[[r]]))
      sd <- sqrt(var)
      random.tmp <- cbind(var, sd)
      colnames(random.tmp) <- c("Variance", "Std. Dev.")
      rownames(random.tmp) <- names.random[[r]]

      group.name <- deparse(x$group.call[[r]][[2L]])
      random[[group.name]] <- random.tmp
    }

    # fixed effects
    est <- x$fit$coefficient.mean
    est.se <- sqrt(pmax(0, diag(x$fit$coefficient.mean.cov)))
    tstat <- est / est.se
    if (x$family$family %in% c("binomial", "poisson")) {
        pval <- 2 * pnorm(-abs(tstat))
    } else {
        pval <- 2 * pt(-abs(tstat), df=x$df.residual) #df here is ad-hoc
    }
    fixed <- cbind(est, est.se, tstat, pval)
    colnames(fixed) <- c("Estimate", "Std. Error", "t value", "Pr(> |t|)")
    rownames(fixed) <- names.fixed

    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

    cat("\nRandom effects:\n")
    for(r in names(random)){
      cat(r,"\n")
      print(random[[r]], digits = digits, na.print = "NA", ...)
    }

    cat("\nFixed effects:\n")
    printCoefmat(fixed, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)

    cat("\n(Dispersion parameter for ", x$family$family,
        " family taken to be ", format(x$dispersion, digits=digits), ")\n",
        sep="")

    cat("\n")
}
