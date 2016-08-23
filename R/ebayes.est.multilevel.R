ebayes.bottom.est <- function
## This is a unit function.
## This function takes in group specific data and global parameter estimates.
## It returns group specific statistics for later MAP inference.
(coefficients, nfixed, subspace, precision, dispersion,
 coefficient.mean, coefficient.cov
 ) {
  coef <- coefficients
  coef.mu <- coefficient.mean
  coef.cov <- coefficient.cov
  r <- length(precision)

  nrandom <- nrow(coef.cov)
  fixed <- seq_len(nfixed)
  random <- nfixed + seq_len(nrandom)

  if (r == 0L) {
    coef.eb <- numeric(nrandom)
    cov.eb <- coef.cov
  } else {
    # implementation trick to avoid 1/li:
    # U (U^T Sigma U + a L^{-1})^{-1} U^T
    #   = U L^{1/2} (L^{1/2} U^T Sigma U L^{1/2} + a I)^{-1} L^{1/2} U^T
    #   = Us (Us^T Sigma Us + a I)^{-1} Us^T
    u <- subspace
    s <- sqrt(precision)
    us <- u %*% diag(s, r, r)
    u1s <- us[fixed,,drop=FALSE]
    u2s <- us[random,,drop=FALSE]

    cov.ii <- t(u2s) %*% coef.cov %*% u2s
    h <- cov.ii + diag(dispersion, r, r)

    w.inv <- pseudo.solve(h)
    w.B <- u2s %*% w.inv %*%  (t(u1s) %*% (coef[fixed] - coef.mu)
			       + t(u2s) %*% coef[random])
    w.A <- u2s %*% w.inv %*% t(u2s)

    coef.eb.B <- coef.cov %*% w.B
    coef.eb.A <- coef.cov %*% w.A

  }


  list(coef.eb.A = coef.eb.A, coef.eb.B = coef.eb.B )
}


ebayes.bottom.est.rec<- function
## This is a recursive function.
## It takes in fit object from mhglm.fit.multilevel and global parameter estimates.
## It returns a tree-structured list of statistics for later MAP inference.
(fit,coef.mean,dispersion,coef.cov.bottom
){
  ng <- length(fit)
  subgp <- seq_len(ng)[ sapply(fit,function(x) class(x) == 'mhglmfit',simplify = 'array') ]
  if(length(subgp)==0){

    nfixed <- ncol(fit$coefficients) - ncol(fit$coefficient.cov)
    coef.eb <- lapply(names(fit$subspace), function(x){ 
			ebayes.bottom.est(fit$coefficients[x,],nfixed,
					  fit$subspace[[x]], 
					  fit$precision[[x]],
					  dispersion,
					  c(coef.mean, rep(0,nfixed - length(coef.mean))),
					  coef.cov.bottom) }) 
    names(coef.eb) <- names(fit$subspace)
    return(coef.eb)

  } else {

    ret <- lapply(fit[subgp],function(x){
		    ebayes.bottom.est.rec(x,coef.mean,dispersion, coef.cov.bottom)
			})
    return(ret)
  }
}

ebayes.aggregate <- function
## This is a unit function.
## This function aggreaget group specific information within leaf nodes.
## coef.eb a list of list, where the sublist contains coef.A and coef.B
## for each subgroup.
## coef.cov.1 is the covariance estimate for parent node. 
## coef.cov.2 is the covariance estimate for leaf node. 
(coef.eb, coef.cov.1, coef.cov.2
){
  ngroups <- length(coef.eb)
  nvars <- ncol(coef.cov.1)

  A.sum <- Reduce("+",lapply(coef.eb,function(x) x$coef.eb.A))
  B.sum <- Reduce("+",lapply(coef.eb,function(x) x$coef.eb.B))

  cov.inv <- pseudo.solve(coef.cov.2 %*% pseudo.solve(coef.cov.1) + A.sum) 
  coef.eb.A <- cov.inv %*% A.sum
  coef.eb.B <- cov.inv %*% B.sum
  list(coef.eb.A = coef.eb.A, coef.eb.B = coef.eb.B)
}


ebayes.aggregate.bottomup<- function
## This is a recursive function. Going from bottom to top. 
## This function takes in a tree-structured list of statistics,
## returns the same tree-structured list, with additional statistics 
## at upper level nodes.
(coef.eb.bottom, coef.cov.all, r
){

  if(!is.null(coef.eb.bottom[[1]]$coef.eb.A)){
    ret <- ebayes.aggregate(coef.eb.bottom, 
			    coef.cov.1 = coef.cov.all[[r-1]], 
			    coef.cov.2 = coef.cov.all[[r]]) 
    coef.eb.bottom <- c(coef.eb.bottom,ret)
  } else {
    coef.eb.bottom <- lapply(coef.eb.bottom, function(x) ebayes.aggregate.bottomup(x,coef.cov.all,r))
  }
  return(coef.eb.bottom)
}

ebayes.est.topdown <- function
## This is a recursive function. Going from top to bottom.
## This function takes in a tree-structured list of statistics, 
## returns the same tree-structured list, with posterior estimates
## at lower level.
(coef.eb.bottom
){
  pavar <- c('coef.eb.A','coef.eb.B','coef.eb','coef.eb.cum')

  if(is.null(coef.eb.bottom[[1]]$coef.eb)){

    ng <- length(coef.eb.bottom)
    paidx <- seq_len(ng)[names(coef.eb.bottom) %in% pavar]
    subidx <- setdiff(seq_len(ng), paidx)

    ret <- lapply(coef.eb.bottom[subidx],function(x){ 
		    x$coef.eb <- x$coef.eb.B - x$coef.eb.A %*% coef.eb.bottom$coef.eb.cum
		    x$coef.eb.cum <- coef.eb.bottom$coef.eb.cum + x$coef.eb
		    x})
    ret <- c(ret,coef.eb.bottom[paidx])

  } else {

    ng <- length(coef.eb.bottom)
    paidx <- seq_len(ng)[names(coef.eb.bottom) %in% pavar]
    subidx <- setdiff(seq_len(ng), paidx)

    ret <- lapply(coef.eb.bottom[subidx],function(x){ x <- ebayes.est.topdown(x);x})
    ret <- c(ret,coef.eb.bottom[paidx])

  }
  return(ret)
}

ebayes.est.print<- function
## This is a recursive function.
## It takes in a tree-structured list of posterior estimates,
## put them into a matrix. 
(coef.eb.bottom,r
){
  if(r==1){
    ret <- lapply(Filter(function(x) !is.array(x), coef.eb.bottom),function(x) { c(x$coef.eb) } )
    ret2 <- Reduce("rbind",ret)

    if(length(ret)==1){
      # if ret only has one item, Reduce will return a vector instead of matrix
      # turn it into a matrix
      ret2 <- matrix(ret2, nrow = 1)
    }

    rownames(ret2) <- names(ret)
    return(ret2)
  } else {
#    ret <-  lapply(Filter(function(x) !is.array(x), coef.eb.bottom),function(x) { ebayes.est.print(x,r-1)})
    ng <- length(coef.eb.bottom)
    subgp <- seq_len(ng)[ sapply(coef.eb.bottom,function(x) !is.array(x) )]
    subgp.name <- names(coef.eb.bottom)[subgp]
    ret <- lapply(subgp, function(x) {
		    est <- ebayes.est.print(coef.eb.bottom[[x]],r-1);
		    rownames(est) <- paste0(subgp.name[x],":",rownames(est)) ; est}) 

    ret2 <- Reduce("rbind",ret)
    return(ret2)
  }
}


ebayes.est.multilevel <- function
## Compute random effects' posterior estimate.
(object){

  coef.mean  <- object$fit$coefficient.mean
  dispersion <- object$dispersion

  coef.cov.all <- object$coef.cov.all
  nlevels <- length(coef.cov.all)
  nrandom <- ncol(coef.cov.all[[1]])

  #----------
  # compute coef.A and coef.B on bottom level
  coef.cov.bottom  <- object$coef.cov.all[[nlevels]]
  coef.eb.bottom <- ebayes.bottom.est.rec(object$fit,coef.mean,dispersion,coef.cov.bottom)

  #----------
  # aggregate coef.A and coef.B from bottom to top
  r <- nlevels
  while (r >1){
    coef.eb.bottom <- ebayes.aggregate.bottomup(coef.eb.bottom, coef.cov.all,r)
    r <- r-1
  }

  #----------
  # compute ebayes.coef from top to bottom
  coef.eb.bottom$coef.eb.cum <- matrix(rep(0,nrandom),ncol=1)

  r <- nlevels
  while (r >0){
    coef.eb.bottom <- ebayes.est.topdown(coef.eb.bottom)
    r <- r-1
  }


  coefficients.eb<- as.list(rep(NULL,nlevels))
  for(r in seq_len(nlevels)){
    coefficients.eb[[r]] <- ebayes.est.print(coef.eb.bottom,r)
  }
  return(coefficients.eb)

}


