fpcr <- function(y, sigmat = NULL, basismat = NULL, halfpens = NULL, 
                 fdobj = NULL, eval.pts = if (!is.null(fdobj)) seq(fdobj$basis$rangeval[1],
                    fdobj$basis$rangeval[2], length.out=401) else NULL, nc, covt = NULL, 
                 mean.signal.term = FALSE, 
                 family = gaussian(), method="REML", gamma = 1, sp=NULL, min.sp=NULL) {
    if (any(is.null(sigmat),is.null(basismat),is.null(halfpens)) &
        is.null(fdobj)) 
        stop("Must specify either (sigmat, basismat, halfpens) or fdobj")
    require(mgcv)
    
    if (!is.null(fdobj)) {
    	sigmat = t(eval.fd(eval.pts, fdobj))
    	basismat = eval.basis(eval.pts, fdobj$basis)
    }

    # Design matrix
    X0 <- if (mean.signal.term) cbind(1, rowMeans(sigmat)) 
          else matrix(1, length(y), 1)    
    if (!is.null(covt)) X0 <- cbind(X0, as.matrix(covt))
    n.unpen.cols <- 1 + mean.signal.term + if (!is.null(covt)) NCOL(covt) else 0
    
    # Decorrelate the signals from the other columns of the design matrix
    sigs.decor <- if (ncol(X0) == 1) scale(sigmat, center = TRUE, scale = FALSE) 
                  else lm(sigmat ~ X0 - 1)$resid
    SB <- sigs.decor %*% basismat
    V.nc <- svd(SB)$v[ , 1 : nc]
    X <- cbind(X0, SB %*% V.nc)

    # Define list of penalties for paraPen
    S <- list()
    npen = if (!is.null(fdobj)) 1 else length(halfpens)
    for (i in 1:npen) {
        S[[i]] <- matrix(0, ncol(X), ncol(X))
        S[[i]][-(1:n.unpen.cols), -(1:n.unpen.cols)] <- 
             if (!is.null(fdobj)) crossprod(V.nc, getbasispenalty(fdobj$basis) %*% V.nc) else crossprod(halfpens[[i]] %*% V.nc) 
    }

    obje = gam(y~X-1, paraPen=list(X=S), family=family, method=method, gamma=gamma, sp=sp, min.sp=min.sp)
    obje$fhat = basismat %*% V.nc %*% obje$coef[-(1:n.unpen.cols)]
    obje$nc = nc
    obje
}