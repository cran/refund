fpcr <- function(y, sigmat = NULL, basismat = NULL, penmats = NULL, 
                 fdobj = NULL, argvals = NULL, nc = NULL, covt = NULL, 
                 mean.signal.term = FALSE, nbasis = NULL, spline.order = NULL, pen.order = NULL,
                 family = gaussian(), method="REML", sp=NULL, ...) {
    if (is.null(sigmat) & is.null(fdobj)) 
        stop("Must specify either sigmat or fdobj")
    if (!is.null(sigmat) & !is.null(fdobj)) 
        stop("Must specify either sigmat or fdobj, but not both")
    require(mgcv)
    
    if (!is.null(fdobj)) {
    	basis = fdobj$basis
    	if (!is.null(nbasis)) if (nbasis!=basis$nbasis) 
    	    warning(paste("'nbasis = ", nbasis, "' overridden since the supplied 'fdobj' has basis dimension ", basis$nbasis, sep=""))
    	if (!is.null(spline.order)) if (spline.order!=norder(basis))     	    warning(paste("'spline.order = ", spline.order, "' overridden since the supplied 'fdobj' has a basis of order ", norder(basis), sep=""))
        if (is.null(argvals)) argvals = seq(basis$rangeval[1],
                    basis$rangeval[2], length.out=401) 
    	sigmat = t(eval.fd(argvals, fdobj))
    }
    
    else {
    	if (is.null(argvals)) argvals = seq(0,1,,ncol(sigmat))
    	if (is.null(nbasis)) nbasis = 40
    	if (is.null(spline.order)) spline.order = 4    	
    	basis = create.bspline.basis(rangeval = c(min(argvals), max(argvals)), nbasis = nbasis, norder = spline.order)
    }

    if (is.null(basismat)) basismat = eval.basis(argvals, basis)
    if (is.null(pen.order)) pen.order = 2
    if (is.null(penmats)) penmats = list(getbasispenalty(basis, pen.order))

    # Design matrix
    X0 <- if (mean.signal.term) cbind(1, rowMeans(sigmat)) 
          else matrix(1, length(y), 1)    
    if (!is.null(covt)) X0 <- cbind(X0, as.matrix(covt))
    n.unpen.cols <- 1 + mean.signal.term + if (!is.null(covt)) NCOL(covt) else 0
    
    # Decorrelate the signals from the other columns of the design matrix
    sigs.decor <- if (ncol(X0) == 1) scale(sigmat, center = TRUE, scale = FALSE) 
                  else lm(sigmat ~ X0 - 1)$resid
    SB <- sigs.decor %*% basismat
    svdSB = svd(SB)
    if (is.null(nc)) nc = min(which(cumsum(svdSB$d) > .99 * sum(svdSB$d)))
    V.nc <- svdSB$v[ , 1 : nc]
    X <- cbind(X0, SB %*% V.nc)

    # Define list of penalties for paraPen
    S <- list()
    npen = length(penmats)
    for (i in 1:npen) {
        S[[i]] <- matrix(0, ncol(X), ncol(X))
        S[[i]][-(1:n.unpen.cols), -(1:n.unpen.cols)] <- 
             crossprod(V.nc, penmats[[i]] %*% V.nc) 
    }

    obje = gam(y~X-1, paraPen=list(X=S), family=family, method=method, sp=sp, ...)
    BV = basismat %*% V.nc
    obje$fhat = BV %*% obje$coef[-(1:n.unpen.cols)]
    obje$se = sqrt(rowSums((BV %*% obje$Vp[-(1:n.unpen.cols),-(1:n.unpen.cols)]) * BV))
    obje$nc = nc
    obje$argvals = argvals
    class(obje) = "fpcr"
    obje
}