fosr.perm.fit <-
function(yfdobj, Z, L=NULL, Z0=NULL, L0=NULL, eval.pts = seq(min(yfdobj$basis$range), max(yfdobj$basis$range), length.out = 201), 
lambda=NULL, lambda0=NULL, multi.sp=FALSE, nperm, prelim, ...) {
	if (length(lambda)>1) stop("'lambda' must be a scalar or NULL")
	if (is.null(Z0) & !is.null(L0)) stop("If 'L0' is given, 'Z0' must be given too")
	df1 = ncol(Z) - 1 
	if (!is.null(L)) df1 = df1 - nrow(L)
	df2 = nrow(Z) - df1 - 1
	if (!is.null(Z0)) df1 = df1 - ncol(Z0) + 1
	if (!is.null(L0)) df1 = df1 + nrow(L0) 

	realfit = fosr(yfdobj, Z, L, eval.pts=eval.pts, lambda=lambda, multi.sp=multi.sp, ...)
	lambda.real = realfit$lambda
	lambda0.real = NULL
	if (!is.null(Z0)) {
		realfit0 = fosr(yfdobj, Z0, L0, eval.pts=eval.pts, lambda=lambda0, multi.sp=multi.sp, ...)
	    lambda0.real = realfit0$lambda
    }
    
    yvec = eval.fd(eval.pts, yfdobj)
    yhatvec = eval.fd(eval.pts, realfit$yhat)
    if (is.null(Z0)) F = (apply(yhatvec, 1, var) / df1) / (apply((yvec - yhatvec)^2, 1, mean) / df2)
    else {
    	yhatvec0 = eval.fd(eval.pts, realfit0$yhat)
    	F = ((apply((yvec - yhatvec0)^2, 1, mean) - apply((yvec - yhatvec)^2, 1, mean)) / df1) / (apply((yvec - yhatvec)^2, 1, mean) / df2)
    }

    lambda.prelim = lambda0.prelim = NULL
    if (prelim>0) {
    	if (length(lambda)==1) stop("If 'lambda' is specified, 'prelim' should be set to 0")
    	lambda.prelim = matrix(NA, prelim, ncol(Z)^multi.sp)
    	if (!is.null(Z0)) lambda0.prelim = matrix(NA, prelim, ncol(Z0)^multi.sp)

    	begin.prelim = proc.time()
    	for (ee in 1:prelim) {
    		yfdobj.perm = yfdobj
    		yfdobj.perm$coefs = yfdobj$coefs[ , sample(ncol(yfdobj$coef))]
    	    lambda.prelim[ee, ] = fosr(yfdobj.perm, Z, L, eval.pts=eval.pts, multi.sp=multi.sp, ...)$lambda 
    	    if (!is.null(Z0)) lambda0.prelim[ee, ] = fosr(yfdobj.perm, Z0, L0, eval.pts=eval.pts, multi.sp=multi.sp, ...)$lambda
    	    if (ee==1) {
    	    	elapsed.time <- max(proc.time() - begin.prelim, na.rm = TRUE)
               if (elapsed.time > 10/prelim) cat("\n***** Estimated computing time for preliminary permutations:", round(prelim * elapsed.time), "seconds *****\n")
            }
    	}
    	lambda = apply(lambda.prelim, 2, median)
    	if (!is.null(Z0)) lambda0 = apply(lambda0.prelim, 2, median)
        cat("***** Computing time for preliminary permutations:", max(proc.time() - begin.prelim, na.rm = TRUE), "seconds *****\n\n")
    }
    
    yfdobj.perm = yfdobj
    F.perm = matrix(NA, nperm, length(eval.pts))
    lambda.perm = matrix(NA, nperm, ncol(Z)^multi.sp)
    lambda0.perm = if (!is.null(Z0)) matrix(NA, nperm, ncol(Z0)^multi.sp) else NULL

    begin.perm = proc.time()    
	for (i in 1:nperm) {
		if (i/20==floor(i/20)) cat('Permutation', i, '\n')
		yfdobj.perm$coefs = yfdobj$coefs[ , sample(ncol(yfdobj$coef))]
	    fit.perm = fosr(yfdobj.perm, Z, L, eval.pts=eval.pts, lambda=lambda, multi.sp=multi.sp, ...)
	    if (!is.null(Z0)) fit.perm0 = fosr(yfdobj.perm, Z0, L0, eval.pts=eval.pts, lambda=lambda0, multi.sp=multi.sp, ...)
	    lambda.perm[i, ] = fit.perm$lambda
	    if (!is.null(Z0)) lambda0.perm[i, ] = fit.perm0$lambda
        yvec.perm = eval.fd(eval.pts, yfdobj.perm)
        yhatvec.perm = eval.fd(eval.pts, fit.perm$yhat)
        if (is.null(Z0)) F.perm[i, ] = (apply(yhatvec.perm, 1, var) / df1) / (apply((yvec.perm - yhatvec.perm)^2, 1, mean) / df2)
        else {
        	yhatvec.perm0 = eval.fd(eval.pts, fit.perm0$yhat)
        	F.perm[i, ] = ((apply((yvec.perm - yhatvec.perm0)^2, 1, mean) - apply((yvec.perm - yhatvec.perm)^2, 1, mean)) / df1) / (apply((yvec.perm - yhatvec.perm)^2, 1, mean) / df2)
        }
    	if (i==1) {
    	    elapsed.time <- max(proc.time() - begin.perm, na.rm = TRUE)
            if (elapsed.time > 10/nperm) cat("\n***** Estimated computing time for permuted-data models:", round(nperm * elapsed.time), "seconds *****\n")
        }
    }
    cat("***** Computing time for permuted-data models:", max(proc.time() - begin.perm, na.rm = TRUE), "seconds *****\n\n")
    ll = list(F=F, F.perm=F.perm, eval.pts=eval.pts, lambda.real=lambda.real, lambda.prelim=lambda.prelim, lambda.perm=lambda.perm, lambda0.real=lambda0.real, lambda0.prelim=lambda0.prelim, lambda0.perm=lambda0.perm)
    class(ll) = "fosr.perm"
    ll
}

