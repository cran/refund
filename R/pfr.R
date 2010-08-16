pfr <-
function(Y, covariates=NULL, funcs, kz=30, kb=30, smooth.cov=FALSE, family="gaussian", ...) {
	require(mgcv)

	kb = min(kz, kb)
	n = length(Y)
	t = seq(0, 1, length = dim(funcs)[2])
	p = ifelse(is.null(covariates), 0, dim(covariates)[2])
	
	## get the eigen decomposition of the smoothed variance matrix
	if(smooth.cov){
		G2 <- var(funcs)
		M <- length(t)
		diag(G2)= rep(NA, M)
		g2 <- as.vector(G2)
		## define a N*N knots for bivariate smoothing
		N <- 10

		## bivariate smoothing using the gamm function
		t1 <- rep(t, each=M)
		t2 <- rep(t, M)
		newdata <- data.frame(t1 = t1, t2 = t2)
		K.0  <- matrix(predict(gamm(as.vector(g2)  ~ te(t1, t2, k = N))$gam, newdata), M, M) # smooth K.0
		K.0 <- (K.0 + t(K.0)) / 2    
		
		eigenDecomp <- eigen(K.0)
	} else {
		varFuncs <- var(funcs)
		eigenDecomp <- eigen(varFuncs)
	}
	
	psi = eigenDecomp$vectors[,1:kz]

	# set the basis to be used for beta(t)
	num=kb-2
	qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
	knots <- quantile(t, qtiles)
	phi = cbind(1, t, sapply(knots, function(k) ((t - k > 0) * (t - k))))

	## set up the design and penalty matrix
	C = funcs %*% eigenDecomp$vectors[ ,1:kz ]
	J = t(eigenDecomp$vectors[,1:kz]) %*% phi
	CJ = C %*% J

	X = cbind(1, covariates, CJ)
	D = diag(c(rep(0, 1+p+2), rep(1, kb-2)))
	
	## fit the model
	gam = gam(Y~X-1, paraPen=list(X=list(D)), family=family, method="REML", ...)

	## get the coefficient and betaHat estimates
	coefs = gam$coef
	fitted.vals <- as.matrix(X[,1:length(coefs)]) %*% coefs
	
	beta.covariates = coefs[1:(p+1)]
	betaHat <- phi %*% coefs[-(1:(p+1))]

	## get the covariance matrix of the estimated functional coefficient
	varBeta=gam$Vp[-(1:(1+p)),-(1:(1+p))]
	varBetaHat=phi%*%varBeta%*%t(phi)
	
	## construct upper and lower bounds for betahat
	Bounds = cbind(betaHat + 1.96*(sqrt(diag(varBetaHat))), 
		betaHat - 1.96*(sqrt(diag(varBetaHat))))
	
	ret <- list(gam, fitted.vals, betaHat, beta.covariates, X, phi, psi, varBetaHat, Bounds)
	names(ret) <- c("gam", "fitted.vals", "betaHat", "beta.covariates", "X", "phi", 
		"psi", "varBetaHat", "Bounds")
	ret
}

