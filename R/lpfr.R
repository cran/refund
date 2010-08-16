lpfr <-
function(Y, subj, covariates=NULL, funcs, kz=30, kb=30, smooth.cov=FALSE, family="gaussian", ...) {
	require(mgcv)

	kb = min(kz, kb)
	n = length(Y)
	t = seq(0, 1, length = dim(funcs)[2])
	p = ifelse(is.null(covariates), 0, dim(covariates)[2])
	N_subj=length(unique(subj))
		
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

	## set up the design and penalty matrices
	C = funcs %*% eigenDecomp$vectors[ ,1:kz ]
	J = t(eigenDecomp$vectors[,1:kz]) %*% phi
	CJ = C %*% J
	
	Z1=matrix(0, nrow=length(Y), ncol=length(unique(subj)))
	for(i in 1:length(unique(subj))){
		Z1[which(subj==unique(subj)[i]),i]=1
	}

	X = cbind(1, covariates, CJ, Z1)
	D1 = diag(c(rep(0, 1+p+2), rep(1, kb-2), rep(0, N_subj)))
	D2 = diag(c(rep(0, 1+p+kb), rep(1, N_subj)))
	
	## fit the model
	fit = gam(Y~X-1, paraPen=list(X=list(D1, D2)), family=family, method="REML", ...)

	## get the coefficient and betaHat estimates
	coefs = fit$coef
	fitted.vals <- as.matrix(X[,1:length(coefs)]) %*% coefs
	
	beta.covariates = coefs[1:(p+1)]
	betaHat <- phi %*% coefs[(2+p):(p+1+kb)]

	ranef = coefs[(p+2+kb):(length(coefs))]

	## get the covariance matrix of the estimated functional coefficient
	varBeta=fit$Vp[(2+p):(p+1+kb),(2+p):(p+1+kb)]
	varBetaHat=phi%*%varBeta%*%t(phi)
	
	## construct upper and lower bounds for betahat
	Bounds = cbind(betaHat + 1.96*(sqrt(diag(varBetaHat))), 
		betaHat - 1.96*(sqrt(diag(varBetaHat))))
	
	ret <- list(fit, fitted.vals, betaHat, beta.covariates, ranef, X, phi, psi, varBetaHat, Bounds)
	names(ret) <- c("fit", "fitted.vals", "betaHat", "beta.covariates", "ranef", "X", "phi", 
		"psi", "varBetaHat", "Bounds")
	ret
}

