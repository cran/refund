fosr <-
function(fdobj, Z, L=NULL, eval.pts=seq(min(fdobj$basis$range),max(fdobj$basis$range), length.out=201), method='OLS', gam.method='REML', lambda=NULL, multi.sp=FALSE, max.iter=1, maxlam=NULL, cv1=FALSE, scale=FALSE) {
	require(fda)
	if (!is.fd(fdobj)) stop('First argument must be a functional data object')
	if (method==substr("OLS", 1, nchar(method))) method = "OLS"
	else if (method==substr("GLS", 1, nchar(method))) method = "GLS"
    if (!method %in% c("OLS", "GLS")) stop("Method must be either OLS or GLS")
	if (method!="OLS" & length(lambda)>1) stop("Vector-valued lambda allowed only if method = OLS")
    if (!is.null(lambda) & multi.sp) stop("Fixed lambda not implemented with multiple penalties")
    if (method=="OLS" & multi.sp) stop("OLS not implemented with multiple penalties")
    require(mgcv)
	
	nsub = ncol(fdobj$coefs)
	bss = fdobj$basis
	nbasis = bss$nbasis
	newfit = NULL

    Z = scale(Z, center=FALSE, scale=scale)	
    q = ncol(Z)

    J = getbasispenalty(bss, 0)  
    svdJ = svd(J)
    J12 = svdJ$u %*% diag(sqrt(svdJ$d)) %*% t(svdJ$u)
    
    if (multi.sp) {
    	pen = vector("list", q)
    	for (j in 1:q) {
    		one1 = matrix(0,q,q); one1[j,j] = 1
    		pen[[j]] = one1 %x% getbasispenalty(bss)
    	}	
    }
    else pen = list(diag(q) %x% getbasispenalty(bss)) 
    
    C = t(fdobj$coefs)
    constr = if (!is.null(L)) L %x% diag(nbasis) else NULL   
    cv = NULL
    
    if (method=="OLS") {
    	if (length(lambda)!=1 | cv1) {  # compute (optimal) LOFO-CV
            lofo = lofocv(C %*% J12, Z %x% J12, S1=pen[[1]], lamvec=lambda, constr=constr, maxlam=maxlam)
            cv = if (is.null(lambda)) lofo$objective else min(lofo[,2]) 
            lambda = if (is.null(lambda)) lofo$min else lofo[which.min(lofo[,2]), 1] 
        }
    }
    
    firstfit = amc(as.vector(tcrossprod(J12, C)), Z %x% J12, gam.method=gam.method, S=pen, C=constr, lambda=lambda)
    B = B.ols = t(matrix(firstfit$coef, ncol=q))
    se = NULL
    
    if (method=="GLS") {
        iter = 0                         
        B.old = 3 * B.ols; newfit = NULL
            
        if (!is.null(lambda) & max.iter>0) warning('Given lambda used for initial fit only')
        
        while (any(abs((B-B.old) / B.old) > .001) & (iter<max.iter)) {            
            iter = iter + 1
            if (max.iter>1) cat('Refit', iter, '\n')
            if (!is.null(newfit)) oldfit = newfit
            else oldfit = firstfit
	        B.old = B
	        resid = as.vector(tcrossprod(J12, C)) - (Z %x% J12) %*% oldfit$coef
            resmat = t(matrix(resid, ncol=nsub))                
            cholW = chol((nsub/(nsub-1)) * solve(cov(resmat)))             
            newfit = amc(as.vector(tcrossprod(cholW %*% J12, C)), Z %x% (cholW %*% J12), gam.method=gam.method, S=pen, C=constr)  
            B = t(matrix(newfit$coef, ncol=q))
        }
    }

    if (method=="OLS" | max.iter==0) {
        covmat = ((nsub-1)/nsub) * cov(t(matrix(as.vector(tcrossprod(J12, C)) - (Z %x% J12) %*% firstfit$coef, ncol=nsub)))
         var.b = firstfit$GinvXT %*% (diag(nsub) %x% covmat) %*% t(firstfit$GinvXT)               
     }
            
    else if (method=="GLS" & max.iter>0) var.b = newfit$Vp
        
    evb = eval.basis(eval.pts, bss)
    se.func = matrix(NA, length(eval.pts), q)
    for (j in 1:q) {
    	W = chol(var.b[(nbasis*(j-1)+1):(nbasis*j), (nbasis*(j-1)+1):(nbasis*j)])
    	se.func[ , j] = sqrt(apply(tcrossprod(W, evb), 2, crossprod))
    }
   
    est.func = eval.fd(eval.pts, fd(t(B), bss))
    fit = if (method=="GLS" & max.iter>0) newfit else firstfit
	roughness = diag(B %*% getbasispenalty(bss) %*% t(B)) 
    
    # Reverse scaling
    skale = attr(Z,"scaled:scale")
    if (!is.null(skale)) {
    	B = t(scale(t(B), center=FALSE, scale=skale))
    	est.func = scale(est.func, center=FALSE, scale=skale)
    	se.func = scale(se.func, center=FALSE, scale=skale)
        roughness = roughness / skale^2  
    }
    
    llist = list(B = B, yhat = fd(t(Z %*% B), bss), est.func = est.func, se.func = se.func, eval.pts = eval.pts, fit=fit, edf=sum(fit$gam$edf), lambda=if (length(fit$gam$sp)>0) fit$gam$sp else fit$gam$full.sp,
cv=cv, roughness = roughness)
    
    class(llist) = "fosr"
    llist    
}

