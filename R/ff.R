ff <-
function(X,
		yind,
		xind=seq(0, 1, l=ncol(X)),
		basistype= c("te", "t2", "s"),
		integration=c("simpson", "trapezoidal"), 
		L=NULL,
		limits=NULL, 
		splinepars=list(bs="ps", m=c(2, 1))
){
	n <- nrow(X)
	nxgrid <- ncol(X)
	stopifnot(all(!is.na(X)))
	
	
	if(!missing(yind))
	if(is.null(dim(yind))){
		yind <- t(t(yind))
	} 
	nygrid <- nrow(yind)
	
	if(is.null(dim(xind))){
		xind <- t(xind)
		stopifnot(ncol(xind) == nxgrid)
		if(nrow(xind)== 1){
			xind <- matrix(as.vector(xind), nrow=n, ncol=nxgrid, byrow=TRUE)
		} 
		stopifnot(nrow(xind) == n)  
	}	
	stopifnot(all.equal(order(xind[1,]), 1:nxgrid), all.equal(order(yind), 1:nygrid))
	
	basistype <- match.arg(basistype)
	integration <- match.arg(integration)
	
	xind.sc <- xind - min(xind)
	xind.sc <- xind.sc/max(xind.sc)
	diffXind <- t(round(apply(xind.sc, 1, diff), 3))
	if(is.null(L) & any(apply(diffXind, 1, function(x) length(unique(x))) != 1) && # gridpoints for any  X_i(s) not equidistant?
			integration=="simpson"){
		warning("Non-equidistant grid detected for ", deparse(substitute(X)), ".\n Changing to trapezoidal rule for integration.")
		integration <- "trapezoidal"
		
	}
	
	if(!is.null(L)){
		stopifnot(nrow(L) == n, ncol(L) == nxgrid)
	} else {
		if(is.null(limits)){
			L <- switch(integration,
					"simpson" = {
						((xind[,nxgrid]-xind[,1])/nxgrid)/3 * 
								matrix(c(1, rep(c(4, 2), length=nxgrid-2), 1), nrow=n, ncol=nxgrid, byrow=TRUE)
					}, 
					"trapezoidal" = {
						diffs <- t(apply(xind, 1, diff))
						.5 * cbind(diffs[,1], t(apply(diffs, 1, filter, filter=c(1,1)))[,-(nxgrid-1)], diffs[,(nxgrid-1)])
					})
		} else {
			stop("<limits> not yet implemented")
		}
	}
	LX <- L*X
	xindname <- paste(deparse(substitute(X)), ".smat", sep="")
	yindname <- paste(deparse(substitute(X)), ".tmat", sep="")
	LXname <- paste("L.", deparse(substitute(X)), sep="")
	
	data <- list(
			xind[rep(1:n, times=nygrid), ], #stack xind nygrid-times
			matrix(rep(yind, times=n), nrow=n*nygrid, ncol=nxgrid), #repeat each entry of yind n times s.t. rows are constant  
			LX[rep(1:n, each=nygrid),])#stack LX nygrid-times
	names(data)  <- c(xindname, yindname, LXname)
	
	splinefun <- as.symbol(basistype) # if(basistype=="te") quote(te) else quote(s)
	frmls <- formals(getFromNamespace(deparse(splinefun), ns="mgcv"))
	frmls <- modifyList(frmls[names(frmls) %in% names(splinepars)], splinepars)
	call <- as.call(c(
					list(splinefun,
							x = as.symbol(substitute(yindname)), 
							z = as.symbol(substitute(xindname)),
							by =as.symbol(substitute(LXname))),
					frmls))
	
	return(list(call=call, data = data, yind=yind, xind=xind[1, ], L=L,
                    xindname=xindname, yindname=yindname, LXname=LXname))
}#end ff()

