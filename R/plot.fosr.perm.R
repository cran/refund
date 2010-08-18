plot.fosr.perm <-
function(x, level=.05, xlabel="", title=NULL,...) {
	if (is.null(level)) {
		if (is.null(x$level)) stop("Must specify level at which to test")
		else testobj = x
    }
    
    else testobj = fosr.perm.test(x, level=level)
    
	eval.pts = testobj$eval.pts
	F = testobj$F
	F.perm = testobj$F.perm
	
    matplot(eval.pts, t(rbind(F, F.perm)), type='l', col='grey', lty=1, ylab="F statistics", xlab=xlabel, main=title,...)
    abline(h=testobj$critval, col=1+1:length(testobj$level), lty=2)
    lines(eval.pts, F, col='blue')
}

