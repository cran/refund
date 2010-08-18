fosr.perm <-
function(yfdobj, Z, L=NULL, Z0=NULL, L0=NULL, 
 eval.pts = seq(min(yfdobj$basis$range), max(yfdobj$basis$range),   
length.out = 201), 
 lambda=NULL, lambda0=NULL, multi.sp=FALSE, nperm, level=.05, plot=TRUE, xlabel="", title=NULL, prelim=15, ...) {
    fpobj1 = fosr.perm.fit(yfdobj, Z=Z, L=L, Z0=Z0, L0=L0, eval.pts=eval.pts, lambda=lambda, lambda0=lambda0, multi.sp=multi.sp, nperm=nperm, prelim=prelim, ...)
    fpobj2 = fosr.perm.test(fpobj1, level=level)
    if (plot) plot(fpobj2, xlabel=xlabel, title=title)
    fpobj2
}

