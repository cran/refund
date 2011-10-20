print.summary.pffr <-
function(x, digits = max(3, getOption("digits") - 3), 
        signif.stars = getOption("show.signif.stars"), ...){
    print(x$family)
    cat("Formula:\n")
    print(x$formula)
    if (length(x$p.coeff)>0)
    { cat("\nConstant coefficients:\n")
        printCoefmat(x$p.table, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    }
    cat("\n")
    if(x$m>0)
    { cat("Smooth terms & functional coefficients:\n")
        printCoefmat(x$s.table, digits = digits, signif.stars = signif.stars, has.Pvalue = TRUE, na.print = "NA",cs.ind=1, ...)
    }
    cat("\nR-sq.(adj) = ",formatC(x$r.sq,digits=3,width=5))
    if (length(x$dev.expl)>0) cat("   Deviance explained = ",formatC(x$dev.expl*100,digits=3,width=4),"%\n",sep="")
    
    if (!is.null(x$method)&&!(x$method%in%c("PQL","lme.ML","lme.REML")))  
        cat(x$method," score = ",formatC(x$sp.criterion,digits=5),sep="")
    
    cat("  Scale est. = ",formatC(x$scale,digits=5,width=8,flag="-"),"  n = ",x$n,"\n",sep="")
    invisible(x)
}

