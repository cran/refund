pffr <-
function(
		formula,
		yind,
		fitter = NA, 
        method="REML",
		bsy.default = list(bs="ps", m=c(2, 1)),# only bs, k, m are propagated...
		...
){
	call <- match.call()
	
	tf <- terms.formula(formula, specials=c("s", "te", "t2", "ff", "c", "sff"))
	trmstrings <- attr(tf, "term.labels")
	terms <- sapply(trmstrings, function(trm) as.call(parse(text=trm))[[1]], simplify=FALSE) 
	frmlenv <- environment(formula)
	
	
	where.c <- attr(tf, "specials")$c - 1    # indices of scalar offset terms
	where.ff <- attr(tf, "specials")$ff - 1  # function-on-function terms
    where.sff <- attr(tf, "specials")$sff - 1  #smooth function-on-function terms
	where.s <- attr(tf, "specials")$s - 1    # smooth terms 
	where.te <- attr(tf, "specials")$te - 1  # tensor product terms
	where.t2 <- attr(tf, "specials")$t2 - 1  # type 2 tensor product terms
	where.par <- which(!(1:length(trmstrings) %in%
                        c(where.c, where.ff, where.sff, where.s, where.te, where.t2))) # indices of linear/factor terms with varying coefficients over yind.
	
	responsename <- attr(tf,"variables")[2][[1]]
	
	newfrml <- paste(responsename, "~", sep="")
	newfrmlenv <- new.env()
	evalenv <- if("data" %in% names(call)) eval(call$data) else NULL
	
	nobs <- nrow(eval(responsename,  env=evalenv, enclos=frmlenv))
	nyindex <- ncol(eval(responsename,  env=evalenv, enclos=frmlenv))
	
	
	if(missing(fitter)){
		fitter <- ifelse(nobs > 1e5, "bam", "gam")
	} 
	fitter <- as.symbol(fitter)
	if(as.character(fitter)=="bam" && !("chunk.size" %in% names(call))){
		call$chunk.size <- max(nobs/5, 10000) 
	}
	if(as.character(fitter)=="gamm4") stopifnot(length(where.te)<1)
    
	
	if(missing(yind)){
		if(length(c(where.ff, where.sff))){
            if(length(where.ff)){
                ffcall <- expand.call(ff, as.call(terms[where.ff][1])[[1]])  
            }  else ffcall <- expand.call(sff, as.call(terms[where.sff][1])[[1]]) 
			if(!is.null(ffcall$yind)){
				yind <- eval(ffcall$yind, env=evalenv, enclos=frmlenv)
				yindname <- deparse(ffcall$yind)
			} else {
				yind <- 1:nyindex
				yindname <- "yindex"	
			}			
		} else {
			yind <- 1:nyindex
			yindname <- "yindex"	
		}
	} else {
		stopifnot(is.vector(yind), is.numeric(yind), 
				length(yind) == nyindex)
		yindname <- deparse(substitute(yind))
	}
	if(length(yindname)>1) yindname <- "yindex"	
	stopifnot(all.equal(order(yind), 1:nyindex))
	
	
	yindvec <- rep(yind, times = nobs)
	yindvecname <- as.symbol(paste(yindname,".vec",sep=""))
	assign(x=deparse(yindvecname), value=yindvec, envir=newfrmlenv)
	
	assign(x=deparse(responsename), value=as.vector(t(eval(responsename, env=evalenv, enclos=frmlenv))), 
			envir=newfrmlenv)
	
	newtrmstrings <- attr(tf, "term.labels")
	
    if(attr(tf, "intercept")){ 
        arglist <- c(name="s", x = as.symbol(yindvecname), bsy.default)
        assign(x= "intcall", value= do.call("call", arglist, envir=newfrmlenv), envir=newfrmlenv)
        newfrmlenv$intcall$x <- as.symbol(yindvecname)  
        
        intstring <- deparse(newfrmlenv$intcall)
        rm(intcall, envir=newfrmlenv)
        
        newfrml <- paste(newfrml, intstring, sep=" ")
        addFint <- TRUE
        names(intstring) <- paste("Intercept(",yindname,")",sep="") 
	} else{
        newfrml <-paste(newfrml, "0", sep="")
        addFint <- FALSE
    } 
	
	if(length(where.c)){ 
		newtrmstrings[where.c] <- sapply(trmstrings[where.c], function(x){
					sub("\\)$", "", sub("^c\\(", "", x)) #c(BLA) --> BLA
				})
	}
    
	if(length(c(where.ff, where.sff))){ 
		ffterms <- lapply(terms[c(where.ff, where.sff)], function(x){
					eval(x, env=evalenv, enclos=frmlenv)
				})
		newtrmstrings[c(where.ff, where.sff)] <- sapply(ffterms, function(x) safeDeparse(x$call))
		lapply(ffterms, function(x){
					lapply(names(x$data), function(nm){
								assign(x=nm, value=x$data[[nm]], envir=newfrmlenv)
								invisible(NULL)
							})
					invisible(NULL)
				})  
        ffterms <- lapply(ffterms, function(x) x[names(x)!="data"])
	} else ffterms <- NULL
    
    
    
	if(length(c(where.s, where.te, where.t2))){ 
		newtrmstrings[c(where.s, where.te, where.t2)] <- 
				sapply(terms[c(where.s, where.te, where.t2)], function(x){
						
							xnew <- x
                            if(deparse(x[[1]]) == "te" && as.character(fitter) == "gamm4") xnew[[1]] <- quote(t2)
							if(deparse(x[[1]]) == "s"){
                                xnew[[1]] <- if(as.character(fitter) != "gamm4") {
                                            quote(te) 
                                        } else quote(t2)
								xnew$d <- if(!is.null(names(xnew))){
                                            c(length(all.vars(xnew[names(xnew)!="xt"])), 1)
                                        } else c(length(all.vars(xnew)), 1)  
							} else {
								if("d" %in% names(x)){ #either expand given d...
									xnew$d <- c(eval(x$d), 1)
								} else {#.. or default to univariate marginal bases
									xnew$d <- rep(1, length(all.vars(x))+1)
								}
							}	
							xnew[[length(xnew)+1]] <- yindvecname
							
							xnew$bs <- if("bs" %in% names(x)){
										if("bs" %in% names(bsy.default)){
											c(eval(x$bs), bsy.default$bs)
										} else {
											c(xnew$bs, "tp")
										}
									} else {
										if("bs" %in% names(bsy.default)){
											c(rep("tp", length(xnew$d)-1), bsy.default$bs)
										} else {
											rep("tp", length(all.vars(x))+1)
										}	
									}
							xnew$m <- if("m" %in% names(x)){
										if("m" %in% names(bsy.default)){
											warning("overriding bsy.default for m in ", deparse(x))
											x$m
										} 
									} else {
										if("m" %in% names(bsy.default)){
											bsy.default$m
										} else {
											NA
										}	
									}
							xnew$k <- if("k" %in% names(x)){
										if("k" %in% names(bsy.default)){
											c(xnew$k, bsy.default$k)
										} 
									} else {
										if("k" %in% names(bsy.default)){
											c(5^length(all.vars(x)), bsy.default$k)
										} else {
											NA
										}	
									}
							safeDeparse(xnew)
						})
	}
    
	if(length(where.par)){ 
		newtrmstrings[where.par] <- sapply(terms[where.par], function(x){
					xnew <- bsy.default
					xnew <- as.call(c(quote(s), yindvecname, by=x, xnew))
					safeDeparse(xnew)
				})
		
	}
    where.notff <- c(where.c, where.par, where.s, where.te, where.t2)
	if(length(where.notff)){
        if("data" %in% names(call)) frmlenv <- list2env(eval(call$data), frmlenv) 
		lapply(terms[where.notff], function(x){
                if(any(sapply(terms[where.c], function(s) s[[1]]==x))){
                    x <- formula(paste("~", gsub("\\)$", "",
                                    gsub("^c\\(", "", deparse(x)))))[[2]]
                } 
                nms <- if(!is.null(names(x))){
                        all.vars(x[names(x) != "xt"]) 
                    }  else all.vars(x)
                 
                
				sapply(nms, function(nm){
							stopifnot(length(get(nm, envir=frmlenv)) == nobs)
							assign(x=nm, 
									value=rep(get(nm, envir=frmlenv), each=nyindex),
									envir=newfrmlenv)	
							invisible(NULL)
						})
				invisible(NULL)
			})
	}
	newfrml <- formula(paste(c(newfrml, newtrmstrings), collapse="+"))
	environment(newfrml) <- newfrmlenv
   
	pffrdata <- list2df(as.list(newfrmlenv))
	
	newcall <- call	
	newcall$formula <- newfrml
	newcall$data <- quote(pffrdata)
	newcall[[1]] <- fitter
	newcall$yind <- NULL
    
    m <- eval(newcall)
    m.smooth <- if(as.character(fitter) %in% c("gamm4","gamm")){
        m$gam$smooth
    } else m$smooth
    
    trmmap <- newtrmstrings
    names(trmmap) <- names(terms)
    if(addFint) trmmap <- c(trmmap, intstring)
    
    labelmap <- as.list(trmmap)
    lbls <- sapply(m.smooth, function(x) x$label)
    if(length(where.par)){
        for(w in where.par) 
            labelmap[[w]] <- {
                     where <- sapply(m.smooth, function(x) x$by) == names(labelmap)[w]
                     sapply(m.smooth[where], function(x) x$label)
                }
        labelmap[-where.par] <- lbls[pmatch(
                sapply(labelmap[-where.par], function(x){
                            tmp <- eval(parse(text=x))
                            if(is.list(tmp)){
                                return(tmp$label)  
                            } else {
                                return(x)   
                            } 
                        }), lbls)]
    } else{
        labelmap[1:length(labelmap)] <-  lbls[pmatch(
                        sapply(labelmap, function(x){
                                  tmp <- eval(parse(text=x))
                                  if(is.list(tmp)){
                                      return(tmp$label)  
                                  } else {
                                      return(x)   
                                  }
                              }), lbls)]
    } 
    if(any(nalbls <- sapply(labelmap, function(x) any(is.na(x))))){
        labelmap[nalbls] <- trmmap[nalbls]
    }
    
    names(m.smooth) <- lbls
    if(as.character(fitter) %in% c("gamm4","gamm")){
        m$gam$smooth <- m.smooth 
    } else{
        m$smooth  <- m.smooth
    } 
    
    ret <-  list(formula=formula, 
            termmap=trmmap, 
            labelmap=labelmap, 
            responsename = responsename,
            nobs=nobs,
            nyindex=nyindex,
            yindname = yindname,
            yind=yind,
            where=list(
                    where.c=where.c,
                    where.ff=where.ff,
                    where.sff=where.sff,
                    where.s=where.s,
                    where.te= where.te,
                    where.t2=where.t2,
                    where.par=where.par
            ),
            ff=ffterms)
    
    if(as.character(fitter) %in% c("gamm4","gamm")){
        m$gam$pffr <- ret
        class(m$gam) <- c("pffr", class(m$gam))
    } else {
        m$pffr <- ret
        class(m) <- c("pffr", class(m))
    }
    return(m)
}# end pffr()	

