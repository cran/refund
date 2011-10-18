\name{pfr}
\alias{pfr}
\title{
Penalized Functional Regression
}
\description{
Implements penalized functional regression (Goldsmith et al., 2011) for generalized linear functional models with scalar outcomes.
}
\usage{
pfr(Y, covariates = NULL, funcs, kz = 30, kb = 30, 
  smooth.cov = FALSE, family = "gaussian", method = "REML", ...)
}
\arguments{
  \item{Y }{Vector of all outcomes over all visits}
  \item{covariates }{Matrix of scalar covariates}
  \item{funcs }{Matrix or list of matrices containing observed functional predictors as rows. NA values are allowed.}
  \item{kz }{Dimension of principal components basis for the observed functional predictors}
  \item{kb }{Dimension of the truncated power series spline basis for the coefficient function}
  \item{smooth.cov }{Logical; do you wish to smooth the covariance matrix of observed functions? Increases computation time, but results in smooth principal components}
  \item{family }{Generalized linear model family}
  \item{method }{Method for estimating the smoothing parameters; defaults to REML}
  \item{... }{Additional arguments passed to the \code{gam} function to fit the regression model.}
}
\details{
  Functional predictors are entered as a matrix or, in the case of multiple functional predictors, as a list of matrices using the \code{funcs} argument. Missing values are allowed in the functional predictors, but it is assumed that they are observed over the same grid. Functional coefficients and confidence bounds are returned as lists in the same order as provided in the \code{funcs} argument, as are principal component and spline bases.
}
\value{
  \item{fit }{The result of the call to \code{gam}}
  \item{fitted.vals }{Predicted outcomes}
  \item{betaHat }{List of estimated coefficient functions}
  \item{beta.covariates }{Parameter estimates for scalar covariates}
  \item{X }{Design matrix used in the model fit}
  \item{phi }{List of truncated power series spline bases for the coefficient functions}
  \item{psi }{List of principal components basis for the functional predictors}
  \item{varBetaHat }{List containing covariance matrices for the estimated coefficient functions}
  \item{Bounds }{List of bounds of a 95\% confidence interval for the estimated coefficient functions}
}
\references{
Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., and Reich, D. (2011).
Penalized functional regression. \emph{Journal of Computational and
Graphical Statistics}, to appear.
}
\author{
Jeff Goldsmith <jgoldsmi@jhsph.edu>
}
\examples{
##################################################################
#########               DTI Data Example                 #########
##################################################################

##################################################################
# use baseline data to regress continuous outcomes on functional 
# predictors (continuous outcomes only recorded for case == 1)
##################################################################

data(DTI)

# subset data as needed for this example
cca = DTI$cca[which(DTI$visit ==1 & DTI$case == 1),]
rcst = DTI$rcst[which(DTI$visit ==1 & DTI$case == 1),]
DTI = DTI[which(DTI$visit ==1 & DTI$case == 1),]


# note there is missingness in the functional predictors
apply(is.na(cca), 2, mean)
apply(is.na(rcst), 2, mean)


# fit two models with single functional predictors and plot the results
fit.cca = pfr(Y=DTI$pasat, funcs = cca, smooth.cov=FALSE)
fit.rcst = pfr(Y=DTI$pasat, funcs = rcst, smooth.cov=FALSE)

par(mfrow = c(1,2))
matplot(cbind(fit.cca$BetaHat[[1]], fit.cca$Bounds[[1]]),
  type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
  main = "CCA")
matplot(cbind(fit.rcst$BetaHat[[1]], fit.rcst$Bounds[[1]]),
  type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
  main = "RCST")


# fit a model with two functional predictors and plot the results
fit.cca.rcst = pfr(Y=DTI$pasat, funcs = list(cca, rcst), 
  smooth.cov=FALSE)

par(mfrow = c(1,2))
matplot(cbind(fit.cca.rcst$BetaHat[[1]], fit.cca.rcst$Bounds[[1]]),
  type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
  main = "CCA")
matplot(cbind(fit.cca.rcst$BetaHat[[2]], fit.cca.rcst$Bounds[[2]]),
  type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
  main = "RCST")


##################################################################
# use baseline data to regress binary case-status outcomes on 
# functional predictors
##################################################################

data(DTI)

# subset data as needed for this example
cca = DTI$cca[which(DTI$visit == 1),]
rcst = DTI$rcst[which(DTI$visit == 1),]
DTI = DTI[which(DTI$visit == 1),]


# fit two models with single functional predictors and plot the results
fit.cca = pfr(Y=DTI$case, funcs = cca, smooth.cov=FALSE, 
  family = "binomial")
fit.rcst = pfr(Y=DTI$case, funcs = rcst, smooth.cov=FALSE, 
  family = "binomial")

par(mfrow = c(1,2))
matplot(cbind(fit.cca$BetaHat[[1]], fit.cca$Bounds[[1]]),
  type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
  main = "CCA")
matplot(cbind(fit.rcst$BetaHat[[1]], fit.rcst$Bounds[[1]]),
  type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
  main = "RCST")


##################################################################
#########              Octane Data Example               #########
##################################################################

data(gasoline)
Y = gasoline$octane
#detach("package:refund")
funcs = gasoline$NIR
wavelengths = as.matrix(2*450:850)

# fit the model using pfr
fit = pfr(Y=Y, funcs=funcs, kz=50, kb=50)

# plot the estimated coefficient function
matplot(wavelengths, cbind(fit$BetaHat[[1]], fit$Bounds[[1]]), 
  type='l', lwd=c(2,1,1), lty=c(1,2,2), xlab = "Wavelengths", 
  ylab = "Coefficient Function", col=1)

}