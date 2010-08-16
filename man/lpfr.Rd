\name{lpfr}
\alias{lpfr}
\title{
Longitudinal Penalized Functional Regression
}
\description{
Implements longitudinal penalized functional regression (Goldsmith et al. 2010) for generalized linear functional models with scalar outcomes and subject-specific random intercepts.
}
\usage{
lpfr(Y, subj, covariates = NULL, funcs, kz = 30, kb = 30, smooth.cov = FALSE, family = "gaussian", ...)
}
\arguments{
  \item{Y }{Vector of all outcomes over all visits}
  \item{subj }{Vector containing the subject number for each observation}
  \item{covariates }{Matrix of scalar covariates}
  \item{funcs }{Matrix of observed functional predictors}
  \item{kz }{Dimension of principal components basis for the observed functional predictors}
  \item{kb }{Dimension of the truncated power series spline basis for the coefficient function}
  \item{smooth.cov }{Logical; do you wish to smooth the covariance matrix of observed functions? Increases computation time, but results in smooth principal components}
  \item{family }{Generalized linear model family}
  \item{... }{Additional arguments passed to the \code{gam} function to fit the regression model.}
}
\value{
  \item{fit }{The result of the call to \code{gam}}
  \item{fitted.vals }{Predicted outcomes}
  \item{betaHat }{Estimated coefficient function}
  \item{beta.covariates }{Parameter estimates for scalar covariates}
  \item{ranef }{Vector of subject-specific random intercepts}
  \item{X }{Design matrix used in the model fit}
  \item{phi }{Truncated power series spline basis for the coefficient function}
  \item{psi }{Principal components basis for the functional predictors}
  \item{varBetaHat }{Covariance matrix for the estimated coefficient function}
  \item{Bounds }{Bounds of a 95\% confidence interval for the estimated coefficient function}
}
\references{
Goldsmith, J., Crainiceanu, C., Caffo, B., and Reich, D. (2010). Longitudinal penalized functional regression.  Available at \url{http://www.bepress.com/jhubiostat/paper216}
}
\author{
Jeff Goldsmith <jgoldsmi@jhsph.edu>
}
\examples{

set.seed(100)

## set a few constants
t=seq(0,10, length=101)		# grid on which functions are observed
N_obs=length(t)				# dimension of the above grid
I=100  						# number of subjects
J = 3						# number of visits
subj=rep(1:I, each=J)		# vector of subject numbers
VarX = .5					# measurement error variance
VarY = 10					# variance on the outcome
VarRanEf = 50				# random effect variance

## define the true beta function
trueBeta = 1*sin(t*pi/5)


## generate true functions
true.funcs <- matrix(0, nrow=I*J, ncol=N_obs)
for(i2 in 1:(I*J)){
	true.funcs[i2,]=true.funcs[i2,]+runif(1, 0, 5)
	true.funcs[i2,]=true.funcs[i2,]+rnorm(1, 1, 0.2)*t
 	for(j2 in 1:10){
		e=rnorm(2, 0, 1/j2^(2))
		true.funcs[i2,]=true.funcs[i2,]+e[1]*sin((2*pi/10)*t*j2)
		true.funcs[i2,]=true.funcs[i2,]+e[2]*cos((2*pi/10)*t*j2) 
	}
}

## add measurement error to the true functions
funcs = true.funcs + rnorm(I*J*N_obs, sd=sqrt(VarX))

## generate random effects: subject specific intercepts
ranef=rep(rnorm(I, 0, sd=sqrt(VarRanEf)), each=J)

## generate outcomes from true gamma function and functions measured without error
Y <- sapply(1:(I*J), function(u) sum(true.funcs[u,]*trueBeta)+ranef[u])+rnorm(I*J, 0, sqrt(VarY))


## fit the model using the lpfr function
fit = lpfr(Y=Y, funcs=funcs, subj=subj)

## plot the true and estimated coefficient function
par(mfrow=c(1,2))
plot(trueBeta, type='l', lwd=2)
points(fit$betaHat, type='l', lwd=2, col="red")

points(fit$Bounds[,1], type = 'l', lty =2, col = "red")
points(fit$Bounds[,2], type = 'l', lty =2, col = "red")

plot(ranef[J*(1:I)], fit$ranef)
}



