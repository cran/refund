\name{pfr}
\alias{pfr}
\title{
Penalized Functional Regression
}
\description{
Implements penalized functional regression (Goldsmith et al. 2011) for generalized linear functional models with scalar outcomes.
}
\usage{
pfr(Y, covariates = NULL, funcs, kz = 30, kb = 30, smooth.cov = FALSE, family = "gaussian", ...)
}
\arguments{
  \item{Y }{Length \code{n} vector of outcomes}
  \item{covariates }{\code{n} by \code{p} matrix of scalar covariates}
  \item{funcs }{\code{n} by \code{N} matrix of observed functional predictors}
  \item{kz }{Dimension of principal components basis for the observed functional predictors}
  \item{kb }{Dimension of the truncated power series spline basis for the coefficient function}
  \item{smooth.cov }{Logical; do you wish to smooth the covariance matrix of observed functions? Increases computation time, but results in smooth principal components}
  \item{family }{Generalized linear model family}
  \item{... }{Additional arguments passed to the \code{gam} function to fit the regression model.}
}
\value{
  \item{gam }{The result of the call to \code{gam}}
  \item{fitted.vals }{Predicted outcomes}
  \item{betaHat }{Estimated coefficient function}
  \item{beta.covariates }{Parameter estimates for scalar covariates}
  \item{X }{Design matrix used in the model fit}
  \item{phi }{Truncated power series spline basis for the coefficient function}
  \item{psi }{Principal components basis for the functional predictors}
  \item{varBetaHat }{Covariance matrix for the estimated coefficient function}
  \item{Bounds }{Bounds of a 95\% confidence interval for the estimated coefficient function}
}
\references{
Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., and Reich, D. (2011). Penalized functional regression. \emph{Journal of Computational and Graphical Statistics}, to appear.
}
\author{
Jeff Goldsmith <jgoldsmi@jhsph.edu>
}
\examples{
##### Octane Data Example #####
data(gasoline)
Y = gasoline$octane
funcs = gasoline$NIR
wavelengths = 2*450:850

# fit the model using pfr
fit = pfr(Y=Y, funcs=funcs, kz=50, kb=50)

# plot the estimated coefficient function
matplot(wavelengths, cbind(fit$betaHat, fit$Bound), type='l', lwd=c(2,1,1), 
    lty=c(1,2,2), xlab = "Wavelengths", ylab = "Coefficient Function", col=1)
}