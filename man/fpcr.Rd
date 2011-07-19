\name{fpcr}
\Rdversion{1.1}
\alias{fpcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional principal component regression}
\description{
  Implements functional principal component regression (Reiss and Ogden, 2007, 2010) for generalized linear models with scalar responses and functional predictors.}
\usage{
fpcr(y, sigmat = NULL, basismat = NULL, penmats = NULL, 
     fdobj = NULL, argvals = NULL, nc = NULL, covt = NULL, 
     mean.signal.term = FALSE, nbasis = NULL, 
     spline.order = NULL, pen.order = NULL,
     family = gaussian(), method="REML", sp=NULL, \dots) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
scalar outcome vector.}
  \item{sigmat}{
matrix of signals: n x N, where n is the length of y and N is the
 number of sites at which each signal is defined.}
  \item{basismat}{
N x K matrix of values of K basis functions at the N sites}
  \item{penmats}{
a list, each of whose components is a K x K matrix defining a penalty on the basis coefficients.}
\item{fdobj}{functional data object (class \code{fd}) giving the functional predictors.}
\item{argvals}{points at which the signals (functional predictors) and the coefficient function are evaluated.  By default, if the functional predictors are given by the n x N  matrix \code{sigmat}, \code{argvals} is set to N equally spaced points from 0 to 1; if they are given by \code{fdobj}, \code{argvals} is set to 401 equally spaced points spanning the domain of the given functions.}
  \item{nc}{
number of principal components. By default, the smallest number accounting for at least 99\% of the variance in the predictors.}
  \item{covt}{
covariates: an n-row matrix, or a vector of length n.}
  \item{mean.signal.term}{
logical: should the mean of each subject's signal be 
included as a covariate?}
  \item{nbasis}{number of basis functions.  Ignored if \code{fdobj} is supplied.  If \code{fdobj} is \emph{not} supplied, this defaults to 40, i.e., 40 B-spline basis functions are used.}
  \item{spline.order}{order of B-splines used, if {fdobj} is not supplied; defaults to 4, i.e., cubic B-splines.}
  \item{pen.order}{order of derivative penalty applied when estimating the coefficient function; defaults to 2.}
  \item{family}{
generalized linear model family.}
  \item{method}{
smoothing parameter selection method, passed to function \code{\link{gam}} in package \code{mgcv}; see the \code{\link{gam}} documentation for details.}
  \item{sp}{
a fixed smoothing parameter; if \code{NULL}, an optimal value is chosen (see \code{method}).}
  \item{\dots}{
other arguments passed to function \code{\link{gam}} in package \code{mgcv}.}
}
\details{
The functional predictors are given either in functional data object form, using argument \code{fdobj} (see the \code{fda} package of Ramsay, Hooker and Graves, 2009, and Method 1 in the example below), or explicitly, using \code{sigmat} (see Method 2 in the example).  In the latter case, arguments \code{basismat} and \code{penmats} can also be used to specify the basis and/or penalty matrices (see Method 3).}
\value{
An object of class \code{gam} (see \code{\link{gamObject}} in the \code{mgcv} package documentation), with  
additional components \code{nc} and \code{argvals} (same as the arguments of the same names) as well as \code{fhat} (coefficient 
 function estimate) and \code{se} (pointwise Bayesian standard error).}
\references{
Ramsay, J. O., Hooker, G., and Graves, S. (2009).
\emph{Functional Data Analysis with R and MATLAB}.
New York: Springer.

Reiss, P. T., and Ogden, R. T. (2007).  Functional principal component regression and functional partial least squares.  \emph{Journal of the American Statistical Association}, 102, 984--996.

Reiss, P. T., and Ogden, R. T. (2010).  Functional generalized linear models with images as predictors.  \emph{Biometrics}, 66, 61--69.

Wood, S. N. (2006). \emph{Generalized Additive Models: An Introduction with R}. Boca Raton, FL: Chapman &
  Hall.}

\author{
Philip Reiss \email{phil.reiss@nyumc.org}}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{##### Octane data example #####
data(gasoline)

# Create the requisite functional data objects
bbasis = create.bspline.basis(c(900, 1700), 40)
wavelengths = 2*450:850
nir = matrix(NA, 401, 60)
for (i in 1:60) nir[ , i] = gasoline$NIR[i, ]
# Why not just take transpose of gasoline$NIR above?  
# Because for some reason it leads to an error in the following statement
gas.fd = smooth.basisPar(wavelengths, nir, bbasis)$fd

# Method 1: Call fpcr with fdobj argument
gasmod1 = fpcr(gasoline$octane, fdobj = gas.fd, nc = 30)
# Method 2: Call fpcr with explicit signal matrix
gasmod2 = fpcr(gasoline$octane, sigmat = gasoline$NIR, nc = 30)
# Method 3: Call fpcr with explicit signal, basis, and penalty matrices
gasmod3 = fpcr(gasoline$octane, sigmat = gasoline$NIR, 
               basismat = eval.basis(wavelengths, bbasis), 
               penmats = list(getbasispenalty(bbasis)), nc = 30)

plot(gasmod1, xlab='Wavelength')

# Check that all 3 calls yield essentially identical estimates
all.equal(gasmod1$fhat, gasmod2$fhat, gasmod3$fhat)
# But note that, in general, you'd have to specify argvals in Method 1  
# to get the same coefficient function values as with Methods 2 & 3.
}

\keyword{models}
\keyword{regression}
\keyword{smooth}

