% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fpca.face.R
\name{fpca.face}
\alias{fpca.face}
\title{Functional principal component analysis with fast covariance estimation}
\usage{
fpca.face(
  Y = NULL,
  ydata = NULL,
  Y.pred = NULL,
  argvals = NULL,
  pve = 0.99,
  npc = NULL,
  var = FALSE,
  simul = FALSE,
  sim.alpha = 0.95,
  center = TRUE,
  knots = 35,
  p = 3,
  m = 2,
  lambda = NULL,
  alpha = 1,
  search.grid = TRUE,
  search.length = 100,
  method = "L-BFGS-B",
  lower = -20,
  upper = 20,
  control = NULL,
  periodicity = FALSE
)
}
\arguments{
\item{Y, ydata}{the user must supply either \code{Y}, a matrix of functions
observed on a regular grid, or a data frame \code{ydata} representing
irregularly observed functions. See Details.}

\item{Y.pred}{if desired, a matrix of functions to be approximated using
the FPC decomposition.}

\item{argvals}{numeric; function argument.}

\item{pve}{proportion of variance explained: used to choose the number of
principal components.}

\item{npc}{how many smooth SVs to try to extract, if \code{NA} (the
default) the hard thresholding rule of Gavish and Donoho (2014) is used (see
Details, References).}

\item{var}{logical; should an estimate of standard error be returned?}

\item{simul}{logical; if \code{TRUE} curves will we simulated using
Monte Carlo to obtain an estimate of the \code{sim.alpha} quantile at each
\code{argval}; ignored if \code{var == FALSE}}

\item{sim.alpha}{numeric; if \code{simul==TRUE}, quantile to estimate at
each \code{argval}; ignored if \code{var == FALSE}}

\item{center}{logical; center \code{Y} so that its column-means are 0? Defaults to
\code{TRUE}}

\item{knots}{number of knots to use or the vectors of knots; defaults to 35}

\item{p}{integer; the degree of B-splines functions to use}

\item{m}{integer; the order of difference penalty to use}

\item{lambda}{smoothing parameter; if not specified smoothing parameter is
chosen using \code{\link[stats]{optim}} or a grid search}

\item{alpha}{numeric; tuning parameter for GCV; see parameter \code{gamma}
in \code{\link[mgcv]{gam}}}

\item{search.grid}{logical; should a grid search be used to find \code{lambda}?
Otherwise, \code{\link[stats]{optim}} is used}

\item{search.length}{integer; length of grid to use for grid search for
\code{lambda}; ignored if \code{search.grid} is \code{FALSE}}

\item{method}{method to use; see \code{\link[stats]{optim}}}

\item{lower}{see \code{\link[stats]{optim}}}

\item{upper}{see \code{\link[stats]{optim}}}

\item{control}{see \code{\link[stats]{optim}}}

\item{periodicity}{Option for a periodic spline basis. Defaults to FALSE.}
}
\value{
A list with components
\enumerate{
\item \code{Yhat} - If \code{Y.pred} is specified, the smooth version of
\code{Y.pred}.   Otherwise, if \code{Y.pred=NULL}, the smooth version of \code{Y}.
\item \code{scores} - matrix of scores
\item \code{mu} - mean function
\item \code{npc} - number of principal components
\item \code{efunctions} - matrix of eigenvectors
\item \code{evalues} - vector of eigenvalues
\item \code{pve} - The percent variance explained by the returned number of PCs
}
if \code{var == TRUE} additional components are returned
\enumerate{
\item \code{sigma2} - estimate of the error variance
\item \code{VarMats} - list of covariance function estimate for each
subject
\item \code{diag.var} - matrix containing the diagonals of each matrix in
\item \code{crit.val} - list of estimated quantiles; only returned if
\code{simul == TRUE}
}
}
\description{
A fast implementation of the sandwich smoother (Xiao et al., 2013)
for covariance matrix smoothing. Pooled generalized cross validation
at the data level is used for selecting the smoothing parameter.
}
\examples{
#### settings
I <- 50 # number of subjects
J <- 3000 # dimension of the data
t <- (1:J)/J # a regular grid on [0,1]
N <- 4 #number of eigenfunctions
sigma <- 2 ##standard deviation of random noises
lambdaTrue <- c(1,0.5,0.5^2,0.5^3) # True eigenvalues

case = 1
### True Eigenfunctions

if(case==1) phi <- sqrt(2)*cbind(sin(2*pi*t),cos(2*pi*t),
                                sin(4*pi*t),cos(4*pi*t))
if(case==2) phi <- cbind(rep(1,J),sqrt(3)*(2*t-1),
                          sqrt(5)*(6*t^2-6*t+1),
                         sqrt(7)*(20*t^3-30*t^2+12*t-1))

###################################################
########     Generate Data            #############
###################################################
xi <- matrix(rnorm(I*N),I,N);
xi <- xi \%*\% diag(sqrt(lambdaTrue))
X <- xi \%*\% t(phi); # of size I by J
Y <- X + sigma*matrix(rnorm(I*J),I,J)

results <- fpca.face(Y,center = TRUE, argvals=t,knots=100,pve=0.99)

# calculate percent variance explained by each PC
 evalues = results$evalues
 pve_vec = evalues * results$npc/sum(evalues)

###################################################
####               FACE                ########
###################################################
Phi <- results$efunctions
eigenvalues <- results$evalues

for(k in 1:N){
  if(Phi[,k] \%*\% phi[,k]< 0)
    Phi[,k] <- - Phi[,k]
}

### plot eigenfunctions
par(mfrow=c(N/2,2))
seq <- (1:(J/10))*10
for(k in 1:N){
  plot(t[seq],Phi[seq,k]*sqrt(J),type="l",lwd = 3,
       ylim = c(-2,2),col = "red",
       ylab = paste("Eigenfunction ",k,sep=""),
       xlab="t",main="FACE")

  lines(t[seq],phi[seq,k],lwd = 2, col = "black")
}
}
\references{
Xiao, L., Li, Y., and Ruppert, D. (2013).
Fast bivariate \emph{P}-splines: the sandwich smoother,
\emph{Journal of the Royal Statistical Society: Series B}, 75(3), 577-599.

Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
Fast covariance estimation for high-dimensional functional data.
\emph{Statistics and Computing}, 26, 409-421.
DOI: 10.1007/s11222-014-9485-x.
}
\seealso{
\code{\link{fpca.sc}}  for another covariance-estimate based
smoothing of \code{Y}; \code{\link{fpca2s}} and \code{\link{fpca.ssvd}}
for two SVD-based smoothings.
}
\author{
Luo Xiao
}
