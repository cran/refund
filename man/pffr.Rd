% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pffr.R
\name{pffr}
\alias{pffr}
\title{Penalized flexible functional regression}
\usage{
pffr(
  formula,
  yind,
  data = NULL,
  ydata = NULL,
  algorithm = NA,
  method = "REML",
  tensortype = c("ti", "t2"),
  bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
  bs.int = list(bs = "ps", k = 20, m = c(2, 1)),
  ...
)
}
\arguments{
\item{formula}{a formula with special terms as for \code{\link[mgcv]{gam}},
with additional special terms \code{\link{ff}(), \link{sff}(),
\link{ffpc}(), \link{pcre}()} and \code{c()}.}

\item{yind}{a vector with length equal to the number of columns of the matrix
of functional responses giving the vector of evaluation points \eqn{(t_1,
\dots ,t_{G})}. If not supplied, \code{yind} is set to
\code{1:ncol(<response>)}.}

\item{data}{an (optional) \code{data.frame} containing the data. Can also be
a named list for regular data. Functional covariates have to be supplied as
<no. of observations> by <no. of evaluations> matrices, i.e. each row is
one functional observation.}

\item{ydata}{an (optional) \code{data.frame} supplying functional responses
that are not observed on a regular grid. See Details.}

\item{algorithm}{the name of the function used to estimate the model.
Defaults to \code{\link[mgcv]{gam}} if the matrix of functional responses
has less than \code{2e5} data points and to \code{\link[mgcv]{bam}} if not.
\code{'\link[mgcv]{gamm}'}, \code{'\link[gamm4]{gamm4}'} and
\code{'\link[mgcv]{jagam}'} are valid options as well. See Details for
\code{'\link[gamm4]{gamm4}'} and \code{'\link[mgcv]{jagam}'}.}

\item{method}{Defaults to \code{"REML"}-estimation, including of unknown
scale. If \code{algorithm="bam"}, the default is switched to
\code{"fREML"}. See \code{\link[mgcv]{gam}} and \code{\link[mgcv]{bam}} for
details.}

\item{tensortype}{which typ of tensor product splines to use. One of
"\code{\link[mgcv]{ti}}" or "\code{\link[mgcv]{t2}}", defaults to
\code{ti}. \code{t2}-type terms do not enforce the more suitable special
constraints for functional regression, see Details.}

\item{bs.yindex}{a named (!) list giving the parameters for spline bases on
the index of the functional response. Defaults to \code{list(bs="ps", k=5,
m=c(2, 1))}, i.e. 5 cubic B-splines bases with first order difference
penalty.}

\item{bs.int}{a named (!) list giving the parameters for the spline basis for
the global functional intercept. Defaults to \code{list(bs="ps", k=20,
m=c(2, 1))}, i.e. 20 cubic B-splines bases with first order difference
penalty.}

\item{...}{additional arguments that are valid for \code{\link[mgcv]{gam}},
\code{\link[mgcv]{bam}}, \code{'\link[gamm4]{gamm4}'} or
\code{'\link[mgcv]{jagam}'}. \code{subset} is not implemented.}
}
\value{
A fitted \code{pffr}-object, which is a
  \code{\link[mgcv]{gam}}-object with some additional information in an
  \code{pffr}-entry. If \code{algorithm} is \code{"gamm"} or \code{"gamm4"},
  only the \code{$gam} part of the returned list is modified in this way.\cr
  Available methods/functions to postprocess fitted models:
  \code{\link{summary.pffr}}, \code{\link{plot.pffr}},
  \code{\link{coef.pffr}}, \code{\link{fitted.pffr}},
  \code{\link{residuals.pffr}}, \code{\link{predict.pffr}},
  \code{\link{model.matrix.pffr}},  \code{\link{qq.pffr}},
  \code{\link{pffr.check}}.\cr If \code{algorithm} is \code{"jagam"}, only
  the location of the model file and the usual
  \code{\link[mgcv]{jagam}}-object are returned, you have to run the sampler
  yourself.\cr
}
\description{
Implements additive regression for functional and scalar covariates and
functional responses. This function is a wrapper for \code{mgcv}'s
\code{\link[mgcv]{gam}} and its siblings to fit models of the general form
\cr \eqn{E(Y_i(t)) = g(\mu(t) + \int X_i(s)\beta(s,t)ds + f(z_{1i}, t) +
f(z_{2i}) + z_{3i} \beta_3(t) + \dots )}\cr with a functional (but not
necessarily continuous) response \eqn{Y(t)}, response function \eqn{g},
(optional) smooth intercept \eqn{\mu(t)}, (multiple) functional covariates
\eqn{X(t)} and scalar covariates \eqn{z_1}, \eqn{z_2}, etc.
}
\section{Details}{
 The routine can estimate \enumerate{ \item linear
  functional effects of scalar (numeric or factor) covariates that vary
  smoothly over \eqn{t} (e.g. \eqn{z_{1i} \beta_1(t)}, specified as
  \code{~z1}), \item nonlinear, and possibly multivariate functional effects
  of (one or multiple) scalar covariates \eqn{z} that vary smoothly over the
  index \eqn{t} of \eqn{Y(t)} (e.g. \eqn{f(z_{2i}, t)}, specified in the
  \code{formula} simply as \code{~s(z2)}) \item (nonlinear) effects of scalar
  covariates that are constant over \eqn{t} (e.g. \eqn{f(z_{3i})}, specified
  as \code{~c(s(z3))}, or \eqn{\beta_3 z_{3i}}, specified as \code{~c(z3)}),
  \item function-on-function regression terms (e.g. \eqn{\int
  X_i(s)\beta(s,t)ds}, specified as \code{~ff(X, yindex=t, xindex=s)}, see
  \code{\link{ff}}). Terms given by \code{\link{sff}} and \code{\link{ffpc}}
  provide nonlinear and FPC-based effects of functional covariates,
  respectively. \item concurrent effects of functional covariates \code{X}
  measured on the same grid as the response  are specified as follows:
  \code{~s(x)} for a smooth, index-varying effect \eqn{f(X(t),t)}, \code{~x}
  for a linear index-varying effect \eqn{X(t)\beta(t)}, \code{~c(s(x))} for a
  constant nonlinear effect \eqn{f(X(t))}, \code{~c(x)} for a constant linear
  effect \eqn{X(t)\beta}. \item Smooth functional random intercepts
  \eqn{b_{0g(i)}(t)} for a grouping variable \code{g} with levels \eqn{g(i)}
  can be specified via \code{~s(g, bs="re")}), functional random slopes
  \eqn{u_i b_{1g(i)}(t)} in a numeric variable \code{u} via \code{~s(g, u,
  bs="re")}). Scheipl, Staicu, Greven (2013) contains code examples for
  modeling correlated functional random intercepts using
  \code{\link[mgcv]{mrf}}-terms. } Use the \code{c()}-notation to denote
  model terms that are constant over the index of the functional response.\cr

  Internally, univariate smooth terms without a \code{c()}-wrapper are
  expanded into bivariate smooth terms in the original covariate and the
  index of the functional response. Bivariate smooth terms (\code{s(), te()}
  or \code{t2()}) without a \code{c()}-wrapper are expanded into trivariate
  smooth terms in the original covariates and the index of the functional
  response. Linear terms for scalar covariates or categorical covariates are
  expanded into varying coefficient terms, varying smoothly over the index of
  the functional response. For factor variables, a separate smooth function
  with its own smoothing parameter is estimated for each level of the
  factor.\cr \cr The marginal spline basis used for the index of the the
  functional response is specified via the \emph{global} argument
  \code{bs.yindex}. If necessary, this can be overriden for any specific term
  by supplying a \code{bs.yindex}-argument to that term in the formula, e.g.
  \code{~s(x, bs.yindex=list(bs="tp", k=7))} would yield a tensor product
  spline over \code{x} and the index of the response in which the marginal
  basis for the index of the response are 7 cubic thin-plate spline functions
  (overriding the global default for the basis and penalty on the index of
  the response given by the \emph{global} \code{bs.yindex}-argument).\cr Use
  \code{~-1 + c(1) + ...} to specify a model with only a constant and no
  functional intercept. \cr

  The functional covariates have to be supplied as a \eqn{n} by <no. of
  evaluations> matrices, i.e. each row is one functional observation. For
  data on a regular grid, the functional response is supplied in the same
  format, i.e. as a matrix-valued entry in \code{data},  which can contain
  missing values.\cr

  If the functional responses are \emph{sparse or irregular} (i.e., not
  evaluated on the same evaluation points across all observations), the
  \code{ydata}-argument can be used to specify the responses: \code{ydata}
  must be a \code{data.frame} with 3 columns called \code{'.obs', '.index',
  '.value'} which specify which curve the point belongs to
  (\code{'.obs'}=\eqn{i}), at which \eqn{t} it was observed
  (\code{'.index'}=\eqn{t}), and the observed value
  (\code{'.value'}=\eqn{Y_i(t)}). Note that the vector of unique sorted
  entries in \code{ydata$.obs} must be equal to \code{rownames(data)} to
  ensure the correct association of entries in \code{ydata} to the
  corresponding rows of \code{data}. For both regular and irregular
  functional responses, the model is then fitted with the data in long
  format, i.e., for data on a grid the rows of the matrix of the functional
  response evaluations \eqn{Y_i(t)} are stacked into one long vector and the
  covariates are expanded/repeated correspondingly. This means the models get
  quite big fairly fast, since the effective number of rows in the design
  matrix is number of observations times number of evaluations of \eqn{Y(t)}
  per observation.\cr

  Note that \code{pffr} does not use \code{mgcv}'s default identifiability
  constraints (i.e., \eqn{\sum_{i,t} \hat f(z_i, x_i, t) = 0} or
  \eqn{\sum_{i,t} \hat f(x_i, t) = 0}) for tensor product terms whose
  marginals include the index \eqn{t} of the functional response.  Instead,
  \eqn{\sum_i \hat f(z_i, x_i, t) = 0} for all \eqn{t} is enforced, so that
  effects varying over \eqn{t} can be interpreted as local deviations from
  the global functional intercept. This is achieved by using
  \code{\link[mgcv]{ti}}-terms with a suitably modified \code{mc}-argument.
  Note that this is not possible if \code{algorithm='gamm4'} since only
  \code{t2}-type terms can then be used and these modified constraints are
  not available for \code{t2}. We recommend using centered scalar covariates
  for terms like \eqn{z \beta(t)} (\code{~z}) and centered functional
  covariates with \eqn{\sum_i X_i(t) = 0} for all \eqn{t} in \code{ff}-terms
  so that the global functional intercept can be interpreted as the global
  mean function.

  The \code{family}-argument can be used to specify all of the response
  distributions and link functions described in
  \code{\link[mgcv]{family.mgcv}}. Note that  \code{family = "gaulss"} is
  treated in a special way: Users can supply the formula for the variance by
  supplying a special argument \code{varformula}, but this is not modified in
  the way that the \code{formula}-argument is but handed over to the fitter
  directly, so this is for expert use only. If \code{varformula} is not
  given, \code{pffr} will use the parameters from argument \code{bs.int} to
  define a spline basis along the index of the response, i.e., a smooth
  variance function over $t$ for responses $Y(t)$.
}

\examples{
###############################################################################
# univariate model:
# Y(t) = f(t)  + \int X1(s)\beta(s,t)ds + eps
set.seed(2121)
data1 <- pffrSim(scenario="ff", n=40)
t <- attr(data1, "yindex")
s <- attr(data1, "xindex")
m1 <- pffr(Y ~ ff(X1, xind=s), yind=t, data=data1)
summary(m1)
plot(m1, pages=1)

\dontrun{
###############################################################################
# multivariate model:
# E(Y(t)) = \beta_0(t)  + \int X1(s)\beta_1(s,t)ds + xlin \beta_3(t) +
#        f_1(xte1, xte2) + f_2(xsmoo, t) + \beta_4 xconst
data2 <- pffrSim(scenario="all", n=200)
t <- attr(data2, "yindex")
s <- attr(data2, "xindex")
m2 <- pffr(Y ~  ff(X1, xind=s) + #linear function-on-function
                xlin  +  #varying coefficient term
                c(te(xte1, xte2)) + #bivariate smooth term in xte1 & xte2, const. over Y-index
                s(xsmoo) + #smooth effect of xsmoo varying over Y-index
                c(xconst), # linear effect of xconst constant over Y-index
        yind=t,
        data=data2)
summary(m2)
plot(m2)
str(coef(m2))
# convenience functions:
preddata <- pffrSim(scenario="all", n=20)
str(predict(m2, newdata=preddata))
str(predict(m2, type="terms"))
cm2 <- coef(m2)
cm2$pterms
str(cm2$smterms, 2)
str(cm2$smterms[["s(xsmoo)"]]$coef)

#############################################################################
# sparse data (80\% missing on a regular grid):
set.seed(88182004)
data3 <- pffrSim(scenario=c("int", "smoo"), n=100, propmissing=0.8)
t <- attr(data3, "yindex")
m3.sparse <- pffr(Y ~ s(xsmoo), data=data3$data, ydata=data3$ydata, yind=t)
summary(m3.sparse)
plot(m3.sparse,pages=1)
}
}
\references{
Ivanescu, A., Staicu, A.-M., Scheipl, F. and Greven, S. (2015).
  Penalized function-on-function regression. Computational Statistics,
  30(2):539--568. \url{https://biostats.bepress.com/jhubiostat/paper254/}

  Scheipl, F., Staicu, A.-M. and Greven, S. (2015). Functional Additive Mixed
  Models. Journal of Computational & Graphical Statistics, 24(2): 477--501.
  \url{ https://arxiv.org/abs/1207.5947}

  F. Scheipl, J. Gertheiss, S. Greven (2016):  Generalized Functional Additive Mixed Models,
  Electronic Journal of Statistics, 10(1), 1455--1492.
  \url{https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-10/issue-1/Generalized-functional-additive-mixed-models/10.1214/16-EJS1145.full}
}
\seealso{
\code{\link[mgcv]{smooth.terms}} for details of \code{mgcv} syntax
  and available spline bases and penalties.
}
\author{
Fabian Scheipl, Sonja Greven
}
