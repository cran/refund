\name{pffr}
\alias{pffr}
\title{Penalized function-on-function regression}
\usage{
  pffr(formula, yind, fitter = NA, method = "REML",
  bsy.default = list(bs = "ps", m = c(2, 1)), ...)
}
\arguments{
  \item{formula}{a formula with special terms as for
  \code{\link[mgcv]{gam}}, with additional special terms
  \code{\link{ff}()} and \code{c()}}

  \item{yind}{a vector with length equal to the number of
  columns of the matrix of functional responses giving the
  vector of evaluation points \eqn{(t_1, \dots ,t_{G})}. If
  \code{formula} contains an \code{\link{ff}}-term which
  specifies \code{yind} this is used. If neither is given,
  \code{yind} is \code{1:ncol(<response>)}.}

  \item{fitter}{the name of the function used to estimate
  the model. Defaults to \code{\link[mgcv]{gam}} if the
  matrix of functional responses has less than \code{2e5}
  data points and to \code{\link[mgcv]{bam}} if not. "gamm"
  (see \code{\link[mgcv]{gamm}}) and "gamm4" (see
  \code{\link[gamm4]{gamm4}}) are valid options as well.}

  \item{method}{Defaults to \code{"REML"}-estimation,
  including of unknown scale. See \code{\link[mgcv]{gam}}
  for details.}

  \item{bsy.default}{a named (!) list giving the parameters
  for spline bases on the index of the functional response.
  Defaults to \code{list(bs="ps", m=c(2, 1))}, i.e. a cubic
  B-spline basis with first order difference penalty. Only
  arguments \code{bs, k, m} are used, see
  \code{\link[mgcv]{s}} for details.}

  \item{...}{additional arguments that are valid for
  \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}.
  \code{weights, subset, offset} are not yet implemented!}
}
\value{
  a fitted \code{pffr}-object, which is a
  \code{\link[mgcv]{gam}}-object with some additional
  information in an \code{pffr}-entry. If \code{fitter} is
  \code{"gamm"} or \code{"gamm4"}, only the \code{$gam}
  part of the returned list is modified in this way.
}
\description{
  Implements additive regression for functional and scalar
  covariates and functional responses. This function is a
  wrapper for \code{mgcv}'s \code{\link[mgcv]{gam}} and its
  siblings to fit models of the general form \cr
  \eqn{E(Y_i(t)) = g(\mu(t) + \int X_i(s)\beta(s,t)ds +
  f(z_{1i}, t) + f(z_{2i}) + z_{3i} \beta_3(t) + \dots}
  )\cr with a functional (but not necessarily continuous)
  response \eqn{Y(t)}, response function \eqn{g},
  (optional) smooth intercept \eqn{\mu(t)}, (multiple)
  functional covariates \eqn{X(t)} and scalar covariates
  \eqn{z_1}, \eqn{z_2}, etc.
}
\section{Details}{
  The routine can estimate \enumerate{ \item (nonlinear,
  and possibly multivariate) effects of (one or multiple)
  scalar covariates \eqn{z} that vary smoothly over the
  index \eqn{t} of \eqn{Y(t)} (e.g. \eqn{f(z_{1i}, t)},
  specified in the \code{formula} simply as \code{~s(z1)}),
  \item (nonlinear) effects of scalar covariates that are
  constant over \eqn{t} (e.g. \eqn{f(z_{2i})}, specified as
  \code{~c(s(z2))}, or \eqn{\beta_2 z_{2i}}, specified as
  \code{~c(z2)}), \item linear functional effects of scalar
  (numeric or factor) covariates that vary smoothly over
  \eqn{t} (e.g. \eqn{z_{3i} \beta_3(t)}, specified as
  \code{~z3}), \item function-on-function regression terms
  (e.g. \eqn{\int X_i(s)\beta(s,t)ds}, specified as
  \code{~ff(X, yindex=t, xindex=s)}, see \code{\link{ff}}).
  } Use the \code{c()}-notation to denote model terms that
  are constant over the index of the functional
  response.\cr

  Internally, univariate smooth terms without a
  \code{c()}-wrapper are expanded into bivariate smooth
  terms in the original covariate and the index of the
  functional response.  Bivariate smooth terms (\code{s(),
  te()} or \code{t2()}) without a \code{c()}-wrapper are
  expanded into trivariate smooth terms in the original
  covariates and the index of the functional response.
  Linear terms for scalar covariates or categorical
  covariates are expanded into varying coefficient terms,
  varying smoothly over the index of the functional
  response.\cr Use \code{~-1 + c(1) + ...} to specify a
  model with a constant instead of a functional intercept.
  \cr

  The functional response and functional covariates have to
  be supplied as n by <no. of evaluations> matrices, i.e.
  each row is one functional observation. The model is then
  fitted with the data in long format, i.e., the rows of
  the matrix of the functional response evaluations
  \eqn{Y_i(t)} are stacked into one long vector and the
  covariates are expanded/repeated correspondingly. This
  means the models get quite big fairly fast, since the
  effective number of rows in the design matrix is number
  of observations times number of evaluations of \eqn{Y(t)}
  per observation.\cr

  Note that only default identifiability constraints as
  implemented in \code{\link[mgcv]{gam.side}} are used,
  i.e. \cr \eqn{\sum_{i,t} \hat f(z_i, t) = 0}. (\code{\link[mgcv]{gam.side}}
  may impose other identifiability constraints as
  well but doesn't give interpretable results.) These may
  not be well suited for models with multiple functional
  effects varying across \eqn{t}.\cr Centering estimated
  effects \eqn{\hat f(z_i, t)} (\code{~s(z)}) across
  observations so that \eqn{\sum_i \hat f(z_i, t) = 0} for
  all \eqn{t}, i.e., transforming \eqn{\hat f(z_i, t)
  \rightarrow \hat f(z_i, t) - n^{-1} \sum_i \hat f(z_i,
  t)} and adding the mean function \eqn{n^{-1} \sum_i \hat
  f(z_i, t)} to the global functional intercept can improve
  the interpretability of results. We also recommend using
  centered scalar covariates for terms like \eqn{z
  \beta(t)} (\code{~z}) and centered functional covariates
  with \eqn{\sum_i X_i(t) = 0} for all \eqn{t} in
  \code{ff}-terms.
}
\author{
  Fabian Scheipl, Sonja Greven
}
\references{
Ivanescu, A. E., Staicu, A.-M., Greven, S., Scheipl, F., and Crainiceanu, C. M. (2011).  Penalized function-on-function regression.  Submitted.
}
\examples{


data(DTI, package="refund")
DTI$case <- factor(DTI$case)
DTI$ID <- factor(DTI$ID)
DTI$visit <- factor(DTI$visit)
rcst.complete <- apply(DTI, 1, function(x) !any(is.na(x)))
DTIcomplete <- DTI[rcst.complete,]

## split into test and training data
set.seed(2213213)
trainind <- sort(sample(1:sum(rcst.complete), 150))
train <- DTIcomplete[trainind,]
test <-  DTIcomplete[-trainind,]
# drop subjects not in the train data
test <- test[test$ID \%in\% unique(train$ID),]

rcstind <- 1:55
# constant random intercepts for ID,
# smooth effect of pasat varying over index of rcst
# constant effect of visit
# m1 adds a functional effect for cca
summary(m0 <- pffr(rcst ~ c(s(ID, bs="re")) + s(pasat) + c(visit),
                yind=rcstind, data=train))
summary(m1 <- pffr(rcst ~ c(s(ID, bs="re")) + s(pasat) +
                c(visit) + ff(cca, yind=rcstind),
                yind=rcstind, data=train))
AIC(m0, m1)
plot(m1, pers=TRUE, pages=1, ticktype="detailed")
str(coefs.m1 <- coef(m1))

#get predicted trajectories for test data
pr.m0 <- predict(m0, newdata=test)
pr.m1 <- predict(m1, newdata=test)
# (riemann integrated) squared prediction error:
(IMSE <- c(m0 = sum((pr.m0-test$rcst)^2),
                    m1 = sum((pr.m1-test$rcst)^2)))
layout(t(1:3))
matplot(t(test$rcst), type="l", lty=1, col=rgb(0,0,0,.1), main="observed rcst",
          ylim=range(test$rcst, pr.m0, pr.m1))
matplot(t(pr.m0),type="l", lty=1, col=rgb(0,0,0,.1), main="predicted: m0",
          ylim=range(test$rcst, pr.m0, pr.m1))
matplot(t(pr.m1),type="l", lty=1, col=rgb(0,0,0,.1), main="predicted: m1",
          ylim=range(test$rcst, pr.m0, pr.m1))

}

