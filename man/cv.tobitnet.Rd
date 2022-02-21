\name{cv.tobitnet}
\alias{cv.tobitnet}
\title{Cross-validation for tobitnet}
\description{
  Performs k-fold cross-validation for the \code{lambda1} penalty parameter for \code{tobitnet}.
}
\usage{
cv.tobitnet(x, y, c = 0, lambda1 = NULL, nfolds = 10, early.stop = TRUE, ...)
}

\arguments{
  \item{x}{numeric predictor matrix with \code{n} rows and \code{nvars} columns, where each row corresponds to an observation and each column corresponds to a predictor}
  \item{y}{numeric response vector}
  \item{c}{left-censoring value for the response. Default is 0.}
  \item{lambda1}{vector of values for the lasso penalty parameter \eqn{\lambda_1}{\lambda1}. If this is not provided, then a \eqn{\lambda_1}{\lambda1} path is automatically generated based on the values of \code{nlambda} and \code{lambda.factor}}
  \item{nfolds}{number of folds. Default is 10.}
  \item{early.stop}{allow \code{tobitnet} to stop early in the \code{lambda1} sequence if the decrease in the deviance between \code{lambda1} values is negligible. Default is \code{TRUE}.}
  \item{...}{additional arguments for \code{tobitnet}}
}

\value{
  An object of class \dQuote{\code{cv.tobitnet}}
 \item{call}{function call}
 \item{cvm}{mean-squared prediction errors from cross-validation - a vector of length \code{length(lambda1)}}
 \item{cvsd}{estimate of the standard error of \code{cvm} - a vector of length \code{length(lambda1)}}
 \item{lambda1}{user-provided or automatically generated \code{lambda1} path}
 \item{lambda2}{user-provided \code{lambda2} value}
 \item{lambda1.min}{the value of \code{lambda1} that gives the minimum \code{cvm}}
 \item{lambda1.1se}{the largest value of \code{lambda1} such that \code{cvm} is within 1 standard error of the minimum.}
}

\details{
  Folds are created using stratified sampling so that the proportion of censored observations in each fold is roughly the same.
}

\author{
Tate Jacobson\cr
Maintainer: Tate Jacobson <jaco2583@umn.edu>
}

\seealso{\code{\link{tobitnet}}, \code{\link{predict.tobitnet}}}

\examples{
n <- 100
p <- 10
x <- vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
y <- 3 + x\%*\%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y <- pmax(y, 0)

#With a user-provided lambda1
cv.tnet1 <- cv.tobitnet(x = x, y = y, lambda1 = seq(0.1, 0.9, by = 0.1))

#With lambda1 automatically generated
cv.tnet2 <- cv.tobitnet(x = x, y = y)
}