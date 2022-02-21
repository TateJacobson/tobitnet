\name{cv.tobitscad}
\alias{cv.tobitscad}
\title{Cross-validation for tobitscad}
\description{
  Performs k-fold cross-validation for the \code{lambda} penalty parameter for \code{tobitscad}.
}
\usage{
cv.tobitscad(x, y, c = 0, a = 3, iter = 3, nlambda = 100, 
  lambda.factor = ifelse(n/2 < nvars, 0.1, 0.05), lambda = NULL, 
  nfolds = 10, early.stop = TRUE,...)
}
\arguments{
  \item{x}{numeric predictor matrix with \code{n} rows and \code{nvars} columns, where each row corresponds to an observation and each column corresponds to a predictor}
  \item{y}{numeric response vector}
  \item{c}{left-censoring value for the response. Default is 0.}
  \item{a}{the tuning parameter for the SCAD penalty. Default is 3.}
  \item{iter}{the number of iterations for the LLA algorithm. Default is 3.}
  \item{nlambda}{if \code{lambda} is not provided, this is the number of \eqn{\lambda} values used for the solution path. Default is 100.}
  \item{lambda.factor}{if \code{lambda} is not provided, this is the factor multiplied by \eqn{\lambda_{\max} }{\lambdamax}, the largest \eqn{\lambda} value in the automatically generated sequence for the solution path, to get \eqn{\lambda_{\min}}{\lambdamin}, the smallest \eqn{\lambda} value for the solution path. Default is 0.1 if \code{n/2 < nvars} and 0.05 otherwise.}
  \item{lambda}{vector of values for the SCAD penalty parameter \eqn{\lambda}. If this is not provided, then a \eqn{\lambda} path is automatically generated based on the values of \code{nlambda} and \code{lambda.factor}.}
  \item{nfolds}{number of folds. Default is 10.}
  \item{early.stop}{allow \code{tobitscad} to stop early in the \code{lambda} sequence if the decrease in the deviance between \code{lambda} values is negligible. Default is \code{TRUE}.}
  \item{...}{additional arguments for \code{tobitscad}}
}
\value{
  An object of class \dQuote{\code{cv.tobitscad}}
 \item{call}{function call}
 \item{cvm}{mean-squared prediction errors from cross-validation - a vector of length \code{length(lambda)}}
 \item{cvsd}{estimate of the standard error of \code{cvm} - a vector of length \code{length(lambda)}}
 \item{lambda}{user-provided or automatically generated \code{lambda} path}
 \item{lambda.min}{the value of \code{lambda} that gives the minimum \code{cvm}}
 \item{lambda.1se}{the largest value of \code{lambda} such that \code{cvm} is within 1 standard error of the minimum.}
}

\details{
  Folds are created using stratified sampling so that the proportion of censored observations in each fold is roughly the same.
}

\author{
Tate Jacobson\cr
Maintainer: Tate Jacobson <jaco2583@umn.edu>
}

\seealso{\code{\link{tobitscad}}, \code{\link{predict.tobitscad}}}

\examples{
n <- 100
p <- 10
x <- vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
y <- 3 + x\%*\%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y <- pmax(y, 0)

#With a user-provided lambda1
cv.tscad1 <- cv.tobitscad(x = x, y = y, lambda = seq(0.1, 0.9, by = 0.1))

#With lambda1 automatically generated
cv.tscad2 <- cv.tobitscad(x = x, y = y)
}