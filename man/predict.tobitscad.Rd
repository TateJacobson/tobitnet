\name{predict.tobitscad}
\alias{predict.tobitscad}
\title{Make predictions with a tobitscad object}

\description{
  This function makes predictions from a fitted \code{tobitscad} object.
}

\usage{
  \method{predict}{tobitscad}(object, newx, lambda = NULL, type = c("censored", "uncensored"), ...)
}

\arguments{
  \item{object}{fitted \code{tobitscad} object}
  \item{newx}{matrix of new predictor values at which to make predictions}
  \item{lambda}{value(s) of the penalty parameter \code{lambda} at which to make predictions. Default is the sequence used to fit the \code{tobitscad} object.}
  \item{type}{the type of prediction to be returned. \code{"censored"} predictions are left-censored at the \code{c} value used to fit the \code{tobitscad} object.}
  \item{...}{additional arguments to update the fitted \code{tobitscad} object.}
}

\value{
  Predicted values from the \code{tobitscad} object evaluated at \code{newx}.  A matrix with \code{nrow(newx)} rows and \code{length(lambda)} columns.
}

\author{
Tate Jacobson \cr
Maintainer: Tate Jacobson \email{tatejacobson1@gmail.com}
}

\seealso{ \code{\link{tobitscad}}, \code{\link{cv.tobitscad}} }

\examples{
n <- 100
p <- 10
x <- vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
y <- 3 + x\%*\%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y <- pmax(y, 0)

tscad <- tobitscad(x = x, y = y)

#Make predictions with lambda path from fitted tobitscad object
newx <- vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
predict(tscad, newx = newx)

#Make predictions with a new lambda path
predict(tscad, newx = newx, lambda = c(0.1, 0.25, 0.5))
}
