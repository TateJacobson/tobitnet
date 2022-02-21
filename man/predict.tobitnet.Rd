\name{predict.tobitnet}
\alias{predict.tobitnet}
\title{Make predictions with a tobitnet object}

\description{
  This function makes predictions from a fitted \code{tobitnet} object.
}

\usage{
  \method{predict}{tobitnet}(object, newx, lambda1 = NULL, type = c("censored", "uncensored"), ...)
}

\arguments{
  \item{object}{fitted \code{tobitnet} object}
  \item{newx}{matrix of new predictor values at which to make predictions}
  \item{lambda1}{value(s) of the penalty parameter \code{lambda1} at which to make predictions. Default is the sequence used to fit the \code{tobitnet} object.}
  \item{type}{the type of prediction to be returned. \code{"censored"} predictions are left-censored at the \code{c} value used to fit the \code{tobitnet} object.}
  \item{...}{additional arguments to update the fitted \code{tobitnet} object.}
}

\value{
  Predicted values from the \code{tobitnet} object evaluated at \code{newx}.  A matrix with \code{nrow(newx)} rows and \code{length(lambda1)} columns.
}

\author{
Tate Jacobson\cr
Maintainer: Tate Jacobson <jaco2583@umn.edu>
}

\seealso{ \code{\link{tobitnet}}, \code{\link{cv.tobitnet}} }
\examples{
n <- 100
p <- 10
x <- vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
y <- 3 + x\%*\%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y <- pmax(y, 0)

tnet <- tobitnet(x = x, y = y)

#Make predictions with lambda1 path from fitted tobitnet object
newx <- vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
predict(tnet, newx = newx)

#Make predictions with a new lambda1 path
predict(tnet, newx = newx, lambda1 = c(0.1, 0.25, 0.5))
}
