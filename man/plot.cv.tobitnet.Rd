\name{plot.cv.tobitnet}
\alias{plot.cv.tobitnet}
\title{Plot the cross-validation results from cv.tobitnet}
\description{
  Plots the cross-validation curve, with 1-standard error bars, from a \code{cv.tobitnet} object as a function of the \code{lambda1} values used.
}
\usage{
  \method{plot}{cv.tobitnet}(x, ...)
}
\arguments{
  \item{x}{fitted \code{cv.tobitnet} object}
  \item{...}{other graphical parameters for \code{plot}}
}
\author{
Tate Jacobson\cr
Maintainer: Tate Jacobson <jaco2583@umn.edu>
}

\seealso{\code{\link{tobitnet}}, \code{\link{cv.tobitnet}}}

\examples{
n <- 100
p <- 10
x <- vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
y <- 3 + x\%*\%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y <- pmax(y, 0)

cv.tnet <- cv.tobitnet(x = x, y = y)
plot(cv.tnet)
}
