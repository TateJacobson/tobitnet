\name{plot.cv.tobitscad}
\alias{plot.cv.tobitscad}
\title{Plot the cross-validation results from cv.tobitscad}
\description{
  Plots the cross-validation curve, with 1-standard error bars, from a \code{cv.tobitscad} object as a function of the \code{lambda} values used.
}
\usage{
  \method{plot}{cv.tobitscad}(x, ...)
}
\arguments{
  \item{x}{fitted \code{cv.tobitscad} object}
  \item{...}{other graphical parameters for \code{plot}}
}
\author{
Tate Jacobson\cr
Maintainer: Tate Jacobson <jaco2583@umn.edu>
}

\seealso{\code{\link{tobitscad}}, \code{\link{cv.tobitscad}}}

\examples{
  n = 100
  p = 10
  x = vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
  y = 3 + x\%*\%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
  y = pmax(y, 0)

  cv.tscad = cv.tobitscad(x = x, y = y)
  plot(cv.tscad)
}
