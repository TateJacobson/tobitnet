\name{tobitscad}
\alias{tobitscad}
\title{Fit the Tobit model with a SCAD penalty}

\description{
Fits the Tobit model with a SCAD penalty using a local linear approximation (LLA) algorithm.  The solution path is computed for a sequence of values for the penalty parameter \eqn{\lambda}.
}

\usage{
tobitscad(x, y, c = 0, a = 3 , iter = 3, nlambda = 100,
            lambda.factor = ifelse(n/2 < nvars, 0.1, 0.05),
            lambda = NULL, eps = 1e-7, standardize = TRUE,
            maxit = 1e6, early.stop = TRUE)

}

\arguments{
  \item{x}{numeric predictor matrix with \code{n} rows and \code{nvars} columns, where each row corresponds to an observation and each column corresponds to a predictor}
  \item{y}{numeric response vector}
  \item{c}{left-censoring value for the response. Default is 0.}
  \item{a}{the tuning parameter for the SCAD penalty. Default is 3.}
  \item{iter}{the number of iterations for the LLA algorithm. Default is 3.}
  \item{nlambda}{if \code{lambda} is not provided, this is the number of \eqn{\lambda} values used for the solution path. Default is 100.}
  \item{lambda.factor}{if \code{lambda} is not provided, this is the factor multiplied by \eqn{\lambda_{\max}}{\lambdamax}, the largest \eqn{\lambda} value in the automatically generated sequence for the solution path, to get \eqn{\lambda_{\min}}{\lambdamin}, the smallest \eqn{\lambda} value for the solution path. Default is 0.1 if \code{n/2 < nvars} and 0.05 otherwise.}
  \item{lambda}{vector of values for the SCAD penalty parameter \eqn{\lambda}. If this is not provided, then a \eqn{\lambda} path is automatically generated based on the values of \code{nlambda} and \code{lambda.factor}.}
  \item{eps}{threshold for the update step size in coordinate descent at which to declare convergence. Default is 1e-7.}
  \item{standardize}{boolean, if \code{TRUE} then \code{x} is standardized before fitting. Default is \code{TRUE}.}
  \item{maxit}{maximum number of coordinate descent cycles to attempt. Default is 1e6.}
  \item{early.stop}{allow \code{tobitscad} to stop early in the \code{lambda} sequence if the decrease in the deviance between \code{lambda} values is negligible. Default is \code{TRUE}.}
}

\value{
An object with S3 class \dQuote{\code{tobitscad}}.
 \item{call}{function call}
 \item{sigma}{vector of estimated sigma values of length \code{length(lambda)}}
 \item{b0}{vector of estimated intercept values of length \code{length(lambda)}}
 \item{beta}{matrix of estimated coefficient values with \code{nvars} rows and \code{length(lambda)} columns}
 \item{c}{user-provided \code{c} argument}
 \item{lambda}{user-provided or automatically generated \code{lambda} path}
 \item{dev}{model deviance}
 \item{nulldev}{null deviance (from the intercept-only model)}
}

\details{
  The \code{tobitscad} function fits the Tobit model with a SCAD penalty using a local linear approximation (LLA) algorithm.
  
  \subsection{Handling the predictors}{
    By default, \code{tobitscad} standardizes the predictors in \code{x} so that their coefficients are all equally penalized.  
    If a predictor has zero variance, we set its coefficient to zero and exclude it from the model fitting procedure.
  }
  
  \subsection{The Tobit model with a SCAD penalty}{
    To model the censored response \eqn{y}, we assume that there is a latent uncensored response  \eqn{y^*}{y*} and that it comes from a linear model.  That is, we assume \eqn{y = \max\{y^*,c\}}{y = max{y*,c}} where \eqn{y^* = x' \beta + \epsilon}{y* = x' \beta + \epsilon}, \eqn{\epsilon \sim N(0, \sigma^2)}{\epsilon ~ N(0, \sigma^2)}.  
We use Olsen's (1978) reparameterization of the Tobit log likelihood, with \eqn{\delta_j = \beta_j/\sigma}{\deltaj = \betaj/\sigma} for \eqn{j = 0, 1, \ldots, p} and \eqn{\gamma = 1/\sigma}. We aim to minimize the following objective function:
    \deqn{R_n(\delta, \gamma) = - \frac{1}{n} \log L_n(\delta, \gamma) + P_{\lambda}(\delta) }{
    Rn(\delta, \gamma) = - 1/n log Ln(\delta, \gamma) + P\lambda(\delta) } 
    where \eqn{L_n(\delta, \gamma)}{Ln(\delta, \gamma)} denotes the reparameterized Tobit likelihood and \eqn{ P_{\lambda}(\delta) }{ P\lambda(\delta) } denotes the SCAD penalty. 
  }
  
  \subsection{LLA algorithm details}{
    Because \eqn{ P_{\lambda}(\delta) }{ P\lambda(\delta) } is a folded concave penalty, \eqn{R_n(\delta, \gamma)}{Rn(\delta, \gamma)} is not convex. As such, we use a local linear approximation (LLA) algorithm to find a local minimizer of \eqn{R_n(\delta, \gamma)}{Rn(\delta, \gamma)}.  
    In the \eqn{m}th LLA step, we aim to minimize the following objective function
      \deqn{ - \frac{1}{n} \log L_n(\delta, \gamma) + \sum_{j = 1}^p \hat{w}_{j}^{(m-1)} |\delta_j|  }{
      - 1/n log Ln(\delta, \gamma) + \sum wj(m-1) |\deltaj|
      }
    where \eqn{\hat{w}_{j}^{(m-1)} = P_{\lambda}'(|\hat{\delta}_j^{(m-1)}|)}{wj(m-1) = P\lambda'(|\deltaj(m-1)|)} and \eqn{ \hat{\delta}_j^{(m-1)} }{\deltaj(m-1)} denotes the estimate of \eqn{ \delta_j }{\deltaj} from the previous LLA step.  Note that we initialize the LLA algorithm with \eqn{\hat{\delta}_j^{(0)} = 0}{\deltaj(0) = 0} for \eqn{j = 1, \ldots, p}.  
    We use the generalized coordinate descent (GCD) algorithm for \code{tobitnet} to minimize the objective in each LLA step.  
    See the \code{\link{tobitnet}} documentation for details on the GCD algorithm, specifically the warning issued if the KKT conditions are not satisfied and \code{early.stop} behavior.
  }
  
  \subsection{Parameter estimates}{
    \code{tobitscad} returns \eqn{\beta} and \eqn{\sigma} rather than \eqn{\delta} and \eqn{\gamma} to faciliate comparison to parameter estimates from other regression models. 
  }

}

\author{
Tate Jacobson\cr
Maintainer: Tate Jacobson <jaco2583@umn.edu>
}

\references{
Olsen, R.J. (1978) Note on the Uniqueness of the Maximum Likelihood Estimator for the
Tobit Model. \emph{Econometrica}, \bold{46(5)}, 1211--1215.

Jacobson, T. and Zou, H. (2022) \emph{High-dimensional Censored Regression via the Penalized Tobit Likelihood.} Manuscript submitted for publication.
}

\seealso{\code{\link{tobitnet}}, \code{\link{predict.tobitscad}}, \code{\link{plot.tobitscad}}, and \code{\link{cv.tobitscad}}}

\examples{
n <- 100
p <- 10
x <- vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
y <- 3 + x\%*\%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y <- pmax(y, 0)

#With a user-provided lambda
tscad1 <- tobitscad(x = x, y = y, lambda = 0.05)

#With lambda1 automatically generated
tscad2 <- tobitscad(x = x, y = y)
}
