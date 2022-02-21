\name{tobitnet}
\alias{tobitnet}
\title{Fit the Tobit model with a weighted elastic net penalty}

\description{
Fits the Tobit model with a weighted elastic net penalty. The solution path is computed for a sequence of values for the lasso penalty parameter \eqn{\lambda_1}{\lambda1} given a user-provided value for the ridge penalty parameter \eqn{\lambda_2}{\lambda2}.
}

\usage{
tobitnet(x, y, c = 0, nlambda = 100, lambda.factor = ifelse(n/2 < nvars, 0.05, 0.01), 
          lambda1 = NULL, lambda2 = 0, pf1 = rep(1, nvars), pf2 = rep(1, nvars), 
          eps = 1e-7, standardize = TRUE, maxit = 1e6, early.stop = TRUE)
}

\arguments{
  \item{x}{numeric predictor matrix with \code{n} rows and \code{nvars} columns, where each row corresponds to an observation and each column corresponds to a predictor}
  \item{y}{numeric response vector}
  \item{c}{left-censoring value for the response. Default is 0.}
  \item{nlambda}{if \code{lambda1} is not provided, this is the number of \eqn{\lambda_1}{\lambda1} values used for the solution path. Default is 100.}
  \item{lambda.factor}{if \code{lambda1} is not provided, this is the factor multiplied by \eqn{\lambda_{\max}}{\lambdamax}, the largest \eqn{\lambda_1}{\lambda1} value in the automatically generated sequence for the solution path, to get \eqn{\lambda_{\min}}{\lambdamin}, the smallest \eqn{\lambda_1}{\lambda1} value for the solution path. Default is 0.05 if \code{n/2 < nvars} and 0.01 otherwise.}
  \item{lambda1}{vector of values for the lasso penalty parameter \eqn{\lambda_1}{\lambda1}. If this is not provided, then a \eqn{\lambda_1}{\lambda1} path is automatically generated based on the values of \code{nlambda} and \code{lambda.factor}.}
  \item{lambda2}{scalar value for the ridge penalty parameter \eqn{\lambda_2}{\lambda2}. Default is 0.}
  \item{pf1}{vector of \code{nvars} penalty factors (weights) for the lasso penalty. Default is 1 for each predictor.}
  \item{pf2}{vector of \code{nvars} penalty factors (weights) for the ridge penalty. Default is 1 for each predictor.}
  \item{eps}{threshold for the update step size in coordinate descent at which to declare convergence. Default is 1e-7.}
  \item{standardize}{boolean, if \code{TRUE} then \code{x} is standardized before fitting. Default is \code{TRUE}.}
  \item{maxit}{maximum number of coordinate descent cycles to attempt. Default is 1e6.}
  \item{early.stop}{allow \code{tobitnet} to stop early in the \code{lambda1} sequence if the decrease in the deviance between \code{lambda1} values is negligible. Default is \code{TRUE}.}
}

\value{
An object with S3 class \dQuote{\code{tobitnet}}.
 \item{call}{function call}
 \item{sigma}{vector of estimated sigma values of length \code{length(lambda1)}}
 \item{b0}{vector of estimated intercept values of length \code{length(lambda1)}}
 \item{beta}{matrix of estimated coefficient values with \code{nvars} rows and \code{length(lambda1)} columns}
 \item{c}{user-provided \code{c} argument}
 \item{lambda1}{user-provided or automatically generated \code{lambda1} path}
 \item{lambda2}{user-provided \code{lambda2} value}
 \item{dev}{model deviance}
 \item{nulldev}{null deviance (from the intercept-only model)}
}

\details{
The \code{tobitnet} function fits the Tobit model with a weighted elastic net penalty.

\subsection{Handling the predictors}{
By default, \code{tobitnet} standardizes the predictors in \code{x} so that their coefficients are all equally penalized.  
If a predictor has zero variance, we set its coefficient to zero and exclude it from the model fitting procedure.  
}

\subsection{The Tobit model with a weighted elastic net penalty}{
To model the censored response \eqn{y}, we assume that there is a latent uncensored response  \eqn{y^*}{y*} and that it comes from a linear model.  That is, we assume \eqn{y = \max\{y^*,c\}}{y = max{y*,c}} where \eqn{y^* = x' \beta + \epsilon}{y* = x' \beta + \epsilon}, \eqn{\epsilon \sim N(0, \sigma^2)}{\epsilon ~ N(0, \sigma^2)}.  
We use Olsen's (1978) reparameterization of the Tobit log likelihood, with \eqn{\delta_j = \beta_j/\sigma}{\deltaj = \betaj/\sigma} for \eqn{j = 0, 1, \ldots, p} and \eqn{\gamma = 1/\sigma}. We aim to minimize the following objective function:
\deqn{R_n(\delta, \gamma) = - \frac{1}{n} \log L_n(\delta, \gamma) + \sum_{j = 1}^p \lambda_1 w_{1j} |\delta_j| + \lambda_2 w_{2j} \frac{\delta_j^2}{2} }{
Rn(\delta, \gamma) = - 1/n log Ln(\delta, \gamma) + \sum \lambda1 w1j |\deltaj| + \lambda2 w2j \deltaj^2/2
} 
where \eqn{L_n(\delta, \gamma)}{Ln(\delta, \gamma)} denotes the reparameterized Tobit likelihood, the \eqn{w_{1j}}{w1j} are the penalty factors for the lasso penalty, and the \eqn{w_{2j}}{w2j} are the penalty factors for the ridge penalty.  
The weighted elastic net has ridge regression, the lasso, the elastic net, the adaptive lasso, and the adaptive elastic net as special cases. 
}

\subsection{GCD algorithm details}{
We use a generalized coordinate descent (GCD) algorithm to minimize \eqn{R_n(\delta, \gamma)}{Rn(\delta, \gamma)}.  For a fixed \eqn{\lambda_2}{\lambda2}, we fit \eqn{R_n(\delta, \gamma)}{Rn(\delta, \gamma)} for a sequence of \eqn{\lambda_1}{\lambda1} values.  We start at the smallest \eqn{\lambda_1}{\lambda1} value such that all of the coefficient estimates are 0, which we denote by \eqn{\lambda_{\max}}{\lambdamax}.  
We set the smallest value in the \eqn{\lambda_1}{\lambda1} sequence to be \eqn{\lambda_{\min} = \lambda_{\max}}{\lambdamin = \lambdamax}\code{*lambda.factor}.  The remaining \code{nlambda}\eqn{- 2} values in the \eqn{\lambda_1}{\lambda1} sequence are evenly spaced between \eqn{\lambda_{\max}}{\lambdamax} and \eqn{\lambda_{\min}}{\lambdamin} on the natural log scale.
  
We use warm starts to speed up computation of the solution path.  That is, we use the estimated parameters for \eqn{\lambda_1[k]}{\lambda1[k]} as initial values for fitting the model with penalty \eqn{\lambda_1[k+1]}{\lambda1[k+1]} for \eqn{k = 1, \ldots,}\code{length(lambda1)}\eqn{-1}.  
We also incorporate active set cycling into the GCD algorithm to speed up computation.  We start by cycling through all \eqn{p} coefficients, \eqn{\delta_j}{\deltaj}.  Any coefficients which remain at zero after this cycle are removed from the active set.  From there, we only cycle through coefficients in the active set until our convergence criterion is met---that is, until \eqn{\max_j (\hat{\delta}_j^{current} - \hat{\delta}_j^{new})^2 <}{max (\deltaj(current) - \deltaj(new))^2 < }\code{eps}.  
As a final check, we run through all \eqn{p} coefficients again to see if the active set changes.  If it changes, we continue to cycle through updates until our convergence criterion is met again.  If it does not change, the algorithm ends and we check the KKT conditions to verify that the algorithm has converged to a minimizer (if the KKT conditions are \emph{not} satisifed, \code{tobitnet} issues a warning but still returns the solution).
  
We compute the model \dQuote{deviance} which is given by \code{dev[k]}\eqn{= -2\log L_n(\delta, \gamma) }{= -2 log Ln(\delta, \gamma)} for each \eqn{\lambda_1[k]}{\lambda1[k]} (we exclude the log-likelihood for the saturated model from our expression for the deviance for technical reasons).  We stop early in the solution path if \eqn{|(dev[k] - dev[k+1])/nulldev| < 1e-5} (provided that \code{early.stop} is set to \code{TRUE}).
}

\subsection{Parameter estimates}{
\code{tobitnet} returns \eqn{\beta} and \eqn{\sigma} rather than \eqn{\delta} and \eqn{\gamma} to faciliate comparison to parameter estimates from other regression models. 
}

}

\author{
Tate Jacobson
}

\references{
Olsen, R.J. (1978) Note on the Uniqueness of the Maximum Likelihood Estimator for the
Tobit Model. \emph{Econometrica}, \bold{46(5)}, 1211--1215.

Jacobson, T. and Zou, H. (2022) \emph{High-dimensional Censored Regression via the Penalized Tobit Likelihood.} Manuscript submitted for publication.
}

\seealso{\code{\link{predict.tobitnet}}, \code{\link{plot.tobitnet}}, and \code{\link{cv.tobitnet}}}

\examples{
n <- 100
p <- 10
x <- vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
y <- 3 + x\%*\%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y <- pmax(y, 0)

#With a user-provided lambda1
tnet1 <- tobitnet(x = x, y = y, lambda1 = 0.05)

#With lambda1 automatically generated
tnet2 <- tobitnet(x = x, y = y)
}
