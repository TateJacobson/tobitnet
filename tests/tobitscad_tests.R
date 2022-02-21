#Tests for tobitscad
library(tobitnet)

set.seed(2022)
n = 100
p = 10
x = vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
y = 3 + x%*%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y = pmax(y, 0)

#Works with user-provided lambda1
tsout1 <- tobitscad(x = x, y = y, lambda = 0.1)
class(tsout1) == "tobitscad"
length(tsout1) == 8

is.numeric(tsout1$sigma)
is.finite(tsout1$sigma)
all(tsout1$sigma > 0)
length(tsout1$sigma) == 1

is.numeric(tsout1$b0)
is.finite(tsout1$b0)
length(tsout1$b0) == 1

is.numeric(tsout1$beta)
all(is.finite(tsout1$beta))
ncol(tsout1$beta) == 1
nrow(tsout1$beta) == p

tsout1$lambda == 0.1

tsout1$c == 0

#Works when x is a data frame
x.df = data.frame(x)
tsout1 <- tobitscad(x = x.df, y = y, lambda = 0.1)
class(tsout1) == "tobitscad"
length(tsout1) == 8

is.numeric(tsout1$sigma)
is.finite(tsout1$sigma)
all(tsout1$sigma > 0)
length(tsout1$sigma) == 1

is.numeric(tsout1$b0)
is.finite(tsout1$b0)
length(tsout1$b0) == 1

is.numeric(tsout1$beta)
all(is.finite(tsout1$beta))
ncol(tsout1$beta) == 1
nrow(tsout1$beta) == p

tsout1$lambda == 0.1

tsout1$c == 0

#Works with different c, lambda, a, iter
tsout1 <- tobitscad(x = x, y = pmax(y,5), c = 5, a = 3.7, iter = 4, lambda = c(0.2, 0.3) )
class(tsout1) == "tobitscad"
length(tsout1) == 8

is.numeric(tsout1$sigma)
all(is.finite(tsout1$sigma))
all(tsout1$sigma > 0)
length(tsout1$sigma) == 2

is.numeric(tsout1$b0)
all(is.finite(tsout1$b0))
length(tsout1$b0) == 2

is.numeric(tsout1$beta)
all(is.finite(tsout1$beta))
ncol(tsout1$beta) == 2
nrow(tsout1$beta) == p

all.equal(tsout1$lambda, c(0.2, 0.3))

tsout1$c == 5

#Default arguments work
tsout1 <- tobitscad(x = x, y = y)
is.numeric(tsout1$sigma)
all(is.finite(tsout1$sigma))
all(tsout1$sigma > 0)
length(tsout1$sigma) == 100

is.numeric(tsout1$b0)
all(is.finite(tsout1$b0))
length(tsout1$b0) == 100

is.numeric(tsout1$beta)
all(is.finite(tsout1$beta))
ncol(tsout1$beta) == 100
nrow(tsout1$beta) == p

length(tsout1$lambda) == 100
is.numeric(tsout1$lambda)
all(is.finite(tsout1$lambda))

tsout1$c == 0

#Works with standardize
tsout1sF <- tobitscad(x = x, y = y, lambda = 0.2, standardize = F)

#Intercept column does not impact model fitting
tsout1b <- tobitscad(x = cbind(rep(1,n) ,x), y = y)
all.equal(tsout1$sigma, tsout1b$sigma)
all.equal(tsout1$b0, tsout1b$b0)
all.equal(tsout1$beta, tsout1b$beta[-1,])
all.equal( nrow(tsout1$beta) +1, nrow(tsout1b$beta))
all.equal( ncol(tsout1$beta), ncol(tsout1b$beta))

#Works with a different beta/sigma combination
y = x%*%c(3, rep(0,4), 3, rep(0,4)) + rnorm(n,0,10)

y = pmax(y, 0)
tsout1 <- tobitscad(x = x, y = y)
is.numeric(tsout1$sigma)
all(is.finite(tsout1$sigma))
all(tsout1$sigma > 0)
length(tsout1$sigma) == 100

is.numeric(tsout1$b0)
all(is.finite(tsout1$b0))
length(tsout1$b0) == 100

is.numeric(tsout1$beta)
all(is.finite(tsout1$beta))
ncol(tsout1$beta) == 100
nrow(tsout1$beta) == p

length(tsout1$lambda) == 100
is.numeric(tsout1$lambda)
all(is.finite(tsout1$lambda))

tsout1$c == 0
