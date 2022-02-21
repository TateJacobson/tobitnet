#Tests for the tobitnet function
library(AER)
library(tobitnet)

set.seed(2022)
n = 100
p = 10
x = vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
y = 3 + x%*%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y = pmax(y, 0)

#Works with user-provided lambda1
tnout1 <- tobitnet(x = x, y = y, lambda1 = 0.05)
class(tnout1) == "tobitnet"
length(tnout1) == 9

is.numeric(tnout1$sigma)
is.finite(tnout1$sigma)
all(tnout1$sigma > 0)
length(tnout1$sigma) == 1

is.numeric(tnout1$b0)
is.finite(tnout1$b0)
length(tnout1$b0) == 1

is.numeric(tnout1$beta)
all(is.finite(tnout1$beta))
ncol(tnout1$beta) == 1
nrow(tnout1$beta) == p

tnout1$lambda1 == 0.05
tnout1$lambda2 == 0

tnout1$c == 0

#Works when x is a data frame
x.df = data.frame(x)
tnout1 <- tobitnet(x = x.df, y = y, lambda1 = 0.05)
class(tnout1) == "tobitnet"
length(tnout1) == 9

is.numeric(tnout1$sigma)
is.finite(tnout1$sigma)
all(tnout1$sigma > 0)
length(tnout1$sigma) == 1

is.numeric(tnout1$b0)
is.finite(tnout1$b0)
length(tnout1$b0) == 1

is.numeric(tnout1$beta)
all(is.finite(tnout1$beta))
ncol(tnout1$beta) == 1
nrow(tnout1$beta) == p

tnout1$lambda1 == 0.05
tnout1$lambda2 == 0

tnout1$c == 0

#Works with different c, lambda1, lambda2
tnout1 <- tobitnet(x = x, y = pmax(y,5), c = 5, lambda1 = 0.1, lambda2 = 0.1)
class(tnout1) == "tobitnet"
length(tnout1) == 9

is.numeric(tnout1$sigma)
is.finite(tnout1$sigma)
tnout1$sigma > 0
length(tnout1$sigma) == 1

is.numeric(tnout1$b0)
is.finite(tnout1$b0)
length(tnout1$b0) == 1

is.numeric(tnout1$beta)
all(is.finite(tnout1$beta))
ncol(tnout1$beta) == 1
nrow(tnout1$beta) == p

tnout1$lambda1 == 0.1
tnout1$lambda2 == 0.1

tnout1$c == 5

#Default arguments work
tnout1 <- tobitnet(x = x, y = y)
is.numeric(tnout1$sigma)
all(is.finite(tnout1$sigma))
all(tnout1$sigma > 0)
length(tnout1$sigma) == 100

is.numeric(tnout1$b0)
all(is.finite(tnout1$b0))
length(tnout1$b0) == 100

is.numeric(tnout1$beta)
all(is.finite(tnout1$beta))
ncol(tnout1$beta) == 100
nrow(tnout1$beta) == p

length(tnout1$lambda1) == 100
is.numeric(tnout1$lambda1)
all(is.finite(tnout1$lambda1))

tnout1$lambda2 == 0

tnout1$c == 0

#Standardize = F does not change estimates when lambda1 = lambda2 = 0
tnout1sT <- tobitnet(x = x, y = y, lambda1 = 0, lambda2 = 0, standardize = T, eps = 1e-20)
tnout1sF <- tobitnet(x = x, y = y, lambda1 = 0, lambda2 = 0, standardize = F, eps = 1e-20)
all.equal(tnout1sT$sigma, tnout1sF$sigma)
all.equal(tnout1sT$b0, tnout1sF$b0)
all.equal(tnout1sT$beta, tnout1sF$beta)

#Intercept column does not impact model fitting
tnout1b <- tobitnet(x = cbind(rep(1,n) ,x), y = y, early.stop = F)
all.equal(tnout1$sigma, tnout1b$sigma)
all.equal(tnout1$b0, tnout1b$b0)
all.equal(tnout1$beta, tnout1b$beta[-1,])
all.equal( nrow(tnout1$beta) +1, nrow(tnout1b$beta))
all.equal( ncol(tnout1$beta), ncol(tnout1b$beta))

#Checking tobitnet estimates against tobit estimates
tnout1 <- tobitnet(x = x, y = y, lambda1 = 0, eps = 1e-20)
tob1 <- tobit(y ~ x, left = 0)
all.equal(drop(tnout1$beta), as.vector(tob1$coefficients)[-1])
all.equal(tnout1$b0, as.vector(tob1$coefficients)[1])
all.equal(tob1$scale, tnout1$sigma)

#Check that standardize = F doesn't throw off estimates
tnout1 <- tobitnet(x = x, y = y, c = 0, lambda1 = 0, lambda2 = 0, eps = 1e-20, standardize = F)
all.equal(drop(tnout1$beta), as.vector(tob1$coefficients)[-1])
all.equal(tnout1$b0, as.vector(tob1$coefficients)[1])
all.equal(tob1$scale, tnout1$sigma)

#Check that c > 0 doesn't throw off estimates
y = pmax(y, 5)
tnout1 <- tobitnet(x = x, y = y, c = 5, lambda1 = 0, lambda2 = 0, eps = 1e-20)
tob1 <- tobit(y ~ x, left = 5)
all.equal(drop(tnout1$beta), as.vector(tob1$coefficients)[-1])
all.equal(tnout1$b0, as.vector(tob1$coefficients)[1])
all.equal(tob1$scale, tnout1$sigma)

#Compare with a different beta/sigma combination
y = x%*%c(3, rep(0,4), 3, rep(0,4)) + rnorm(n,0,10)

y = pmax(y, 0)
tnout1 <- tobitnet(x = x, y = y, c = 0, lambda1 = 0, lambda2 = 0, eps = 1e-20)
tob1 <- tobit(y ~ x, left = 0)
all.equal(drop(tnout1$beta), as.vector(tob1$coefficients)[-1])
all.equal(tnout1$b0, as.vector(tob1$coefficients)[1])
all.equal(tob1$scale, tnout1$sigma)
