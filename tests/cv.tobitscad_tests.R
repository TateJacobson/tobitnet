#Tests for cv.tobitscad
library(tobitnet)

set.seed(2022)
n = 100
p = 10
x = vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )

y = 3 + x%*%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y = pmax(y, 0)

#Default arguments work
tscv1 = cv.tobitscad(x, y)
class(tscv1) == "cv.tobitscad"
length(tscv1) == 7

is.numeric(tscv1$cvm)
all(is.finite(tscv1$cvm))
length(tscv1$cvm) == 100

is.numeric(tscv1$cvsd)
all(is.finite(tscv1$cvsd))
length(tscv1$cvsd) == 100

is.numeric(tscv1$lambda)
all(is.finite(tscv1$lambda))
all(tscv1$lambda > 0)
length(tscv1$lambda) == 100

is.numeric(tscv1$lambda.min)
is.finite(tscv1$lambda.min)
tscv1$lambda.min > 0
length(tscv1$lambda.min) == 1

is.numeric(tscv1$lambda.1se)
is.finite(tscv1$lambda.1se)
tscv1$lambda.1se > 0
length(tscv1$lambda.1se) == 1

#Works when x is a data frame
x.df = data.frame(x)
tscv1 = cv.tobitscad(x.df, y)
class(tscv1) == "cv.tobitscad"
length(tscv1) == 7

is.numeric(tscv1$cvm)
all(is.finite(tscv1$cvm))
length(tscv1$cvm) == 100

is.numeric(tscv1$cvsd)
all(is.finite(tscv1$cvsd))
length(tscv1$cvsd) == 100

is.numeric(tscv1$lambda)
all(is.finite(tscv1$lambda))
all(tscv1$lambda > 0)
length(tscv1$lambda) == 100

is.numeric(tscv1$lambda.min)
is.finite(tscv1$lambda.min)
tscv1$lambda.min > 0
length(tscv1$lambda.min) == 1

is.numeric(tscv1$lambda.1se)
is.finite(tscv1$lambda.1se)
tscv1$lambda.1se > 0
length(tscv1$lambda.1se) == 1

#Works with user-provided lambda
tscv1 = cv.tobitscad(x,y, lambda = c(0.1, 0.2))
class(tscv1) == "cv.tobitscad"
length(tscv1) == 7

is.numeric(tscv1$cvm)
all(is.finite(tscv1$cvm))
length(tscv1$cvm) == 2

is.numeric(tscv1$cvsd)
all(is.finite(tscv1$cvsd))
length(tscv1$cvsd) == 2

all.equal(tscv1$lambda,c(0.1, 0.2))

is.numeric(tscv1$lambda.min)
is.finite(tscv1$lambda.min)
tscv1$lambda.min > 0
length(tscv1$lambda.min) == 1

is.numeric(tscv1$lambda.1se)
is.finite(tscv1$lambda.1se)
tscv1$lambda.1se > 0
length(tscv1$lambda.1se) == 1

#Works with type.measure = "deviance"
tscv1 = cv.tobitscad(x,y, lambda = c(0.1, 0.2), type.measure = "deviance")
class(tscv1) == "cv.tobitscad"
length(tscv1) == 7

is.numeric(tscv1$cvm)
all(is.finite(tscv1$cvm))
length(tscv1$cvm) == 2

is.numeric(tscv1$cvsd)
all(is.finite(tscv1$cvsd))
length(tscv1$cvsd) == 2

all.equal(tscv1$lambda,c(0.1, 0.2))

is.numeric(tscv1$lambda.min)
is.finite(tscv1$lambda.min)
tscv1$lambda.min > 0
length(tscv1$lambda.min) == 1

is.numeric(tscv1$lambda.1se)
is.finite(tscv1$lambda.1se)
tscv1$lambda.1se > 0
length(tscv1$lambda.1se) == 1

#Works with type.measure = "mae"
tscv1 = cv.tobitscad(x,y, lambda = c(0.1, 0.2), type.measure = "mae")
class(tscv1) == "cv.tobitscad"
length(tscv1) == 7

is.numeric(tscv1$cvm)
all(is.finite(tscv1$cvm))
length(tscv1$cvm) == 2

is.numeric(tscv1$cvsd)
all(is.finite(tscv1$cvsd))
length(tscv1$cvsd) == 2

all.equal(tscv1$lambda,c(0.1, 0.2))

is.numeric(tscv1$lambda.min)
is.finite(tscv1$lambda.min)
tscv1$lambda.min > 0
length(tscv1$lambda.min) == 1

is.numeric(tscv1$lambda.1se)
is.finite(tscv1$lambda.1se)
tscv1$lambda.1se > 0
length(tscv1$lambda.1se) == 1

#Works with different c, nfolds
y = pmax(y,3)

tscv1 = cv.tobitscad(x, y, c = 3, lambda = c(0.1, 0.2), nfolds = 5)
class(tscv1) == "cv.tobitscad"
length(tscv1) == 7

is.numeric(tscv1$cvm)
all(is.finite(tscv1$cvm))
length(tscv1$cvm) == 2

is.numeric(tscv1$cvsd)
all(is.finite(tscv1$cvsd))
length(tscv1$cvsd) == 2

all.equal(tscv1$lambda,c(0.1, 0.2))

is.numeric(tscv1$lambda.min)
is.finite(tscv1$lambda.min)
tscv1$lambda.min > 0
length(tscv1$lambda.min) == 1

is.numeric(tscv1$lambda.1se)
is.finite(tscv1$lambda.1se)
tscv1$lambda.1se > 0
length(tscv1$lambda.1se) == 1

#Works with constant column in some folds
tscv1 = cv.tobitscad(cbind(c(rep(1, n-1) ,0), x), y, c = 3, lambda = c(0.1, 0.2), nfolds = 5)
class(tscv1) == "cv.tobitscad"
length(tscv1) == 7

is.numeric(tscv1$cvm)
all(is.finite(tscv1$cvm))
length(tscv1$cvm) == 2

is.numeric(tscv1$cvsd)
all(is.finite(tscv1$cvsd))
length(tscv1$cvsd) == 2

all.equal(tscv1$lambda,c(0.1, 0.2))

is.numeric(tscv1$lambda.min)
is.finite(tscv1$lambda.min)
tscv1$lambda.min > 0
length(tscv1$lambda.min) == 1

is.numeric(tscv1$lambda.1se)
is.finite(tscv1$lambda.1se)
tscv1$lambda.1se > 0
length(tscv1$lambda.1se) == 1
