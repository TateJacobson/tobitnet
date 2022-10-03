#Tests for cv.tobitnet
library(tobitnet)

set.seed(2022)
n = 100
p = 10
x = vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )

y = 3 + x%*%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y = pmax(y, 0)

#Default arguments work
tncv1 = cv.tobitnet(x,y)
class(tncv1) == "cv.tobitnet"
length(tncv1) == 8

is.numeric(tncv1$cvm)
all(is.finite(tncv1$cvm))
length(tncv1$cvm) == 100

is.numeric(tncv1$cvsd)
all(is.finite(tncv1$cvsd))
length(tncv1$cvsd) == 100

is.numeric(tncv1$lambda1)
all(is.finite(tncv1$lambda1))
all(tncv1$lambda1 > 0)
length(tncv1$lambda1) == 100

tncv1$lambda2 == 0

is.numeric(tncv1$lambda1.min)
is.finite(tncv1$lambda1.min)
tncv1$lambda1.min > 0
length(tncv1$lambda1.min) == 1

is.numeric(tncv1$lambda1.1se)
is.finite(tncv1$lambda1.1se)
tncv1$lambda1.1se > 0
length(tncv1$lambda1.1se) == 1

#Works when x is a data frame
x.df = data.frame(x)
tncv1 = cv.tobitnet(x.df, y)
class(tncv1) == "cv.tobitnet"
length(tncv1) == 8

is.numeric(tncv1$cvm)
all(is.finite(tncv1$cvm))
length(tncv1$cvm) == 100

is.numeric(tncv1$cvsd)
all(is.finite(tncv1$cvsd))
length(tncv1$cvsd) == 100

is.numeric(tncv1$lambda1)
all(is.finite(tncv1$lambda1))
all(tncv1$lambda1 > 0)
length(tncv1$lambda1) == 100

tncv1$lambda2 == 0

is.numeric(tncv1$lambda1.min)
is.finite(tncv1$lambda1.min)
tncv1$lambda1.min > 0
length(tncv1$lambda1.min) == 1

is.numeric(tncv1$lambda1.1se)
is.finite(tncv1$lambda1.1se)
tncv1$lambda1.1se > 0
length(tncv1$lambda1.1se) == 1

#Works with user-provided lambda1, lambda2
tncv1 = cv.tobitnet(x,y, lambda1 = c(0.1, 0.2), lambda2 = 0.1)
class(tncv1) == "cv.tobitnet"
length(tncv1) == 8

is.numeric(tncv1$cvm)
all(is.finite(tncv1$cvm))
length(tncv1$cvm) == 2

is.numeric(tncv1$cvsd)
all(is.finite(tncv1$cvsd))
length(tncv1$cvsd) == 2

all.equal(tncv1$lambda1,c(0.1, 0.2))
tncv1$lambda2 == 0.1

is.numeric(tncv1$lambda1.min)
is.finite(tncv1$lambda1.min)
tncv1$lambda1.min > 0
length(tncv1$lambda1.min) == 1

is.numeric(tncv1$lambda1.1se)
is.finite(tncv1$lambda1.1se)
tncv1$lambda1.1se > 0
length(tncv1$lambda1.1se) == 1

#Works with type.measure = "deviance"
tncv1 = cv.tobitnet(x, y, lambda1 = c(0.1, 0.2), lambda2 = 0.1, type.measure = "deviance")
class(tncv1) == "cv.tobitnet"
length(tncv1) == 8

is.numeric(tncv1$cvm)
all(is.finite(tncv1$cvm))
length(tncv1$cvm) == 2

is.numeric(tncv1$cvsd)
all(is.finite(tncv1$cvsd))
length(tncv1$cvsd) == 2

all.equal(tncv1$lambda1,c(0.1, 0.2))
tncv1$lambda2 == 0.1

is.numeric(tncv1$lambda1.min)
is.finite(tncv1$lambda1.min)
tncv1$lambda1.min > 0
length(tncv1$lambda1.min) == 1

is.numeric(tncv1$lambda1.1se)
is.finite(tncv1$lambda1.1se)
tncv1$lambda1.1se > 0
length(tncv1$lambda1.1se) == 1

#Works with type.measure = "mae
tncv1 = cv.tobitnet(x, y, lambda1 = c(0.1, 0.2), lambda2 = 0.1, type.measure = "mae")
class(tncv1) == "cv.tobitnet"
length(tncv1) == 8

is.numeric(tncv1$cvm)
all(is.finite(tncv1$cvm))
length(tncv1$cvm) == 2

is.numeric(tncv1$cvsd)
all(is.finite(tncv1$cvsd))
length(tncv1$cvsd) == 2

all.equal(tncv1$lambda1,c(0.1, 0.2))
tncv1$lambda2 == 0.1

is.numeric(tncv1$lambda1.min)
is.finite(tncv1$lambda1.min)
tncv1$lambda1.min > 0
length(tncv1$lambda1.min) == 1

is.numeric(tncv1$lambda1.1se)
is.finite(tncv1$lambda1.1se)
tncv1$lambda1.1se > 0
length(tncv1$lambda1.1se) == 1

#Works with different c, nfolds
y = pmax(y,3)

tncv1 = cv.tobitnet(x, y, c = 3, lambda1 = c(0.1, 0.2), lambda2 = 0.1, nfolds = 5)
class(tncv1) == "cv.tobitnet"
length(tncv1) == 8

is.numeric(tncv1$cvm)
all(is.finite(tncv1$cvm))
length(tncv1$cvm) == 2

is.numeric(tncv1$cvsd)
all(is.finite(tncv1$cvsd))
length(tncv1$cvsd) == 2

all.equal(tncv1$lambda1,c(0.1, 0.2))
tncv1$lambda2 == 0.1

is.numeric(tncv1$lambda1.min)
is.finite(tncv1$lambda1.min)
tncv1$lambda1.min > 0
length(tncv1$lambda1.min) == 1

is.numeric(tncv1$lambda1.1se)
is.finite(tncv1$lambda1.1se)
tncv1$lambda1.1se > 0
length(tncv1$lambda1.1se) == 1

#Works with constant column in some folds
tncv1 = cv.tobitnet(cbind(c(rep(1, n-1) ,0), x) , y, c = 3, lambda1 = c(0.1, 0.2), lambda2 = 0.1, nfolds = 5)
class(tncv1) == "cv.tobitnet"
length(tncv1) == 8

is.numeric(tncv1$cvm)
all(is.finite(tncv1$cvm))
length(tncv1$cvm) == 2

is.numeric(tncv1$cvsd)
all(is.finite(tncv1$cvsd))
length(tncv1$cvsd) == 2

all.equal(tncv1$lambda1,c(0.1, 0.2))
tncv1$lambda2 == 0.1

is.numeric(tncv1$lambda1.min)
is.finite(tncv1$lambda1.min)
tncv1$lambda1.min > 0
length(tncv1$lambda1.min) == 1

is.numeric(tncv1$lambda1.1se)
is.finite(tncv1$lambda1.1se)
tncv1$lambda1.1se > 0
length(tncv1$lambda1.1se) == 1
