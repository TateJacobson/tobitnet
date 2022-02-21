#Tests for predict.tobitscad
library(tobitnet)

set.seed(2022)
n = 100
p = 10
x = vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
x.test = vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )

y = -5 + x%*%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y = pmax(y, 1)

tsout1 = tobitscad(x = x, y = y, c = 1)

#Default arguments work
tsout_pred1 = predict(tsout1, newx = x.test)

#Output is as expected
is.numeric(tsout_pred1)
all( is.finite(tsout_pred1) )
nrow(tsout_pred1) == nrow(x.test)
ncol(tsout_pred1) == length(tsout1$lambda)
#c from tobitnet passes to predict.tobitnet
min(tsout_pred1) == 1

#Accepts a new lambda1
tsout_pred2 = predict(tsout1, newx = x.test, lambda = c(0.1,0.5))
is.numeric(tsout_pred2)
nrow(tsout_pred2) == nrow(x.test)
ncol(tsout_pred2) == 2

#type = "uncensored" works
tsout_pred2 = predict(tsout1, newx = x.test, lambda = c(0.1,0.5), type = "uncensored")
is.numeric(tsout_pred2)
nrow(tsout_pred2) == nrow(x.test)
ncol(tsout_pred2) == 2
#Predicted values can fall below c = 1
min(tsout_pred2) < 1

#Including an intercept column does not impact results
tsout1b = tobitscad(x = cbind(rep(1,n) ,x), y = y, c = 1)

tsout_pred1 = predict(tsout1, newx = x)
tsout_pred1b = predict(tsout1b, newx = cbind(rep(1,n) ,x))
all.equal(tsout_pred1, tsout_pred1b)
all.equal(dim(tsout_pred1), c(n, 100))
all.equal(dim(tsout_pred1), dim(tsout_pred1b))

#Accepts a vector for newx
x.vec = rnorm(10, 1, 1.5)
tsout_pred2 = predict(tsout1, newx = x.vec)
is.numeric(tsout_pred2)
nrow(tsout_pred2) == 1
ncol(tsout_pred2) == length(tsout1$lambda)
