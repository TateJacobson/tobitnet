#Tests for predict.tobitnet
library(tobitnet)

set.seed(2022)
n = 100
p = 10
x = vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )
x.test = vapply(1:p, function(i) rnorm(n, mean = 1, sd = 1.5), FUN.VALUE = numeric(n) )

y = -5 + x%*%c(5, 1, 2, 0.5, 0.1, rep(0, p - 5)) + rnorm(n,0,1)
y = pmax(y, 1)

tnout1 = tobitnet(x = x, y = y, c = 1)

#Default arguments work
tnout_pred1 = predict(tnout1, newx = x.test)

#Output is as expected
is.numeric(tnout_pred1)
all( is.finite(tnout_pred1) )
nrow(tnout_pred1) == nrow(x.test)
ncol(tnout_pred1) == length(tnout1$lambda1)
#c from tobitnet passes to predict.tobitnet
min(tnout_pred1) == 1

#Accepts a new lambda1
tnout_pred2 = predict(tnout1, newx = x.test, lambda1 = c(0.1,0.5))
is.numeric(tnout_pred2)
nrow(tnout_pred2) == nrow(x.test)
ncol(tnout_pred2) == 2

#type = "uncensored" works
tnout_pred2 = predict(tnout1, newx = x.test, lambda1 = c(0.1,0.5), type = "uncensored")
is.numeric(tnout_pred2)
nrow(tnout_pred2) == nrow(x.test)
ncol(tnout_pred2) == 2
#Predicted values can fall below c = 1
min(tnout_pred2) < 1

#Including an intercept column does not impact results
tnout1b = tobitnet(x = cbind(rep(1,n) ,x), y = y, c = 1)

tnout_pred1 = predict(tnout1, newx = x)
tnout_pred1b = predict(tnout1b, newx = cbind(rep(1,n) ,x))
all.equal(tnout_pred1, tnout_pred1b)
all.equal(dim(tnout_pred1), c(n, 100))
all.equal(dim(tnout_pred1), dim(tnout_pred1b))

#Accepts a vector for newx
x.vec = rnorm(10, 1, 1.5)
tnout_pred2 = predict(tnout1, newx = x.vec)
is.numeric(tnout_pred2)
nrow(tnout_pred2) == 1
ncol(tnout_pred2) == length(tnout1$lambda1)
