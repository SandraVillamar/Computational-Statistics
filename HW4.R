#PROBLEM 1

#standardizing data X
standardizeX <- function(X)
{
  mX = matrix(colMeans(X), nr=dim(X)[1], nc=dim(X)[2], byrow=T)
  X0 = X - mX
  sd = sqrt(colMeans(X0^2))
  X0 = X0%*%diag(1/sd)
  output = list(X0=X0, diag=sd)
  return(output)
}

#subgradient descent for optimizing linear MAE
sub_gd <- function(X,y,eta)
{
  n = dim(X)[1]
  d = dim(X)[2]
  standX = standardizeX(X)
  X0 = cbind(rep(1,n), standX$X0)
  
  b0 = numeric(d+1)
  res0 = y
  b1 = rep(1, d+1)
  iter = 0
  
  repeat{
    if(max(abs(b1-b0)) < 10^(-5) | iter > 20000) break
    
    b1 = b0
    subgrad = sign(res0)
    grad = -t(X0) %*% subgrad/n
    b0 = b0 - eta*grad
    res0 = y - X0%*%b0
    iter = iter + 1
  }
  
  b0[-1] = b0[-1]/standX$diag
  b0[1] = b0[1] - sum(colMeans(X)*b0[-1])
  output = list(beta=b0, res=res0, iter=iter)
  return(output)
  
}

#-----------------------------------------------------------------
#PROBLEM 2

df = read.csv("education.csv")
X = data.matrix(df[,4:6])
y = data.matrix(df[,7])
#subgradient descent for MAE
maefit <- sub_gd(X,y,1)
#MSE
lsefit <- lm(y~X)

#-----------------------------------------------------------------
#PROBLEM 3
#PART 1
build_poly <- function(x,d)
{
  X = matrix(nrow=length(x), ncol=d)
  for(j in 1:d)
  {
    X[,j] = x^(j)
  }
  X = cbind(1, X)
  return(X)
}

#PART 2
df = read.csv("polynomial.csv")
x = data.matrix(df[,1])
y = data.matrix(df[,2])

#lse for d = 1,3,7,12
X1 <- build_poly(x,1)
coef1 <- lm(y~X1)$coefficients
coef1 <- coef1[-2]

X3 <- build_poly(x,3)
coef3 <- lm(y~X3)$coefficients
coef3 <- coef3[-2]

X7 <- build_poly(x,7)
coef7 <- lm(y~X7)$coefficients
coef7 <- coef7[-2]

X12 <- build_poly(x,12)
coef12 <- lm(y~X12)$coefficients
coef12 <- coef12[-2]

plot(x,y)
lines(sort(x), fitted(lm(y~X1))[order(x)], col='blue')
lines(sort(x), fitted(lm(y~X3))[order(x)], col='green')
lines(sort(x), fitted(lm(y~X7))[order(x)], col='purple')
lines(sort(x), fitted(lm(y~X12))[order(x)], col='red')

#PART 3
yhat1 = X1%*%coef1
MSE1 = sum((y-yhat1)^2) / 50
RMSE1 = sqrt(2*MSE1)

yhat3 = X3%*%coef3
MSE3 = sum((y-yhat3)^2) / 50
RMSE3 = sqrt(2*MSE3)

yhat7 = X7%*%coef7
MSE7 = sum((y-yhat7)^2) / 50
RMSE7 = sqrt(2*MSE7)

yhat12 = X12%*%coef12
MSE12 = sum((y-yhat12)^2) / 50
RMSE12 = sqrt(2*MSE12)