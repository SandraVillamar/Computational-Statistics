df <- read.csv(file="hitters.csv", header=TRUE, row.names=1)
dataMatrix <- data.matrix(df)
y <- dataMatrix[,19] #salary col
X <- dataMatrix[,-19] #263x19 matrix of predictors

#PROBLEM 1: Part 1--------------------------------------------------------------------------

#M = matrix of "best" variables in order M1 to M19
M <- matrix(0, nr=263, nc=19)

#store min RSS of each model M0 to M19
rssFinal <- numeric(20)

#M0 = sample mean of y
M0 <- mean(y)
#RSS for model M0
rssFinal[1] <- ( norm((y-M0),type="2") )^2

#col 1-19 = M1 to M19 = col's of X in order of best fit
for(k in 1:19)
{
  #create vector for RSS of each possible col to add to model
  rss <- numeric(20-k)
  for(j in 1:(20-k))
  {
    #calculate RSS for the possible models with k variables
    if(k == 1) {
      possibleModel <- as.matrix(X[,j])
    }
    else if((20-k) > 1) { 
      possibleModel <- cbind(M[,1:k-1], as.matrix(X[,j]))
    }
    else {
      #adding last col left in X since k = 19
      X <- as.matrix(X)
      possibleModel <- cbind(M[,1:k-1], X)
    }
    
    rss[j] = (norm( as.numeric((lm(y~possibleModel))$residuals), type="2" ))^2 
  }
  #store min rss for model with k variables
  rssFinal[k+1] <- min(rss) 
  #add "best" col to M
  M[,k] = X[,which.min(rss)] 
  #remove "best" col from matrix X
  if(ncol(X) > 1) {
    X <- X[,-which.min(rss)] 
  }
}

#calculate BIC for all models
bic <- 263*log(rssFinal/263) + log(263)*c(0:19)
bestModelBIC <- min(bic)
#num of variables in the best model
bestModelNum <- which.min(bic) - 1
#beta of best model
beta <- (lm(y~M[,1:bestModelNum]))$coefficients

#plot BIC vs. # of variables
plot(c(0:19), bic, main="Forward Step Selection", xlab="# of Variables", ylab="BIC")
lines(c(0:19), bic, type="l")


#PROBLEM 1: Part 2--------------------------------------------------------------------------
#reset X
X <- dataMatrix[,-19] #263x19 matrix of predictors
y <- dataMatrix[,19] #salary col

#matrix of order of col removal from 19 variables in X
MBack <- matrix(0, nr=263, nc=18)

#store min RSS of each model Mk
rssFinalBack <- numeric(20)

#M0 = sample mean of y
M0 <- mean(y)
#RSS for model M0
rssFinalBack[1] <- ( norm((y-M0), type="2") )^2
#RSS for model M19
rssFinalBack[20] <- ( norm(as.numeric((lm(y~X))$residuals), type="2") )^2

for(k in 0:17)
{
  #create vector for RSS of each possible col to remove from model
  rssBack <- numeric(19-k)
  for(j in 1:(19-k))
  {
    #calculate RSS for each model with 18-k variables
    possibleModel <- as.matrix(X[,-j])
    rssBack[j] = (norm( as.numeric((lm(y~possibleModel))$residuals), type="2" ))^2
  }
  #store min rss for model with 18-k variables
  rssFinalBack[19-k] <- min(rssBack) 
  #add "best" col to remove from X
  MBack[,k+1] = X[,which.min(rssBack)] 
  #remove col from matrix X leaving cols with "best" RSS
  X <- X[,-which.min(rssBack)] 
}
#add last col remaining from X
MBack <- cbind(MBack, matrix(X))

#calculate BIC for all models
bicBack <- 263*log(rssFinalBack/263) + log(263)*c(0:19)
bestModelBack <- min(bicBack)
#num of variables in the best model
bestModelNumBack <- 19 - which.min(bicBack) + 2
#beta of best model
beta <- (lm(y~M[,8:19]))$coefficients

#plot BIC vs. # of variables
plot(c(20:1), bicBack, main="Backward Step Selection", xlab="# of Variables", ylab="BIC")
lines(c(20:1), bicBack, type="l")


#PROBLEM 2: Part 1--------------------------------------------------------------------------
#reset X and y
df <- read.csv(file="hitters.csv", header=TRUE, row.names=1)
dataMatrix <- data.matrix(df)
X <- dataMatrix[,-19] #263x19 matrix of predictors
y <- matrix(dataMatrix[,19]) #salary col

#choose grid of lambda values: 
grid <- 10^seq(10,-2,length=100)

#standardizing data X
standardizeX <- function(X)
{
  mX <- matrix(colMeans(X), nr=dim(X)[1], nc=dim(X)[2], byrow=T)
  X0 <- X - mX  #centered
  sd <- sqrt(colMeans(X0^2))
  X0 <- X0%*%diag(1/sd)  #standardized
  output <- list(X0=X0, diag=sd)
  return(output)
}

#create 20x100 matrix of ridge coef for each lambda
makeRidgeCoefMatrix <- function(X, y, grid) {
  #stand X and center y
  X0 <- standardizeX(X)$X0
  sd <- standardizeX(X)$diag
  ycent <- y - mean(y)
  
  #calculate omega with grid of lambdas
  omega <- matrix(0, nr=19, nc=100)
  for(i in 1:100) {
    omega[,i] <- solve(t(X0)%*%X0 + grid[i]*diag(19)) %*% t(X0)%*%ycent
  }
  
  #recover to beta with original X and y
  beta <- matrix(0, nr=20, nc=100)
  for(i in 1:19) {
    beta[i+1,] <- omega[i,] / sd[i]
  }
  beta[1,] <- mean(y) - t(colMeans(X)) %*% beta[2:20,]
  
  return(beta)
}

#create matrix
coefMatrix <- makeRidgeCoefMatrix(X,y,grid)


#PROBLEM 2: Part 2--------------------------------------------------------------------------

#choose grid of lambda values: 
grid <- 10^seq(10,-2,length=100)

#randomly choose 10-fold
set.seed(1)
perm <- sample(263, 263)
#matrix s.t. row i = list of indices of rows in X that are in fold i
fold <- matrix(0, nr=10, nc=27)
for(k in 1:10) {
  #first 7 folds have 26 obs, so last entry is 0
  if(k <= 7) {
    fold[k,] <- c( perm[ (1+26*(k-1)):(26*k) ], 0)
  }
  #last 3 folds have 27 obs
  else {
    fold[k,] <- perm[(183+27*(k-8)):(209+27*(k-8))]
  }
}

#cross validation for each lambda
crossv <- numeric(100) 
#for lambda(1) to lambda(100)
for(i in 1:100) {
  #prediction error for each fold
  pe <- numeric(10)
  #from fold 1 to 10
  for(k in 1:10) {
    #create sd vector
    sd <- numeric(19)
    #assign Xt = training set: X without fold k
    #Xv = validation set
    if(k <= 7) {
      #get X0 and ycent for each fold
      X0 <- standardizeX(X[-fold[k,-27], ])$X0
      sd <- standardizeX(X[-fold[k,-27], ])$diag
      ycent <- y[-fold[k,-27]] - mean(y[-fold[k,-27]])
      #fold 1-7 have 26 obs, so remove last element in fold[i,]
      Xt <- X0
      yt <- ycent
      Xv <- X[fold[k,-27], ]
      yv <- y[fold[k,-27]]
    }
    else {
      #get X0 and ycent for each fold
      X0 <- standardizeX(X[-fold[k,], ])$X0
      sd <- standardizeX(X[-fold[k,], ])$diag
      ycent <- y[-fold[k,]] - mean(y[-fold[k,]])
      #fold 8-10 have 27 obs
      Xt <- X0
      yt <- ycent
      Xv <- X[fold[k,], ]
      yv <- y[fold[k,]]
    }
    
    #perform ridge reg on Xt with lambda(i)
    omegaRidge = solve(t(Xt)%*%Xt + grid[i]*diag(19)) %*% t(Xt)%*%yt
    
    #recover to beta with original X and y
    betaRidge <- numeric(20)
    for(j in 1:19) {
      betaRidge[j+1] <- omegaRidge[j] / sd[j]
    }
    betaRidge[1] <- mean(yt) - t(colMeans(Xt)) %*% betaRidge[2:20]
    
    #prediction error of betaRidge on Xv
    pe[k] = mean((yv - betaRidge[1] - Xv%*%betaRidge[2:20])^2)
  }
  #cross valid for lambda(i)
  crossv[i] <- mean(pe)
}

#plot CV versus lambda
plot(log(grid), crossv, main="CV Error vs. Lambda", xlab="log(Lambda)", ylab="CV Error")
lines(log(grid), crossv, type="l")
lambdaBest <- grid[which.min(crossv)]
lambdaBestInd <- which.min(crossv)


#PROBLEM 2: Part 3--------------------------------------------------------------------------

#get ridge reg coef with lambdaBest from Part 2 on full data
betaMatrix <- makeRidgeCoefMatrix(X,y,grid)
betaBest <- betaMatrix[,lambdaBestInd]


#PROBLEM 3: Part 1--------------------------------------------------------------------------

#package for which.maxn function
install.packages("doBy")
library(doBy)

#graHTP algorithm: ycent = centered y, X0 = stand X
grahtp <- function(ycent, X0, k, step) {
  
  #row and col of X0
  n = dim(X0)[1]
  p = dim(X0)[2]
  #initialize variable beta_0 = 0
  bt <- numeric(p)
  
  repeat {
    #initialize beta_t+1
    bt1 <- numeric(p)
    #calculate gradient of bt
    grad <- -(1/n)*t(X0)%*%(ycent-X0%*%bt)
    #get bt1
    bt1 <- bt - step*grad
    #get indices of largest k elements in bt1
    k_ind <- which.maxn(abs(bt1), k)
    #submatrix with k_ind cols
    Xk <- X0[,k_ind]
    #calculate bt1 with sparsity k
    bt1[k_ind] <- solve(t(Xk)%*%Xk) %*% t(Xk)%*%ycent
    bt1[-k_ind] <- 0
    
    #stopping condition
    if(norm(bt1-bt, type="2") < 10^-4) break
    
    #otherwise continue iteration
    bt <- bt1
  }
  
  #output
  return(bt1)
}


#PROBLEM 3: Part 2--------------------------------------------------------------------------

#reset X and y
df <- read.csv(file="hitters.csv", header=TRUE, row.names=1)
dataMatrix <- data.matrix(df)
X <- dataMatrix[,-19] #263x19 matrix of predictors
y <- matrix(dataMatrix[,19]) #salary col

#standardize X, center y
X0 <- standardizeX(X)$X0
sd <- standardizeX(X)$diag
ycent <- mean(y) - y

#matrix where col(k) = best beta for sparsity k
bestBetas <- matrix(0, nr=20, nc=19)
#store rss for each model Mk
rss <- numeric(19)

#get best omegas, recover to beta
for(k in 1:19) {
  bestOmega <- grahtp(ycent, X0, k, step=.1)
  
  for(j in 1:19) {
    bestBetas[j+1,k] <- bestOmega[j] / sd[j]
  }
  bestBetas[1,k] <- mean(y) - t(colMeans(X)) %*% bestBetas[2:20,k]
  
  #calculate rss for each model
  rss[k] <- sum((y - bestBetas[1,k] - X%*%bestBetas[2:20,k])^2)
}

#calculate BIC for all models
bic <- 263*log(rss/263) + log(263)*c(1:19)
bestModelBIC <- min(bic)
#num of variables in the best model
bestModelNum <- which.min(bic)

#best 1 var Model
beta_1var <- bestBetas[,1]
