
rm(list=ls())
library("mvtnorm")
#2.1############################################################################
WomenWork <- read.table("WomenWork.dat", header = TRUE)
glmModel <- glm(Work ~ 0 + ., data = WomenWork, family = binomial)
summary(glmModel)

#2.2############################################################################
m     <- dim(WomenWork)[2]
X     <- as.matrix(WomenWork[,2:m])
y     <- as.matrix(WomenWork[,1])
nPara <- dim(X)[2]

tau   <- 10
mu    <- as.vector(rep(0,nPara))
Sigma <- tau^2*diag(nPara)

LogPostLogistic <- function(betaVect,y,X,mu,Sigma){
  
  nPara <- length(betaVect);
  linPred <- X%*%betaVect;
  
  # evaluating the log-likelihood                                    
  logLik <- sum( linPred*y -log(1 + exp(linPred)));
  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  
  # evaluating the prior
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);
  
  # add the log prior and log-likelihood together to get log posterior
  return(logLik + logPrior)
}

initVal <- as.vector(rnorm(dim(X)[2]))
logPost <- LogPostLogistic
OptimResults<-optim(initVal,logPost,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)
OptimResults$par





