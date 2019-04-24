WomenWork = read.table("WomenWork.dat", header = TRUE)
y = WomenWork$Work
X = WomenWork[,2:9]
X = as.matrix(X)
n = nrow(WomenWork) # nr of observations

glmModel <- glm(Work ~ 0 + ., data = WomenWork, family = binomial)
summary(glmModel)


LogPostLogistic <- function(beta, y, X, mu, Sigma){
  m <- dim(X)[2]
  
}