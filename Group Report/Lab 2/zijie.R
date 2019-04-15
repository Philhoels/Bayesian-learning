#1.1############################################################################
rm(list=ls())
library(ggplot2)
library(gridExtra)
set.seed(123456)
data = read.table("TempLinkoping.txt",header = TRUE)

Y = data$temp
X = cbind(1, data$time, data$time**2)
n = nrow(X)

plot(data$time, data$temp,ylim = c(-40,40))

# hyperparameters
      # mean of beta, don't need change
mu0 = c(-10,100,-100)  
      # variance
v0 = 10 # 10 instead of 4
sigma20 = 0.2 # 0.2 instead of 1
omg0 = 0.01*diag(3) 

# joint conjugate prior
prior <- data.frame(matrix(ncol=3,nrow=100))
for (i in 1:100) {
  var = LaplacesDemon::rinvchisq(1,v0,sigma20)
  beta = MASS::mvrnorm(1, mu0, var*solve(omg0))
  prior[i,1:3] = beta
  lines(x=data$time,y=beta%*%t(X),col="pink")
}

###### Zijie idea:
# goal: use all the pink lines cover all points
# so I smaller the variance of beta, which makes the
# pink lines have the same tendancy.
######

#1.2############################################################################

# posterior
k = 3 #  number of regression coefficients
beta_hat = solve(t(X)%*%X)%*%t(X)%*%Y
mun = solve(t(X)%*%X+omg0)%*%(t(X)%*%X%*%beta_hat+omg0%*%mu0)
omgn = t(X)%*%X+omg0
vn = v0+n
sigma2n =( v0*sigma20+(t(Y)%*%Y+t(mu0)%*%omg0%*%mu0-t(mun)%*%omgn%*%mun ) )/vn

# marginal posterior
###### Zijie idea:
# The arguments of marginal posterior of conjugate prior 
# are guessed form marginal posterior of uniform prior (page 6)
######
data_hist = as.data.frame(mvtnorm::rmvt(n=100, delta=mun, sigma=as.numeric(sigma2n)*solve(omgn),df=n-k))
colnames <- c("beta0","beta1","beta2")
colnames(data_hist) <- colnames

f <- function(colname){
  ggplot(data_hist, aes_string(x=colname))+
    geom_histogram(aes(y=..density..),
                   colour="black",
                   fill="white",
                   bins = 30)+
    geom_density(alpha=.2, fill="#FF6666")
}
plot(arrangeGrob(grobs=lapply(colnames, f)))

# overlay a curve for the posterior median
beta_median = matrixStats::colMedians(as.matrix(data_hist))
pred12 = beta_median%*%t(X)

# credible intervals
    # each column represents all possible values at a fixed time
preds <- as.matrix(data_hist)%*%t(X)
pred_interval <- data.frame(nrow=ncol(preds),nrow=2)
colnames(pred_interval) <- c("lower","upper")
for(i in 1:ncol(preds)){
  data_t <- preds[,i]
  pred_interval[i,] <- quantile(data_t,probs = c(0.025,0.975))
}
data12 <- cbind(data,t(pred12),pred_interval)
ggplot(data12)+
  geom_point(aes(x=time, y=temp))+
  geom_line(aes(x=time, y=t(pred12)),color="Blue",size=1)+
  geom_line(aes(x=time, y=lower),color="Red",linetype="dashed",size=1)+
  geom_line(aes(x=time, y=upper),color="Red",linetype="dashed",size=1)

###### Zijie idea:
# NO, the cridible intervals do not contain most of points. They shouldn't.
# Credible intervals describe the robustness of beta, not of model.
# There exists variance when we predict Y using our regression model, 
# and the corresponding intervals describing the condience of model is 
# predective interval. 
######

#1.3############################################################################

pred_highest <- c()
for(i in 1:366){
  pred_highest[i] <- max(preds[,i])
}
data13 <- cbind(data,t(pred12),pred_interval, pred_highest)
ggplot(data13)+
  geom_point(aes(x=time, y=temp))+
  geom_line(aes(x=time, y=t(pred12)),color="Blue",size=1)+
  geom_line(aes(x=time, y=lower),color="Red",linetype="dashed",size=1)+
  geom_line(aes(x=time, y=upper),color="Red",linetype="dashed",size=1)+
  geom_line(aes(x=time,y=pred_highest),color="Cyan",size=1)

  #1.4############################################################################




