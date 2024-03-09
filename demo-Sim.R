###   Simulation under unobserved confounders     ###
# devtools::install_github("socket778/XBCF")
# devtools::install_github("xnie/rlearner")
## required package: XBCF, grf, rlearner, glmnet, sparseGAM
rm(list=ls())
set.seed(123)

library(rlearner)
library(gam)
library(MCMCpack)
source("BCS-default.R")


# settings 
n <- 300  # the number of samples
p <- 20    #  number of covariates 

# covariates 
Z <- cbind(rnorm(n), rnorm(n))
X <- cbind(rnorm(n), rnorm(n), rnorm(n), rbinom(n, 1, 0.5), sample(1:3, n, replace=T))
add.p <- p - 5
if(add.p>0){
  for(j in 1:add.p){
    X <- cbind(X, rnorm(n))
  }
}

# treatment effect
tau <- 1+2*X[,2]*X[,5]+0.5*Z[,1]

# prognostic function
mu <- (-6) + (5-3*X[,5]) + 6*abs(X[,3]-1) + Z[,2]^2

# treatment assignment 
Tr <- rbinom(n, 1, 1/(1+exp(-0.5*Z[,1])))

# observed values
sig <- 1
Y <- mu + Tr*tau + sig*rnorm(n)


## Estimation of HTE
# Causal forest
folds <- sort(seq(n)%%10) + 1    # 10-fold cross validation
fit_CF <- grf::causal_forest(X=X, Y=Y, W=Tr, num.trees=2000, tune.parameters="all", clusters=folds)
Pred_CF <- predict(fit_CF, X, estimate.variance=T)
Est_CF <- Pred_CF[,1]
CI_CF <- rbind(Pred_CF[,1]-1.96*sqrt(Pred_CF[,2]), Pred_CF[,1]+1.96*sqrt(Pred_CF[,2]))

# X-learner
XL_fit <- xboost(as.matrix(X), Tr, Y)
Est_XL <- predict(XL_fit, as.matrix(X))

# R-learner
RL_fit <- rboost(as.matrix(X), Tr, Y)
Est_RL <- predict(RL_fit, as.matrix(X))

# Bayesian causal forest
bcf_fit <- XBCF::XBCF(y=as.matrix(Y), z=as.matrix(Tr), x_con=as.matrix(X), x_mod=as.matrix(X), 
                      pcat_con=3, pcat_mod=3, pihat=as.matrix(rep(0.5,n)), num_sweeps=1500, 
                      burnin=500, n_trees_con=100, n_trees_mod=100, parallel=T)
tau_pos <- XBCF::predictTauDraws(bcf_fit, X)     # (n, mc)-matrix of posterior samples of tau (HTE) 
Est_BCF <- apply(tau_pos, 1, mean)
CI_BCF <- apply(tau_pos, 1, quantile, prob=c(0.025, 0.975)) 

# Bayesian causal synthesis 
BCS_fit <- BCS_default(Y=Y, X=X, Tr=Tr, m=10, mc=1500, burn=500)
Est_BCS <- apply(BCS_fit$HTE, 2, mean)
CI_BCS <- apply(BCS_fit$HTE, 2, quantile, prob=c(0.025, 0.975))


## RMSE
Est <- cbind(Est_CF, Est_XL, Est_RL, Est_BCF, Est_BCS)
sqrt( apply((Est-tau)^2, 2, mean) ) 


## Coverage probability
100*mean(CI_CF[1,]<tau & CI_CF[2,]>tau)
100*mean(CI_BCF[1,]<tau & CI_BCF[2,]>tau)
100*mean(CI_BCS[1,]<tau & CI_BCS[2,]>tau)

# Average length
mean(CI_CF[2,]-CI_CF[1,])
mean(CI_BCF[2,]-CI_BCF[1,])
mean(CI_BCS[2,]-CI_BCS[1,])


# plot 
plot(tau, Est_BCS, ylim=range(Est_BCS, Est_BCF, Est_CF, CI_BCS), pch=20, cex=1.5,
     xlab="HTE (true)", ylab="HTE (estimated)")
points(tau, Est_BCF, col=2, pch=2)
points(tau, Est_CF, col=4, pch=17)
abline(0, 1)
for(i in 1:n){
  lines(x=rep(tau[i], 2), y=CI_BCS[,i], lty=1)
}
legend("topleft", legend=c("BCS", "BCF", "CF"), pch=c(20,2,17), col=c(1,2,4), lty=c(1,NA,NA))

