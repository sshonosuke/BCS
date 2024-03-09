###--------------------------------------------------------------------------### 
###           Demonstration of Bayesian Causal Synthesis (BCS)               ### 
###                         using ACIC data                                  ###                
###--------------------------------------------------------------------------###
rm(list=ls())
set.seed(1)

library(MASS)
library(grf)
library(rlearner)
source("BCS-default.R")


## Data generation 
# Covariate 
load("Covariate.RData")
names(X) <- paste0("X", c(1, 3, 10, 14, 15, 21, 24, 43))

# Settings
j <- 6    # scenario (parameter)     # 1-8
n <- 800    # number of sub-samples (train + test) 
N <- dim(X)[1]
X <- X[sort(sample(1:N, n)),]

# Parameter scenarios
xi_set <- c(1/3, 2)
eta_set <- c(0.25, 1.25)
kap_set <- cbind(c(0.5, 0), c(3, -1))
sel <- expand.grid(c(1,2), c(1,2), c(1,2))[,3:1]

xi <- xi_set[sel[j,1]]
eta <- eta_set[sel[j,2]]
kap1 <- kap_set[1, sel[j,3]]
kap2 <- kap_set[2, sel[j,3]]

# True functions 
f <- X$X1 + X$X43 + 0.3*(X$X10 - 1)
Pi <- ( 1+exp(kap1*f+kap2) )^(-1)
mu <- -sin( pnorm(Pi) ) + X$X43 
tau <- xi*( X$X3*X$X24 + (X$X14-1) + (X$X15-1) )
sigma <- 0.4 + (X$X21-1)/15
sigy <- eta*sqrt( var(mu+Pi*tau) )
Tr <- rbinom(n, 1, Pi)

# non-additive scenario
Yt <- mu + tau*Tr + sigy*rnorm(n)
a <- mean(Yt)
b <- 1.25*sd(Yt)
m1 <- (mu + tau - a) / b
m0 <- (mu - a) / b
Tau_true <- 13*pnorm(m1/sqrt(sigy^2/b^2+1)) - 13*pnorm(m0/sqrt(sigy^2/b^2+1))
Y <- 13*pnorm(Yt, a, sqrt(b))

# split train and test data
n_train <- 500
n_test <- n - n_train
Tau_true_train <- Tau_true[1:n_train]
Tau_true_test <- Tau_true[-(1:n_train)]

X_test <- X[-(1:n_train),]

Y <- Y[1:n_train]
X <- X[1:n_train,]
Tr <- Tr[1:n_train]
Pi <- Pi[1:n_train]


## Estimation of HTE
# Causal forest
folds <- sort(seq(n_train)%%10) + 1    # 10-fold cross validation
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
                      pcat_con=3, pcat_mod=3, pihat=as.matrix(Pi), num_sweeps=1500, 
                      burnin=500, n_trees_con=100, n_trees_mod=100, parallel=T)
tau_pos <- XBCF::predictTauDraws(bcf_fit, X)     # (n, mc)-matrix of posterior samples of tau (HTE) 
Est_BCF <- apply(tau_pos, 1, mean)
CI_BCF <- apply(tau_pos, 1, quantile, prob=c(0.025, 0.975)) 


# BCS (default version)
BCS_fit <- BCS_default(Y=Y, X=X, Tr=Tr, PS=Pi, m=10, mc=1500, burn=500, prediction=T, X_test=X_test)
Est_BCS <- apply(BCS_fit$HTE, 2, mean)
CI_BCS <- apply(BCS_fit$HTE, 2, quantile, prob=c(0.025, 0.975))
Est_BCS_test <- apply(BCS_fit$HTE_test, 2, mean)
CI_BCS_test <- apply(BCS_fit$HTE_test, 2, quantile, prob=c(0.025, 0.975))



## RMSE
Est <- cbind(Est_CF, Est_XL, Est_RL, Est_BCF, Est_BCS)
sqrt( apply((Est-Tau_true_train)^2, 2, mean) ) 
sqrt( mean((Est_BCS_test-Tau_true_test)^2) )   # for test sample 


## Coverage probability
100*mean(CI_CF[1,]<Tau_true_train & CI_CF[2,]>Tau_true_train)
100*mean(CI_BCF[1,]<Tau_true_train & CI_BCF[2,]>Tau_true_train)
100*mean(CI_BCS[1,]<Tau_true_train & CI_BCS[2,]>Tau_true_train)
100*mean(CI_BCS_test[1,]<Tau_true_test & CI_BCS_test[2,]>Tau_true_test)   # for test sample

# Average length
mean(CI_CF[2,]-CI_CF[1,])
mean(CI_BCF[2,]-CI_BCF[1,])
mean(CI_BCS[2,]-CI_BCS[1,])
mean(CI_BCS_test[2,]-CI_BCS_test[1,])   # for test samples


