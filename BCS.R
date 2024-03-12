###------------------------------------------------------------------### 
###          Code for Bayesian Causal Synthesis                      ###
###   (with user specified Gaussian distribution for HTE)            ###
###------------------------------------------------------------------###

## INPUT
# Y: n-dimensional response vector 
# X: (n,p)-covariate matrix 
# Tr: n-dimensional vector of treatment indicators (1:treatment, 0:alternative)
# PS: n-dimensional vector of propensity scores (optinal) 
# mF: (n,J)-matrix of point estimation of heterogeneous treatment effects (HTE)
# vF: matrix of (asymptotic) variance associated with HTE estimates
# m: number of nearest neighbors in NNGP (default is 10)
# mc: length of Markov Chain Monte Carlo
# burn: length of burn-in period
# X_test: covariate matrix in test sample (optional)
# mF_test: HTE estimates in test samples (optional)
# mF_test: associated variance of HTE estimates in test samples (optional)
# Remark 1: exponential covariance function is used in NNGP.  
# Remark 2: prior means of beta (model weight) is 1/J.

## OUTPUT
# posterior samples of beta (model weights), sigma (error variance), 
# tau and phi (scale and range parameters in NNGP)

BCS <- function(Y, X, Tr, PS=NULL, mF, vF, m=10, mc=2000, burn=500, Phi.set=NULL, 
                        print=T, prediction=F, X_test=NULL, mF_test=NULL, vF_test=NULL){
  ## preparation
  n <- length(Y)     # sample size
  if(dim(X)[2]==1){ X <- cbind(0, X) }
  p <- dim(X)[2]
  MF <- cbind(1, mF)
  VF <- cbind(NA, vF)
  J <- dim(MF)[2] 
  if(is.null(PS)){ PS <- rep(0.5, n) }
  prior_var <- 1     # hyperparameter in inverse gamma prior for variance parameters 
  prior_ratio <- 1    
  tau_prior <- cbind(c(prior_var, rep(prior_ratio*prior_var, J-1)), rep(prior_var, J))
  Beta_center <- 1/(J-1)
  
  ## correlation function 
  CF <- function(dist.mat, phi){ exp(-dist.mat/phi) }   
  
  ## neighbor location for prognostic term 
  Z <- cbind(X, PS)
  Z <- Z + 0.01*rnorm(prod(dim(Z)))   # add small noise to prevent the same covariates
  Z_mm <- apply(Z, 2, mean)
  Z_ss <- apply(Z, 2, sd)
  Z <- t( (t(Z) - Z_mm)/Z_ss )
  dd2 <- as.matrix( dist(Z) )
  NN2 <- matrix(NA, n, m)
  for(i in 2:n){
    if(i<=m){  NN2[i, 1:(i-1)] <- sort(order(dd2[i,1:(i-1)]))  }
    if(i>m){   NN2[i,] <- sort(order(dd2[i,1:(i-1)])[1:m])  }
  }
  
  ## neighbor location for treatment effect term 
  X <- Z[,1:p]    # scaled covariate matrix
  dd1 <- as.matrix( dist(X) )
  NN1 <- matrix(NA, n, m)
  for(i in 2:n){
    if(i<=m){  NN1[i, 1:(i-1)] <- sort(order(dd1[i,1:(i-1)]))  }
    if(i>m){   NN1[i,] <- sort(order(dd1[i,1:(i-1)])[1:m])  }
  }
  
  ## set of children indices 
  In <- function(x, i){ i %in% x }
  UU1 <- list()
  UU2 <- list()
  for(i in 1:n){
    UU1[[i]] <- (1:n)[apply(NN1, 1, In, i=i)]
    UU2[[i]] <- (1:n)[apply(NN2, 1, In, i=i)]
  }
  
  ## set of spatial range parameters
  L <- length(Phi.set)
  if(is.null(Phi.set)){
    L <- 10
    Phi.range <- c(0, median(dist(Z)))
    Phi.set <- seq(Phi.range[1], Phi.range[2], length=L+1)[-1]
  }
  
  ## set of prior parameters in NNGP for each candidate value of Phi (treatment effect term)
  BB.set1 <- list()    # coefficients for neighbor locations
  FF.set1 <- list()    # conditional variances 
  for(l in 1:L){
    BB1 <- matrix(NA, n, m)
    FF1 <- c()
    FF1[1] <- 1
    for(i in 2:n){
      if(i<=m){
        mat <- CF(as.matrix( dist(as.matrix(X[c(NN1[i,1:(i-1)], i),]))), Phi.set[l])
        C1 <- solve(mat[1:(i-1), 1:(i-1)])
        C2 <- mat[i, 1:(i-1)]
        FF1[i] <- mat[i, i] - t(C2)%*%C1%*%C2   # conditional variance  
        BB1[i,1:(i-1)] <- as.vector(C2%*%C1)   # coefficient for neighbour variables 
      }
      if(i>m){
        mat <- CF(as.matrix(dist(X[c(NN1[i,], i),])), Phi.set[l])
        C1 <- solve(mat[1:m, 1:m])
        C2 <- mat[m+1, 1:m]
        FF1[i] <- mat[m+1, m+1] - t(C2)%*%C1%*%C2   # conditional variance  
        BB1[i,] <- as.vector(C2%*%C1)   # coefficient for neighbour variables 
      }
    }
    BB.set1[[l]] <- BB1
    FF.set1[[l]] <- FF1
  }
  
  ## set of prior parameters in NNGP for each candidate value of Phi (prognostic term)
  BB.set2 <- list()    # coefficients for neighbor locations
  FF.set2 <- list()    # conditional variances 
  for(l in 1:L){
    BB2 <- matrix(NA, n, m)
    FF2 <- c()
    FF2[1] <- 1
    for(i in 2:n){
      if(i<=m){
        mat <- CF(as.matrix( dist(Z[c(NN2[i,1:(i-1)], i),])), Phi.set[l])
        C1 <- solve(mat[1:(i-1), 1:(i-1)])
        C2 <- mat[i, 1:(i-1)]
        FF2[i] <- mat[i, i] - t(C2)%*%C1%*%C2   # conditional variance  
        BB2[i,1:(i-1)] <- as.vector(C2%*%C1)    # coefficients for neighbour variables 
      }
      if(i>m){
        mat <- CF(as.matrix(dist(Z[c(NN2[i,], i),])), Phi.set[l])
        C1 <- solve(mat[1:m, 1:m])
        C2 <- mat[m+1, 1:m]
        FF2[i] <- mat[m+1, m+1] - t(C2)%*%C1%*%C2   # conditional variance  
        BB2[i,] <- as.vector(C2%*%C1)               # coefficients for neighbour variables 
      }
    }
    BB.set2[[l]] <- BB2
    FF.set2[[l]] <- FF2
  }
  
  ## matrices/arrays to store posterior samples 
  Beta.pos <- array(NA, c(mc, n, J))    # model weight
  Mu.pos <- matrix(NA, mc, n)           # baseline term
  Phi.pos <- matrix(NA, mc, J+1)        # range parameters in NNGP
  Tau.pos <- matrix(NA, mc, J+1)        # scale parameters in NNGP
  Sig.pos <- c()                        # error variance  
  TE.pos <- matrix(NA, mc, n)           # HTE for observed sample
  
  ## initial values
  Beta <- cbind(rep(0, n), matrix(Beta_center, n, J-1))   # model weight
  Mu <- rep(0, n)                         # baseline
  Phi.index <- rep(round(L/2), J+1)       # range parameters 
  Tau <- rep(1, J+1)                      # scale parameters 
  Sig <- 0.5                              # error variance 
  TP <- MF                                # latent values of HTE
  
  ## Markov Chain Monte Carlo  
  BB1 <- array(NA, c(n, m, J))
  FF1 <- matrix(NA, n, J)
  BB2 <- matrix(NA, n, m)
  FF2 <- rep(NA, n)
  for(r in 1:mc){
    TP_Tr <- TP*Tr
    # update Beta (update centered model weight)
    sY <- Y - Mu - apply(TP_Tr[,-1], 1, sum)*Beta_center
    sBeta <- t( t(Beta) - c(0, rep(Beta_center, J-1)) )   # centered model weight 
    for(k in 1:J){
      BB1[,,k] <- BB.set1[[Phi.index[k]]]
      FF1[,k] <- FF.set1[[Phi.index[k]]]
    }
    tFF <- t( Tau[1:J]^2*t(FF1) )
    for(i in 1:n){
      Us <- UU1[[i]]
      if(length(Us)==0){ 
        V <- TP_Tr[i,]%*%t(TP_Tr[i,])/Sig^2 + diag(1/tFF[i,])
        prior.mean <- apply(na.omit(BB1[i,,]*sBeta[NN1[i,],]), 2, sum) / tFF[i,] 
        pos.mean <- TP_Tr[i,]*sY[i]/Sig^2 + prior.mean
      }
      if(length(Us)>0){
        sNN <- matrix(NN1[Us,], length(Us), m)     # set of parents ID including "i" (matrix object)
        prior.prec <- matrix(0, J, J)
        for(k in 1:J){
          prior.prec[k, k] <- 1/tFF[i, k] + sum( na.omit((BB1[Us,,k]^2/tFF[Us, k])[sNN==i]) )
        }
        V <- TP_Tr[i,]%*%t(TP_Tr[i,])/Sig^2 + prior.prec
        length.Us <- length(Us)
        prior.mean <- c()
        for(k in 1:J){
          mBeta <- matrix(sBeta[as.vector(sNN), k], length.Us, m)
          mBeta[sNN==i] <- 0
          mBeta[is.na(mBeta)] <- 0
          sBB <- BB1[Us,,k]
          sBB[is.na(sBB)] <- 0
          bb <- na.omit( (t(sBB))[t(sNN)==i] )
          sBB[sNN==i] <- 0
          a <- sBeta[Us, k] - apply(mBeta*sBB, 1, sum)
          prior.mean[k] <- sum(bb*a/tFF[Us,k])
        }
        prior.mean <- prior.mean + apply(na.omit(BB1[i,,]*sBeta[NN1[i,],]), 2, sum) / tFF[i,] 
        pos.mean <- TP_Tr[i,]*sY[i]/Sig^2 + prior.mean
      }
      IV <- solve(V)
      sBeta[i,] <- mvrnorm(1, IV%*%pos.mean, IV) 
      Beta[i,] <- sBeta[i,] + c(0, rep(Beta_center, J-1)) 
    }
    Beta.pos[r,,] <- Beta
    if(mean(Beta)>1000){ break }
    
    # update Mu (baseline term)
    sY <- Y - apply(TP_Tr*Beta, 1, sum)
    BB2 <- BB.set2[[Phi.index[J+1]]]
    FF2 <- FF.set2[[Phi.index[J+1]]]
    tFF <- Tau[J+1]^2*FF2
    for(i in 1:n){
      Us <- UU2[[i]]
      if(length(Us)==0){ 
        V <- 1/Sig^2 + 1/tFF[i]
        prior.mean <- sum( na.omit(BB2[i,]*Mu[NN2[i,]]) ) / tFF[i] 
        pos.mean <- sY[i]/Sig^2 + prior.mean
      }
      if(length(Us)>0){
        sNN <- matrix(NN2[Us,], length(Us), m)     # set of parents ID including "i" (matrix object)
        prior.prec <- 1/tFF[i] + sum( na.omit((BB2[Us,]^2/tFF[Us])[sNN==i]) )
        V <- 1/Sig^2 + prior.prec
        length.Us <- length(Us)
        mMu <- matrix(Mu[as.vector(sNN)], length.Us, m)
        mMu[sNN==i] <- 0
        mMu[is.na(mMu)] <- 0
        sBB <- BB2[Us,]
        sBB[is.na(sBB)] <- 0
        bb <- na.omit( (t(sBB))[t(sNN)==i] )
        sBB[sNN==i] <- 0
        a <- Mu[Us] - apply(mMu*sBB, 1, sum)
        prior.mean <- sum(bb*a/tFF[Us]) + sum(na.omit(BB2[i,]*Mu[NN2[i,]])) / tFF[i] 
        pos.mean <- sY[i]/Sig^2 + prior.mean
      }
      IV <- 1/V
      Mu[i] <- rnorm(1, pos.mean*IV, sqrt(IV))
    }
    Mu.pos[r,] <- Mu
    
    # update TP (latent value of HTE)
    for(j in 2:J){
      A <- 1/(Tr^2*Beta[,j]^2/Sig^2 + 1/VF[,j])
      resid <- Y - Mu - apply(as.matrix(Tr*TP[,-j]*Beta[,-j]), 1, sum)
      B <- Tr*Beta[,j]*resid/Sig^2 + MF[,j]/VF[,j]
      TP[,j] <- rnorm(n, A*B, sqrt(A))
    }
    
    # update Sig
    resid <- Y - Mu - Tr*apply(TP*Beta, 1, sum)
    Sig <- sqrt( rinvgamma(1, prior_var+n/2, prior_var+sum(resid^2)/2) )
    Sig.pos[r] <- Sig
    
    # update Tau (for model weights)
    sBeta <- t( t(Beta) - c(0, rep(Beta_center, J-1)) )
    prior.sBeta <- matrix(NA, n, J)
    for(k in 1:J){
      prior <- BB1[,,k]*matrix(sBeta[NN1,k], n, m)
      prior[is.na(prior)] <- 0
      prior.sBeta[,k] <- apply(prior, 1, sum)
    }
    ss <- apply((sBeta - prior.sBeta)^2/FF1, 2, sum)
    Tau[1:J] <- sqrt( rinvgamma(J, tau_prior[,1]+n/2, tau_prior[,2]+ss/2) )
    
    # update Tau (for baseline term)
    prior <- BB2*matrix(Mu[NN2], n, m)
    prior[is.na(prior)] <- 0
    prior.Mu <- apply(prior, 1, sum)
    ss <- sum( (Mu - prior.Mu)^2/FF2 )
    Tau[J+1] <- sqrt( rinvgamma(1, prior_var+n/2, prior_var+ss/2) )
    Tau.pos[r, ] <- Tau
    
    # update Phi (for model weights)
    for(k in 1:J){
      new.index <- Phi.index[k] + sample(c(1, -1), 1) 
      if(new.index<1){ new.index <- 1 }
      if(new.index>L){ new.index <- L }
      new.BB <- BB.set1[[new.index]]
      new.FF <- FF.set1[[new.index]]
      new.prior.sBeta <- new.BB*matrix(sBeta[NN1,k], n, m)
      new.prior.sBeta[is.na(new.prior.sBeta)] <- 0
      new.prior.sBeta <- apply(new.prior.sBeta, 1, sum)
      L1 <- sum( dnorm(x=sBeta[,k], mean=prior.sBeta[,k], sd=Tau[k]*sqrt(FF1[,k]), log=T) )
      L2 <- sum( dnorm(x=sBeta[,k], mean=new.prior.sBeta, sd=Tau[k]*sqrt(new.FF), log=T) )
      pp <- min(1, exp(L2-L1))
      if(runif(1)<pp){ Phi.index[k] <- new.index }
    }
    
    # update Phi (for baseline term)
    new.index <- Phi.index[J+1] + sample(c(1, -1), 1) 
    if(new.index<1){ new.index <- 1 }
    if(new.index>L){ new.index <- L }
    new.BB <- BB.set2[[new.index]]
    new.FF <- FF.set2[[new.index]]
    new.prior.Mu <- new.BB*matrix(Mu[NN2], n, m)
    new.prior.Mu[is.na(new.prior.Mu)] <- 0
    new.prior.Mu <- apply(new.prior.Mu, 1, sum)
    L1 <- sum( dnorm(x=Mu, mean=prior.Mu, sd=Tau[J+1]*sqrt(FF2), log=T) )
    L2 <- sum( dnorm(x=Mu, mean=new.prior.Mu, sd=Tau[J+1]*sqrt(new.FF), log=T) )
    pp <- min(1, exp(L2-L1))
    if(runif(1)<pp){ Phi.index[J+1] <- new.index }
    Phi.pos[r,] <- Phi.set[Phi.index]
    
    # update TE (treatment effect for observed sample)
    TE.pos[r,] <- apply(TP*Beta, 1, sum)
    
    # print
    if(round(10*r/mc)==(10*r/mc) & print ){ 
      print( paste0(round(100*r/mc), "% of MCMC has been done!") ) 
    }
  }
  
  ## omit burn-in samples
  om <- 1:burn
  Beta.pos <- Beta.pos[-om,,]
  Mu.pos <- Mu.pos[-om,]
  Sig.pos <- Sig.pos[-om]
  Phi.pos <- Phi.pos[-om,]
  Tau.pos <- Tau.pos[-om,]
  TE.pos <- TE.pos[-om,]
  
  ## Posterior samples of HTE for test sample (optional)
  if(prediction){
    MC <- mc - burn
    MF_test <- cbind(1, mF_test)
    VF_test <- cbind(NA, vF_test)
    tau2.pos <- Tau.pos^2
    X_test <- t((t(X_test)-Z_mm[1:p])/Z_ss[1:p])
    n_test <- dim(X_test)[1]
    Pred.pos <- matrix(NA, MC, n_test)
    for(i in 1:n_test){
      dd <- sqrt( apply((X_test[i,] - t(X))^2, 2, sum) )
      NN <- sort(order(dd)[1:m]) 
      dd2 <- as.matrix(dist(X[NN,]))
      dd3 <- sqrt( apply((X_test[i,] - t(X[NN,]))^2, 2, sum) )
      fBeta <- matrix(NA, MC, J)
      for(k in 1:J){
        for(r in 1:MC){
          C1 <- solve( tau2.pos[r,k]*exp(-dd2/Phi.pos[r,k]) )
          C2 <- tau2.pos[r,k]*exp(-dd3/Phi.pos[r,k])
          fFF <- tau2.pos[r,k] - t(C2)%*%C1%*%C2   
          fBB <- as.vector(C2%*%C1) 
          if(k==1){ 
            mm <- sum( fBB*Beta.pos[r, NN, k] )
            fBeta[r,k] <- rnorm(1, mm, sqrt(fFF)) 
          }
          if(k>1){ 
            mm <- sum( fBB*(Beta.pos[r, NN, k]-Beta_center) )
            fBeta[r,k] <- Beta_center + rnorm(1, mm, sqrt(fFF)) 
          }
        }
      }
      FP_test <- rep(1, MC)
      for(j in 2:J){
        FP_test <- cbind(FP_test, rnorm(MC, MF_test[i,j], sqrt(VF_test[i,j])))
      }
      Pred.pos[,i] <- apply(fBeta*FP_test, 1, sum)
    }
  }
  
  ## Summary
  if(prediction==F){ Pred.pos <- NULL }
  result <- list(beta=Beta.pos, mu=Mu.pos, sig=Sig.pos, phi=Phi.pos, tau=Tau.pos, 
                 TE=TE.pos, pred=Pred.pos)
  return(result)
}


