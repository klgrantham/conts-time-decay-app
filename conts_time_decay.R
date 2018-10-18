# Compilation of functions for continuous-time decaying correlation Shiny app
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(ltsa)
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)

# Equations for variance of treatment effect
# Mean and individual level

vartheta_mean <- function(Vi, Xmat){
  # Calculates the variance of the treatment effect, theta, for a model at the
  # cluster mean level with a particular treatment schedule
  #
  # Inputs:
  # Xmat - a K x T matrix of the treatment schedule (note: all elements either 0 or 1)
  # Vi - a T x T variance matrix for one cluster
  
  K <- nrow(Xmat)
  Tp <- ncol(Xmat)
  Xvec <- as.vector(t(Xmat))
  Vi_inv <- solve(Vi)
  var <- 1/(t(Xvec) %*% (diag(1,K) %x% Vi_inv) %*% Xvec - 
              colSums(Xmat) %*% Vi_inv %*% (matrix(colSums(Xmat),nrow=Tp, ncol=1))/K)
  return(var)
}

vartheta_ind_vec <- function(r, Tp, m, rho0, Xmat){
  # Calculates the variance of the treatment effect, theta, for a model at the
  # individual level with a particular treatment schedule
  #
  # Inputs:
  # Xmat - a vector of K x T matrices of the treatment schedule (note: all elements either 0 or 1)
  # Virow1 - first row of length Tm of the covariance matrix for one cluster
  
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  Virow1 <- as.vector(cbind((sig2E + sig2CP), sig2CP*(r^(matrix(1:(Tp*m-1), nrow=1, ncol=Tp*m-1)/m))))
  Vi <- toeplitz(Virow1)
  Vi_inv <- TrenchInverse(Vi)
  vars <- laply(Xmat, vartheta, Vi_inv)
  return(vars)
}

vartheta <- function(Xmat, Vi_inv) {
  # Returns variance of treatment effect estimator for an inverse covariance
  # matrix and a design matrix
  
  K <- nrow(Xmat)
  Tp <- ncol(Xmat)
  m <- nrow(Vi_inv)/Tp
  
  Q <- Xmat %x% t(rep(1,m))
  B <- colSums(Xmat) %x% rep(1,m)
  C <- diag(Tp) %x% rep(1,m)
  term1 <- sum(diag(Q %*% Vi_inv %*% t(Q))) # Previously: t(D) %*% (diag(1,K) %x% Vi_inv) %*% D, where D <- Xvec %x% rep(1,m)
  term2 <- t(B) %*% Vi_inv %*% C
  term3 <- solve(t(C) %*% Vi_inv %*% C)
  term4 <- t(C) %*% Vi_inv %*% B
  var <- 1/(term1 - (1/K)*term2 %*% term3 %*% term4)
  return(var)
}

# Variance matrices under different models
# Hussey & Hughes, discrete time decay, continuous time decay

HHVi <- function(rho0, Tp, m){
  # Constructs the variance matrix for a single cluster, Vi, under the
  # Hussey & Hughes model (2007), at either the cluster mean level or
  # at the individual level
  #
  # Inputs:
  # Tp - number of time periods
  # m - number of individuals per cluster
  # rho0 - proportion of total variation attributed to cluster random effects

  totalvar <- 1
  tau2 <- rho0*totalvar
  sig2E <- totalvar - tau2
  sig2 <- sig2E/m # subject-specific variance at mean level
  Vi <- diag(sig2,Tp) + tau2*matrix(1, nrow=Tp, ncol=Tp)
  return(Vi)
}

expdecayVi <- function(r, Tp, m, rho0){
  # Constructs the variance matrix for a single cluster, Vi, under the
  # exponential decay model (Kasza et al 2017), at either the cluster-
  # period mean level or at the individual level
  #
  # Inputs:
  # Tp - number of time periods
  # m - number of individuals per cluster
  # rho0 - proportion of total variation attributed to cluster-period random effects

  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/m
  Vi <- diag(sig2,Tp) +
        sig2CP*(r^abs(matrix(1:Tp, nrow=Tp, ncol=Tp, byrow=FALSE) -
                      matrix(1:Tp, nrow=Tp, ncol=Tp, byrow=TRUE)))    
  return(Vi)
}


# Design matrices (or treatment schedules) for different trial designs
# Stepped wedge (SW), cluster randomised crossover (CRXO), parallel, parallel w/baseline

SWdesmat <- function(Tp, N=Tp-1){
  Xsw <- matrix(data=0, ncol=Tp, nrow=(Tp-1))
  for(i in 1:(Tp-1)){
    Xsw[i,(i+1):Tp] <- 1
  }
  if(N%%(Tp-1) != 0) stop('N must be a multiple of Tp-1')
  Xsw %x% rep(1, N/(Tp-1))
}

crxodesmat <- function(Tp, N=Tp){
  if(N%%2 != 0) stop('N must be even')
  Xcrxo <- matrix(data=0, ncol=Tp, nrow=2)
  Xcrxo[1, seq(1,Tp,2)] <- 1
  Xcrxo[2, seq(2,Tp,2)] <- 1
  Xcrxo %x% rep(1, N/2)
}

plleldesmat <- function(Tp, N=Tp){
  if(N%%2 != 0) stop('N must be even')
  Xpllel <- matrix(data=0, ncol=Tp, nrow=2)
  Xpllel[1,] <- 1
  Xpllel %x% rep(1, N/2)
}

pow <- function(vars, effsize, siglevel=0.05){
  z <- qnorm(siglevel/2)
  pow <- pnorm(z + sqrt(1/vars)*effsize)
  return(pow)
}

ErhoUC <- function(r, Tp, m, rho0_CD, N){
  if(r==1) ErhoUC <- rho0_CD
  else{
    sig2E <- 1 - rho0_CD
    quadsum <- rho0_CD*(Tp*m - (2*r^(1/m)*(1 - r^Tp - m*Tp*(1 - r^(1/m))))/(r^(1/m) - 1)^2)
    dblsum <- rho0_CD*(m - (2*r^(1/m)*(1 - r - m*(1 - r^(1/m))))/(r^(1/m) - 1)^2)
    Esig2CP <- (1/(N*Tp*m-N-Tp+1))*(((N*m-1)/(Tp*m^2))*quadsum + (1/m^2)*dblsum - N*rho0_CD)
    Esig2E <- sig2E + (1/(N*Tp*m-N-T+1))*(N*Tp*m*rho0_CD - ((N-1)/(Tp*m))*quadsum - (Tp/m)*dblsum)
    ErhoUC <- Esig2CP/(Esig2CP + Esig2E)
  }
  return(ErhoUC)
}

# Generate variance results for user-defined trial configuration, with progress bar
generate_var_results_prog <- function(Tp, N, m, rho0_CD, updateProgress = NULL) {
  # Calculates the variance of the treatment effect estimator under the models:
  #    continuous time (ct), discrete time (dt), Hussey & Hughes (HH)
  # with trial designs:
  #    stepped wedge (SW),
  #    cluster randomised crossover (CRXO),
  #    parallel (pllel)
  #
  # Inputs:
  # Tp - number of time periods in the trial
  # m - number of subjects measured in each time period
  # rho0 - base correlation between a pair of subjects
  #
  # Example usage: vals <- generate_var_results(Tp=4, m=50, rho0=0.035)
  
  # Set vector of r values (Decay = 1-r)
  rs <- seq(0.5, 1, 0.01)
  # Specify the covariance matrices under the different models
  if (is.function(updateProgress)) {
    updateProgress()
  }
  dtvarmat <- llply(rs, expdecayVi, Tp, m, rho0_CD)
  
  # Calculate expected value of rho0_UC
  rho0_UC <- laply(rs, ErhoUC, Tp, m, rho0_CD, N)
  HHvarmat <- llply(rho0_UC, HHVi, Tp, m)
  
  # Get the variances of the treatment effect estimator under the
  # different models and designs
  SWXmat <- SWdesmat(Tp, N)
  crxoXmat <- crxodesmat(Tp, N)
  pllelXmat <- plleldesmat(Tp, N)
  Xmats <- list(SWXmat, crxoXmat, pllelXmat)
  
  if (is.function(updateProgress)) {
    updateProgress()
  }
  ctres <- laply(rs, vartheta_ind_vec, Tp=Tp, m=m, rho0=rho0_CD, Xmat=Xmats)

  if (is.function(updateProgress)) {
    updateProgress()
  }
  varvals <- data.frame(decay=1-rs,
                        rhoUC = rho0_UC,
                        ctSW = ctres[,1],
                        ctcrxo = ctres[,2],
                        ctpllel = ctres[,3],
                        dtSW = laply(dtvarmat, vartheta_mean, Xmat=SWXmat),
                        dtcrxo = laply(dtvarmat, vartheta_mean, Xmat=crxoXmat),
                        dtpllel = laply(dtvarmat, vartheta_mean, Xmat=pllelXmat),
                        HHSW = laply(HHvarmat, vartheta_mean, Xmat=SWXmat),
                        HHcrxo = laply(HHvarmat, vartheta_mean, Xmat=crxoXmat),
                        HHpllel = laply(HHvarmat, vartheta_mean, Xmat=pllelXmat))
  return(varvals)
}