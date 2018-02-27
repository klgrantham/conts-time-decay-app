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
  var <- 1/(t(Xvec) %*% (diag(1,K) %x% solve(Vi)) %*% Xvec - 
              colSums(Xmat) %*% solve(Vi) %*% (matrix(colSums(Xmat),nrow=Tp, ncol=1))/K)
  return(var)
}

vartheta_ind <- function(Vi, Xmat, Toeplitz=TRUE){
  # Calculates the variance of the treatment effect, theta, for a model at the
  # individual level with a particular treatment schedule
  #
  # Inputs:
  # Xmat - a K x T matrix of the treatment schedule (note: all elements either 0 or 1)
  # Vi - a Tm x Tm variance matrix for one cluster
  
  K <- nrow(Xmat)
  Tp <- ncol(Xmat)
  m <- nrow(Vi)/Tp
  # If continuous time matrix, use Toeplitz inversion algorithm
  if(Toeplitz){ # Could check for Vi[1,2]!=Vi[1,3] but w/ simulated times won't be Toeplitz
    Vi_inv <- TrenchInverse(Vi)
  } else{
    Vi_inv <- solve(Vi)
  }
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

vartheta_ind_inv <- function(Vi_inv, Xmat, Toeplitz=TRUE){
  # Calculates the variance of the treatment effect, theta, for a model at the
  # individual level with a particular treatment schedule
  #
  # Inputs:
  # Vi_inv - a Tm x Tm inverse covariance matrix for one cluster
  # Xmat - a K x T matrix of the treatment schedule (note: all elements either 0 or 1)
  
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

HHVi <- function(Tp, m, rho0, meanlvl=TRUE){
  # Constructs the variance matrix for a single cluster, Vi, under the
  # Hussey & Hughes model (2007), at either the cluster mean level or
  # at the individual level
  #
  # Inputs:
  # Tp - number of time periods
  # m - number of individuals per cluster
  # rho0 - proportion of total variation attributed to cluster random effects
  # meanlvl - boolean for whether to construct the variance matrix for the
  #           cluster mean level (default) or the individual level
  
  totalvar <- 1
  tau2 <- rho0*totalvar
  sig2E <- totalvar - tau2
  sig2 <- sig2E/m # subject-specific variance at mean level
  if(meanlvl==TRUE){
    Vi <- diag(sig2,Tp) + tau2*matrix(1, nrow=Tp, ncol=Tp)
  }
  else{
    Vi <- diag(sig2E,Tp*m) + tau2*matrix(1, nrow=Tp*m, ncol=Tp*m)
  }
  return(Vi)
}

expdecayVi <- function(r, Tp, m, rho0, meanlvl=TRUE){
  # Constructs the variance matrix for a single cluster, Vi, under the
  # exponential decay model (Kasza et al 2017), at either the cluster-
  # period mean level or at the individual level
  #
  # Inputs:
  # Tp - number of time periods
  # m - number of individuals per cluster
  # rho0 - proportion of total variation attributed to cluster-period random effects
  # meanlvl - boolean for whether to construct the variance matrix for the
  #           cluster-period mean level (default) or the individual level
  
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/m
  if(meanlvl==TRUE){
    Vi <- diag(sig2,Tp) +
      sig2CP*(r^abs(matrix(1:Tp, nrow=Tp, ncol=Tp, byrow=FALSE) -
                      matrix(1:Tp, nrow=Tp, ncol=Tp, byrow=TRUE)))    
  }
  else{
    Vi <- diag(sig2E,Tp*m) +
      sig2CP*(r^abs(matrix(rep(1:Tp, each=m), nrow=Tp*m, ncol=Tp*m, byrow=FALSE) -
                      matrix(rep(1:Tp, each=m), nrow=Tp*m, ncol=Tp*m, byrow=TRUE)))
  }
  return(Vi)
}

expdecayVicont <- function(r, Tp, m, rho0, meanlvl=FALSE){
  # Constructs the variance matrix for a single cluster, Vi, under the
  # exponential decay model in continuous itme, at either the cluster-
  # period mean level or at the individual level
  #
  # Inputs:
  # Tp - number of time periods
  # m - number of individuals per cluster
  # rho0 - proportion of total variation attributed to cluster-period random effects
  # meanlvl - boolean for whether to construct the variance matrix for the
  #           cluster-period mean level or the individual level (default)
  
  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP
  sig2 <- sig2E/m
  if(meanlvl==TRUE){
    if(r==1){
      Vi <- expdecayVi(r=r, Tp=Tp, m=m, rho0=rho0, meanlvl=meanlvl) # Reverts to Hussey and Hughes
    }
    else{
      vars <- (1/m^2)*(m + 2*((r^(1/m))*(-m*(r^(1/m)) + m + r - 1))/((r^(1/m)) - 1)^2) # to be multiplied by sig2CP and added to sig2
      covars <- (1/m^2)*(((r^((1/m)-1))*(r-1)^2)/((r^(1/m)) - 1)^2) # to be multiplied by sig2CP and r^{|j-j'|}
      mat <- matrix(covars, nrow=Tp, ncol=Tp) - diag(covars,Tp) + diag(vars,Tp)
      Vi <- diag(sig2,Tp) +
        (sig2CP*(r^abs(matrix(1:Tp, nrow=Tp, ncol=Tp, byrow=FALSE) -
                         matrix(1:Tp, nrow=Tp, ncol=Tp, byrow=TRUE)))*mat)  
    }
  }
  else{
    Vi <- diag(sig2E,Tp*m) +
      sig2CP*(r^(abs(matrix(rep(1:(Tp*m)), nrow=Tp*m, ncol=Tp*m, byrow=FALSE) -
                       matrix(rep(1:(Tp*m)), nrow=Tp*m, ncol=Tp*m, byrow=TRUE))/m))
  }
  return(Vi)
}

# Design matrices (or treatment schedules) for different trial designs
# Stepped wedge (SW), cluster randomised crossover (CRXO), parallel, parallel w/baseline

SWdesmat <- function(Tp){
  Xsw <- matrix(data=0, ncol=Tp, nrow=(Tp-1))
  for(i in 1:(Tp-1)){
    Xsw[i,(i+1):Tp] <- 1
  }
  return(Xsw)
}

crxodesmat <- function(Tp){
  if(Tp%%2==0) {
    nclust <- Tp
  } else {
    nclust <- Tp-1
  }
  Xcrxo <- matrix(data=0, ncol=Tp, nrow=nclust)
  Xcrxo[1:nclust/2, seq(1,Tp,2)] <- 1
  Xcrxo[(nclust/2 + 1):nclust, seq(2,Tp,2)] <- 1
  return(Xcrxo)
}

plleldesmat <- function(Tp){
  if(Tp%%2==0) {
    nclust <- Tp
  } else {
    nclust <- Tp-1
  }
  Xpllel <- matrix(data=0, ncol=Tp, nrow=nclust)
  Xpllel[1:nclust/2,] <- 1
  return(Xpllel)
}

pllelbasedesmat <- function(Tp){
  if(Tp%%2==0) {
    nclust <- Tp
  } else {
    nclust <- Tp-1
  }
  Xpllelbase <- matrix(data=0, ncol=Tp, nrow=nclust)
  Xpllelbase[1:nclust/2, 2:Tp] <- 1
  return(Xpllelbase)
}

# Generate variance results for user-defined trial configuration

generate_var_results <- function(Tp, m, rho0){
  # Calculates the variance of the treatment effect under the models:
  #    continuous time (ct), discrete time (dt), Hussey & Hughes (HH)
  # with trial designs:
  #    stepped wedge (SW), cluster randomised crossover (CRXO),
  #    parallel (pllel), parallel with baseline (pllelbase)
  #
  # Inputs:
  # Tp - number of time periods in the trial
  # m - number of subjects measured in each time period
  # rho0 - base correlation between a pair of subjects
  #
  # Example usage: vals <- comparison_plots(Tp=4, m=50, rho0=0.035)
  
  rs <- seq(0.5, 1, 0.01) # Note: May need to expand range from 0 to 1?
  # Specify the variance matrices under the different models
  ctvarmat <- llply(rs, expdecayVicont, Tp, m, rho0, meanlvl=FALSE)
  dtvarmat <- llply(rs, expdecayVi, Tp, m, rho0, meanlvl=TRUE)
  HHvarmat <- HHVi(Tp, m, rho0, meanlvl=TRUE)
  
  # Get the variances of the treatment effect under the
  # different models and designs
  # Note: Still need to expand the Xmats according to nclust?
  # Scale the non-SW variances by (Tp/(Tp-1)) to account for uneven clusters across designs
  scalefactor <- Tp/(Tp-1)
  varvals <- data.frame(decay = 1-rs,
                        ctSW = laply(ctvarmat, vartheta_ind, Xmat=SWdesmat(Tp)),
                        ctcrxo = scalefactor*laply(ctvarmat, vartheta_ind, Xmat=crxodesmat(Tp)),
                        ctpllel = scalefactor*laply(ctvarmat, vartheta_ind, Xmat=plleldesmat(Tp)),
                        ctpllelbase = scalefactor*laply(ctvarmat, vartheta_ind, Xmat=pllelbasedesmat(Tp)),
                        dtSW = laply(dtvarmat, vartheta_mean, Xmat=SWdesmat(Tp)),
                        dtcrxo = scalefactor*laply(dtvarmat, vartheta_mean, Xmat=crxodesmat(Tp)),
                        dtpllel = scalefactor*laply(dtvarmat, vartheta_mean, Xmat=plleldesmat(Tp)),
                        dtpllelbase = scalefactor*laply(dtvarmat, vartheta_mean, Xmat=pllelbasedesmat(Tp)),
                        HHSW = vartheta_mean(Vi=HHvarmat, Xmat=SWdesmat(Tp)),
                        HHcrxo = scalefactor*vartheta_mean(Vi=HHvarmat, Xmat=crxodesmat(Tp)),
                        HHpllel = scalefactor*vartheta_mean(Vi=HHvarmat, Xmat=plleldesmat(Tp)),
                        HHpllelbase = scalefactor*vartheta_mean(Vi=HHvarmat, Xmat=pllelbasedesmat(Tp)))
  return(varvals)
}

generate_var_results_prog <- function(Tp, m, rho0, updateProgress = NULL){
  # Calculates the variance of the treatment effect under the models:
  #    continuous time (ct), discrete time (dt), Hussey & Hughes (HH)
  # with trial designs:
  #    stepped wedge (SW), cluster randomised crossover (CRXO),
  #    parallel (pllel), parallel with baseline (pllelbase)
  #
  # Inputs:
  # Tp - number of time periods in the trial
  # m - number of subjects measured in each time period
  # rho0 - base correlation between a pair of subjects
  #
  # Example usage: vals <- comparison_plots(Tp=4, m=50, rho0=0.035)
  
  rs <- seq(0.5, 1, 0.01) # Note: May need to expand range from 0 to 1?
  # Specify the variance matrices under the different models
  if (is.function(updateProgress)) {
    updateProgress()
  }
  ctvarmat <- llply(rs, expdecayVicont, Tp, m, rho0, meanlvl=FALSE)
  
  if (is.function(updateProgress)) {
    updateProgress()
  }
  dtvarmat <- llply(rs, expdecayVi, Tp, m, rho0, meanlvl=TRUE)
  
  if (is.function(updateProgress)) {
    updateProgress()
  }
  HHvarmat <- HHVi(Tp, m, rho0, meanlvl=TRUE)
  
  # Get the variances of the treatment effect under the
  # different models and designs
  # Note: Still need to expand the Xmats according to nclust?
  # Scale the non-SW variances by (Tp/(Tp-1)) to account for uneven clusters across designs
  scalefactor <- Tp/(Tp-1)
  
  # Split up computations and include checks for updateProgress()
  if (is.function(updateProgress)) {
    updateProgress()
  }
  ctSW <- laply(ctvarmat, vartheta_ind, Xmat=SWdesmat(Tp))
  
  if (is.function(updateProgress)) {
    updateProgress()
  }
  ctcrxo <- scalefactor*laply(ctvarmat, vartheta_ind, Xmat=crxodesmat(Tp))

  if (is.function(updateProgress)) {
    updateProgress()
  }
  ctpllel <- scalefactor*laply(ctvarmat, vartheta_ind, Xmat=plleldesmat(Tp))
  
  if (is.function(updateProgress)) {
    updateProgress()
  }
  dtSW <- laply(dtvarmat, vartheta_mean, Xmat=SWdesmat(Tp))
  dtcrxo <- scalefactor*laply(dtvarmat, vartheta_mean, Xmat=crxodesmat(Tp))
  dtpllel <- scalefactor*laply(dtvarmat, vartheta_mean, Xmat=plleldesmat(Tp))
  
  if (is.function(updateProgress)) {
    updateProgress()
  }
  HHSW <- vartheta_mean(Vi=HHvarmat, Xmat=SWdesmat(Tp))
  HHcrxo <- scalefactor*vartheta_mean(Vi=HHvarmat, Xmat=crxodesmat(Tp))
  HHpllel <- scalefactor*vartheta_mean(Vi=HHvarmat, Xmat=plleldesmat(Tp))
  
  varvals <- data.frame(decay = 1-rs, ctSW = ctSW, ctcrxo = ctcrxo, ctpllel = ctpllel,
                        dtSW = dtSW, dtcrxo = dtcrxo, dtpllel = dtpllel,
                        HHSW = HHSW, HHcrxo = HHcrxo, HHpllel = HHpllel)
  return(varvals)
}


compare_designs <- function(df.long, ylabel){
  names(df.long)[dim(df.long)[2]] <- "value" # Assumes last column to be plotted
  p <- ggplot(data=df.long, aes(x=decay, y=value, colour=Design, linetype=Design)) +
    geom_line(size=1.5) +
    scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF"),
                       labels = c("CRXO", "Parallel", "SW")) +
    scale_linetype_manual(values = c("solid", "dotdash", "dashed"),
                          labels = c("CRXO", "Parallel", "SW")) +
    xlab("Decay (1 - r)") +
    ylab(ylabel) +
#    labs(title=bquote(paste(.(Tp), " periods, ", .(m), " subjects/cluster-period"))) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, size=14),
          axis.title=element_text(size=16), axis.text=element_text(size=16),
          legend.key.width = unit(1.5, "cm"),
          legend.title=element_text(size=16), legend.text=element_text(size=16),
          legend.position="bottom")
  return(p)
}

plotvar_ct <- function(df.results){
  # Variance, continuous time, all designs
  ctvarvals <- df.results %>%
    select(decay, starts_with('ct')) %>%
    filter(decay<=0.5)
  ctvarvals_long <- gather(data=ctvarvals, key=Design, value=Variance,
                           ctSW:ctpllel, convert=TRUE)

  compare_designs(ctvarvals_long, "Variance")
}

plotrelvar_ctvdt <- function(df.results){
  # Relative variance, continuous vs discrete time, all designs
  ctvdtvarvals <- df.results %>%
    mutate(ratioSW=ctSW/dtSW, ratiocrxo=ctcrxo/dtcrxo, ratiopllel=ctpllel/dtpllel) %>%
    select(decay, starts_with('ratio')) %>%
    filter(decay<=0.5)
  ctvdtvarvals_long <- gather(data=ctvdtvarvals, key=Design, value=Relative_variance,
                              ratioSW:ratiopllel, convert=TRUE)

  compare_designs(ctvdtvarvals_long, "Relative variance") + geom_hline(aes(yintercept=1))
}

plotrelvar_ctvHH <- function(df.results){
  # Relative variance, continuous time vs HH, all designs
  ctvHHvarvals <- df.results %>%
    mutate(ratioSW=ctSW/HHSW, ratiocrxo=ctcrxo/HHcrxo, ratiopllel=ctpllel/HHpllel) %>%
    select(decay, starts_with('ratio')) %>%
    filter(decay<=0.5)
  ctvHHvarvals_long <- gather(data=ctvHHvarvals, key=Design, value=Relative_variance,
                              ratioSW:ratiopllel, convert=TRUE)
  
  compare_designs(ctvHHvarvals_long, "Relative variance") + geom_hline(aes(yintercept=1))
}  
