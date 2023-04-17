################################################################################
###################### ALGORITHM 2 -- Binary Search ############################
################################################################################

library(tidyr)       # pipes
library(dplyr)       # pipes
library(ggplot2)     # plots
library(lme4)        # fit multilevel model
library(mgcv)        # extracting vcov matrices
library(bain)        # Bayesian estimation
library(MASS)        # multinorm - already included in lme4?

bayesian.SSD2 <- function(n.steps = 20, m = 1000, eff.size = .8, t.points = c(0,1,2,3,4), BF.thresh=.8) {
  
  start <- Sys.time()

  n.steps <- 20  # number of different sample sizes to evaluate
  m <- 100       # number of datasets
  eff.size <- .8  # effect size (Cohen's d): .2 (small), .5 (medium), .8 (large)
  BF.thresh=.8
  
  t.points <- c(0,1,2,3,4) # time of measurements
  n <- length(t.points)    # number of measurements per subject
  
  prop.bf10.H0 <- numeric(n.steps) # storage for proportions of BFs exceeding a certain threshold
  prop.bf01.H0 <- numeric(n.steps) # storage for proportions of BFs exceeding a certain threshold
  medbf10.H0 <- numeric(n.steps) # storage for median BFs
  medbf01.H0 <- numeric(n.steps) # storage for median BFs
  medbf.1u.H0 <- numeric(n.steps) # storage for median BFs
  medbf.1c.H0 <- numeric(n.steps) # storage for median BFs
  medpmp.a1.H0 <- numeric(n.steps) # storage for median PMPs
  medpmp.b1.H0 <- numeric(n.steps) # storage for median PMPs
  medpmp.c1.H0 <- numeric(n.steps) # storage for median PMPs
  
  prop.bf10.H1 <- numeric(n.steps) # storage for proportions of BFs exceeding a certain threshold
  prop.bf01.H1 <- numeric(n.steps) # storage for proportions of BFs exceeding a certain threshold
  medbf10.H1 <- numeric(n.steps) # storage for median BFs
  medbf01.H1 <- numeric(n.steps) # storage for median BFs
  medbf.1u.H1 <- numeric(n.steps) # storage for median BFs
  medbf.1c.H1 <- numeric(n.steps) # storage for median BFs
  medpmp.a1.H1 <- numeric(n.steps) # storage for median PMPs
  medpmp.b1.H1 <- numeric(n.steps) # storage for median PMPs
  medpmp.c1.H1 <- numeric(n.steps) # storage for median PMPs
  
  ss.seq <- rep(NA, n.steps)
  
  modelsH0 <- rep(list(rep(vector("list", m))),n.steps)
  modelsH1 <- rep(list(rep(vector("list", m))),n.steps)
  
  set.seed(1234)      # for reproducibility
  #seeds <- vector("list", n.steps) 
  
  sigmasq.u0 <- 1   # variance of individual deviation from treatment intercept 
  sigmasq.u1 <- 1   # variance of individual deviation from treatment slope
  sigma.u0.u1 <- 0  # covariance between sigmasq.u0 and sigmasq.u1. If positive, then individuals with higher(lower) initial values tend to have a higher (lower) rate of change over time.
  sigmasq.e <- 1    # error variance
  
  bf10.H0 <- rep(NA,m) # storage for BFs
  bf01.H0 <- rep(NA,m) # storage for BFs
  bf.1u.H0 <- rep(NA,m) # storage for BFs
  bf.1c.H0 <- rep(NA,m) # storage for BFs
  pmp.a1.H0 <- rep(NA,m) # storage for PMPs
  pmp.b1.H0 <- rep(NA,m) # storage for PMPs
  pmp.c1.H0 <- rep(NA,m) # storage for PMPs
  
  bf10.H1 <- rep(NA,m) # storage for BFs
  bf01.H1 <- rep(NA,m) # storage for BFs
  bf.1u.H1 <- rep(NA,m) # storage for BFs
  bf.1c.H1 <- rep(NA,m) # storage for BFs
  pmp.a1.H1 <- rep(NA,m) # storage for PMPs
  pmp.b1.H1 <- rep(NA,m) # storage for PMPs
  pmp.c1.H1 <- rep(NA,m) # storage for PMPs
  
  Nmin <- 20
  Nmax <- 1000
  
  N <- numeric(n.steps)
  
  for(j in 1:n.steps){
    
    N[j] <- round((Nmin + Nmax)/2, digits = 0)
    
    # create data vectors
    y <- rep(NA, N[j])    # data storage
    t <- rep(t.points, N[j])
    id <- rep(seq_len(N[j]), each=n)
    treat <- as.numeric(as.character(gl(n=2, k=n, length=N[j]*n, labels=c(0,1))))
    dat0 <- data.frame(id, treat, t)
    
    beta0 <- 0 # average y at t0 
    beta1 <- 0 # average increase for x=0
    beta2H0 <- 0                           # average difference in slopes between conditions under H0
    beta2H1 <- eff.size * sqrt(sigmasq.u1) # average difference in slopes between conditions under H1
    
    for (i in 1:m) {
      #seeds[[i]] <- .Random.seed
      multinorm <- mvrnorm(n=N[j], mu=c(0,0), matrix(c(sigmasq.u0, sigma.u0.u1, sigma.u0.u1, sigmasq.u1), nrow=2, ncol=2)) # draw individual deviation from treatment intercept and slope from a multivariate normal distribution with mean 0.
      u0 <- rep(multinorm[,1], each=n)
      u1 <- rep(multinorm[,2], each=n)
      e <- rnorm(N[j]*n, 0, sqrt(sigmasq.e))
      
      yH0 <- beta0 + u0 + beta1*t + beta2H0*treat*t + u1*t + e # data-generating mechanism under H0
      yH1 <- beta0 + u0 + beta1*t + beta2H1*treat*t + u1*t + e # data-generating mechanism under H1
      datH0 <- data.frame(dat0, yH0)
      datH1 <- data.frame(dat0, yH1)
      
      #inter <- lme(y ~ t + t:treat, random =~ t | id, data = dat, control = lmeControl(opt="optim"))
      modelsH0[[j]][[i]] <- lmer(yH0 ~ t + t:treat + (t | id), data = datH0, control = lmerControl(calc.derivs = F))
      modelsH1[[j]][[i]] <- lmer(yH1 ~ t + t:treat + (t | id), data = datH1, control = lmerControl(calc.derivs = F))
      #if (isSingular(models[[j]][[i]])) {next} # uncomment if singular models are not to be included
      
      
      estH0 <- modelsH0[[j]][[i]]@beta[3]
      estH1 <- modelsH1[[j]][[i]]@beta[3]
      names(estH0) <- c("t:treat")
      names(estH1) <- c("t:treat")
      
      sigH0 <- list(as.matrix(vcov(modelsH0[[j]][[i]])[3,3]))
      sigH1 <- list(as.matrix(vcov(modelsH1[[j]][[i]])[3,3]))
      
      resultH0 <- bain(estH0, hypotheses <- "t:treat>0;t:treat=0", n = N[j], Sigma = sigH0, group_parameters = 1, joint_parameters = 0) 
      resultH1 <- bain(estH1, hypotheses <- "t:treat>0;t:treat=0", n = N[j], Sigma = sigH1, group_parameters = 1, joint_parameters = 0)
      
      bf10.H0[i] <- resultH0$BFmatrix[1,2]
      bf01.H0[i] <- resultH0$BFmatrix[2,1]
      bf.1u.H0[i] <- resultH0$fit$BF.u[1]
      bf.1c.H0[i] <- resultH0$fit$BF.c[1]
      pmp.a1.H0[i] <- resultH0$fit$PMPa[1]
      pmp.b1.H0[i] <- resultH0$fit$PMPb[1]
      pmp.c1.H0[i] <- resultH0$fit$PMPc[1]
      
      bf10.H1[i] <- resultH1$BFmatrix[1,2]
      bf01.H1[i] <- resultH1$BFmatrix[2,1]
      bf.1u.H1[i] <- resultH1$fit$BF.u[1]
      bf.1c.H1[i] <- resultH1$fit$BF.c[1]
      pmp.a1.H1[i] <- resultH1$fit$PMPa[1]
      pmp.b1.H1[i] <- resultH1$fit$PMPb[1]
      pmp.c1.H1[i] <- resultH1$fit$PMPc[1]
    }
    
    medbf10.H0[j] <- median(bf10.H0, na.rm=T)
    prop.bf10.H0[j] <- length(bf10.H0[bf10.H0>BFthresh])/m
    prop.bf01.H0[j] <- length(bf01.H0[bf01.H0>BFthresh])/m
    medbf01.H0[j] <- median(bf01.H0, na.rm=T)
    medbf.1u.H0[j] <- median(bf.1u.H0, na.rm=T)
    medbf.1c.H0[j] <- median(bf.1c.H0, na.rm=T)
    medpmp.a1.H0[j] <- median(pmp.a1.H0, na.rm=T)
    medpmp.b1.H0[j] <- median(pmp.b1.H0, na.rm=T)
    medpmp.c1.H0[j] <- median(pmp.c1.H0, na.rm=T)
    
    medbf10.H1[j] <- median(bf10.H1, na.rm=T)
    prop.bf10.H1[j] <- length(bf10.H1[bf10.H1>BFthresh])/m
    prop.bf01.H1[j] <- length(bf01.H1[bf01.H1>BFthresh])/m
    medbf01.H1[j] <- median(bf01.H1, na.rm=T)
    medbf.1u.H1[j] <- median(bf.1u.H1, na.rm=T)
    medbf.1c.H1[j] <- median(bf.1c.H1, na.rm=T)
    medpmp.a1.H1[j] <- median(pmp.a1.H1, na.rm=T)
    medpmp.b1.H1[j] <- median(pmp.b1.H1, na.rm=T)
    medpmp.c1.H1[j] <- median(pmp.c1.H1, na.rm=T)
    
    print(N[j])
    
    if(N[j]==Nmin+1) {break}
    
    ifelse(prop.bf01.H0[j]>.8 & prop.bf10.H1[j]>.8, Nmax <- N[j], Nmin <- N[j])
      
  }
  
  print(Sys.time() - start) 
  print(N[j])
  
  # par(mfrow=c(2,4))
  # plot(x=N, y=medbf10[1:length(N)], type="l", xlab="N", ylab="median BF10")
  # plot(x=N, y=medbf01[1:length(N)], type="l", xlab="N", ylab="median BF01")
  # plot(x=N, y=medbf.u[1:length(N)], type="l", xlab="N", ylab="median BF unconstr.")
  # plot(x=N, y=medbf.c[1:length(N)], type="l", xlab="N", ylab="median BF complement")
  # 
  # plot(x=N, y=mpmp.a[1:length(N)], type="l", xlab="N", ylab="mean PMPa")
  # plot(x=N, y=mpmp.b[1:length(N)], type="l", xlab="N", ylab="mean PMPb")
  # plot(x=N, y=mpmp.c[1:length(N)], type="l", xlab="N", ylab="mean PMPc")
  # plot(x=N, y=prop.bf10[1:length(N)], type="l", xlab="N", ylab="proportion of BFs larger 3")
  # abline(h=.8, col="red")
}



