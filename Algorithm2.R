
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

bayesian.SSD2 <- function(true.hyp = 1, n.steps = 20, m = 1000, eff.size = .8, t.points = c(0,1,2,3,4)) {
  
  start <- Sys.time()

  true.hyp <- 1   # which hypothesis is true? H0: b=0; H1: b>0
  n.steps <- 100  # number of different sample sizes to evaluate
  m <- 1000       # number of datasets
  eff.size <- .2  # effect size (Cohen's d): .2 (small), .5 (medium), .8 (large)
  
  t.points <- c(0,1,2,3,4) # time of measurements
  n <- length(t.points)    # number of measurements per subject
  
  mbf10 <- numeric(n.steps) # storage for mean BFs 
  mbf01 <- numeric(n.steps) # storage for mean BFs 
  mbf.u <- numeric(n.steps) # storage for mean BFs 
  mbf.c <- numeric(n.steps) # storage for mean BFs
  mpmp.a <- numeric(n.steps) # storage for mean PMPs
  mpmp.b <- numeric(n.steps) # storage for mean PMPs
  mpmp.c <- numeric(n.steps) # storage for mean PMPs
  
  quantbf10 <- numeric(n.steps) # storage for quantiles in vectors containing BFs
  prop.bf10 <- numeric(n.steps) # storage for proportions of BFs exceeding a certain threshold
  prop.bf01 <- numeric(n.steps) # storage for proportions of BFs exceeding a certain threshold
  medbf10 <- numeric(n.steps) # storage for median BFs
  medbf01 <- numeric(n.steps) # storage for median BFs
  medbf.u <- numeric(n.steps) # storage for median BFs
  medbf.c <- numeric(n.steps) # storage for median BFs
  medpmp.a <- numeric(n.steps) # storage for median PMPs
  medpmp.b <- numeric(n.steps) # storage for median PMPs
  medpmp.c <- numeric(n.steps) # storage for median PMPs
  
  ss.seq <- rep(NA, n.steps)
  
  models <- rep(list(rep(vector("list", m))),n.steps)
  
  set.seed(1234)      # for reproducibility
  #seeds <- vector("list", n.steps) 
  
  sigmasq.u0 <- 1   # variance of individual deviation from treatment intercept 
  sigmasq.u1 <- 1   # variance of individual deviation from treatment slope
  sigma.u0.u1 <- 0  # covariance between sigmasq.u0 and sigmasq.u1. If positive, then individuals with higher(lower) initial values tend to have a higher (lower) rate of change over time.
  sigmasq.e <- 1    # error variance
  
  bf10 <- rep(NA,m) # storage for BFs
  bf01 <- rep(NA,m) # storage for BFs
  bf.u <- rep(NA,m) # storage for BFs
  bf.c <- rep(NA,m) # storage for BFs
  pmp.a <- rep(NA,m) # storage for PMPs
  pmp.b <- rep(NA,m) # storage for PMPs
  pmp.c <- rep(NA,m) # storage for PMPs
  
  Nmin <- 10
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
    
    beta0 <- rep(0, N[j]*n) # average y at t0 
    beta1 <- rep(0, N[j]*n) # average increase for x=0
    beta2 <- eff.size * sqrt(sigmasq.u1) * rep(true.hyp==1, N[j]*n) # average difference in slopes between conditions
    
    for (i in 1:m) {
      #seeds[[i]] <- .Random.seed
      multinorm <- mvrnorm(n=N[j], mu=c(0,0), matrix(c(sigmasq.u0, sigma.u0.u1, sigma.u0.u1, sigmasq.u1), nrow=2, ncol=2)) # draw individual deviation from treatment intercept and slope from a multivariate normal distribution with mean 0.
      u0 <- rep(multinorm[,1], each=n)
      u1 <- rep(multinorm[,2], each=n)
      e <- rnorm(N[j]*n, 0, sqrt(sigmasq.e))
      
      y <- beta0 + u0 + beta1*t + beta2*treat*t + u1*t + e
      dat <- data.frame(dat0, y)
      
      #inter <- lme(y ~ t + t:treat, random =~ t | id, data = dat, control = lmeControl(opt="optim"))
      models[[j]][[i]] <- lmer(y ~ t + t:treat + (t | id), data = dat, control = lmerControl(calc.derivs = F))
      #if (isSingular(models[[j]][[i]])) {next} # uncomment if singular models are not to be included
      
      est <- models[[j]][[i]]@beta[3]
      names(est) <- c("t:treat")
      sig <- list(as.matrix(vcov(models[[j]][[i]])[3,3]))
      
      result <- bain(est, hypotheses <- "t:treat>0;t:treat=0", n = N[j], Sigma = sig, group_parameters = 1, joint_parameters = 0) # double check whether group and joint pars are correct
      
      bf10[i] <- result$BFmatrix[1,2]
      bf01[i] <- result$BFmatrix[2,1]
      bf.u[i] <- result$fit$BF.u[1]
      bf.c[i] <- result$fit$BF.c[1]
      pmp.a[i] <- result$fit$PMPa[1]
      pmp.b[i] <- result$fit$PMPb[1]
      pmp.c[i] <- result$fit$PMPc[1]
    }
    
    mbf10[j] <- mean(bf10, na.rm=T)
    mbf01[j] <- mean(bf01, na.rm=T)
    mbf.u[j] <- mean(bf.u, na.rm=T)
    mbf.c[j] <- mean(bf.c, na.rm=T)
    mpmp.a[j] <- mean(pmp.a, na.rm=T)
    mpmp.b[j] <- mean(pmp.b, na.rm=T)
    mpmp.c[j] <- mean(pmp.c, na.rm=T)
    
    medbf10[j] <- median(bf10, na.rm=T)
    quantbf10[j] <- quantile(bf10, probs=0.2, na.rm=T)
    prop.bf10[j] <- length(bf10[bf10>3])/m
    prop.bf01[j] <- length(bf01[bf01>3])/m
    medbf01[j] <- median(bf01, na.rm=T)
    medbf.u[j] <- median(bf.u, na.rm=T)
    medbf.c[j] <- median(bf.c, na.rm=T)
    medpmp.a[j] <- median(pmp.a, na.rm=T)
    medpmp.b[j] <- median(pmp.b, na.rm=T)
    medpmp.c[j] <- median(pmp.c, na.rm=T)
    
    print(N[j])
    if(N[j]==Nmin+1) {break}
    ifelse (prop.bf10[j]>.8, Nmax <- N[j], Nmin <- N[j])
      
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



