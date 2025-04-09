#===============================================================================
# Function for data generation and Bayes Factor calculation for a single data set
#===============================================================================

# Hypotheses to be tested:
# H0: beta2 = 0
# H1: beta2 > 0
# where beta2 is the coefficient of interaction of time and treatment condition (treatment effect)

# The function "getbf" uses the following arguments:
# N = the total sample size (number of subjects)
# t.points = position of the measurement occasions in time 
# var.u0 = intercept variance
# var.u1 = slope variance
# cov = covariance between intercept and slope variance
# var.e = error variance
# eff.size = effect size defined as beta/sqrt(var.u1), where beta is the coefficient of interaction
# fraction = fraction of information used to specify prior, b = fraction/N
# Neff = if "worst": effective sample size = N, if "best": effective sample size = N*n, 
# where n = number of measurement occasions
# log.grow = indicates whether to use logarithmic (TRUE) or linear growth (FALSE)
# hyp = for which hypotheses should the BF be calculated ("H0", "H1", "both" or "h0", "h1", "b")

# Note 1: this function requires the function "get_neff" to be loaded into the global environment
# Note 2: this function requires the packages "lme4" and "MASS" to be installed

#-------------------------------------------------------------------------------

getbf <- function(N, t.points, var.u0, var.u1, cov, var.e, eff.size, fraction, log.grow, hyp, beta1){
  
  n <- length(t.points)  # number of measurement occasions
  ifelse(log.grow==F,  # if logarithmic growth is used, take log of t.points
         t <- rep(t.points, N),  # create time variable t
         ifelse(min(t.points)==0,  # if the first timepoint is zero, we add 1 to all timepoints because log(0) is undefined
                t <- rep(log(t.points+1), N), 
                t <- rep(log(t.points), N)  # otherwise, just use log(t.points)
         )
  )
  id <- rep(seq_len(N), each=n)  # create ID variable
  treat <- as.numeric(as.character(gl(n=2, k=n, length=N*n, labels=c(0,1))))  # create treatment variable
  dat0 <- data.frame(id, treat, t)  # combine into data frame
  beta2.H1 <- eff.size * sqrt(var.u1)  # create coefficient of interaction under H1 from effect size; beta2=0|H0
  
  multinorm <- MASS::mvrnorm(n=2*N, mu=c(0,0), matrix(c(var.u0, cov, cov, var.u1), nrow=2, ncol=2))  # draw random effects
  
  # generate data under H0 -----------------------------------------------------
  u0.H0 <- rep(multinorm[1:(nrow(multinorm)/2),1], each=n)  # random intercepts for H0
  u1.H0 <- rep(multinorm[1:(nrow(multinorm)/2),2], each=n)  # random slopes for H0
  e.H0 <- rnorm(N*n, 0, sqrt(var.e))  # error variance for H0
  y.H0 <- beta1*t + u0.H0 + 0*treat*t + u1.H0*t + e.H0  # create outcome variable y under H0
  dat.H0 <- data.frame(dat0, y.H0)  # add y under H0 to data frame
  
  # generate data under H1 -----------------------------------------------------
  u0.H1 <- rep(multinorm[(nrow(multinorm)/2+1):(nrow(multinorm)),1], each=n)  # random intercepts for H1
  u1.H1 <- rep(multinorm[(nrow(multinorm)/2+1):(nrow(multinorm)),2], each=n)  # random slopes for H1
  e.H1 <- rnorm(N*n, 0, sqrt(var.e))  # error variance for H1
  y.H1 <- beta1*t + u0.H1 + beta2.H1*treat*t + u1.H1*t + e.H1  # create outcome variable y under H1
  dat.H1 <- data.frame(dat0, y.H1)  # add y under H1 to data frame
  
  if(hyp == "both" | hyp == "b"){  # calculate BF both both H0 and H1
    
    # fit MLM to dataset under H0
    models.H0 <- lme4::lmer(y.H0 ~ t + treat + t:treat + (1 + t | id), data = dat.H0, REML=F, control = lme4::lmerControl(calc.derivs = F))  # fit MLM model under H0
    est.H0 <- models.H0@beta[4]  # extract estimate for coefficient of interaction under H0
    sig.H0 <- vcov(models.H0)[4,4]  # extract residual variance under H0  
    
    # calculate Neff and b
    n_eff.H0 <- get_neff(model=models.H0)
    b.H0 <- fraction/n_eff.H0
    
    # calculate fits and complexities under H0
    comp0.H0 <- dnorm(0, mean=0, sd=sqrt(sig.H0/b.H0))
    fit0.H0 <- dnorm(0, mean=est.H0, sd=sqrt(sig.H0))
    comp1.H0 <- 1-pnorm(0, mean=0, sd=sqrt(sig.H0/b.H0))
    fit1.H0 <- 1-pnorm(0, mean=est.H0, sd=sqrt(sig.H0))
    
    # calculate BFs under H0
    BFu.H0 <- fit0.H0/comp0.H0
    BFc.H0 <- BFu.H0
    BFu1.H0 <- fit1.H0/comp1.H0
    BFs.H0 <- BFu.H0/BFu1.H0
    pmp.a.H0 <- BFu.H0/(BFu.H0 + BFu1.H0)
    
    # do the same for H1 -------------------------------------------------------
    
    # fit MLM to dataset under H1
    models.H1 <- lme4::lmer(y.H1 ~ t + treat + t:treat + (1 + t | id), data = dat.H1, REML=F, control = lme4::lmerControl(calc.derivs = F))  # fit MLM model under H1
    est.H1 <- models.H1@beta[4]  # extract estimate for coefficient of interaction under H1
    sig.H1 <- vcov(models.H1)[4,4]  # extract residual variance under H0
    
    # calculate Neff and b
    n_eff.H1 <- get_neff(model=models.H1)
    b.H1 <- fraction/n_eff.H1
    
    # calculate fits and complexities under H1
    comp1.H1 <- 1-pnorm(0, mean=0, sd=sqrt(sig.H1/b.H1))
    fit1.H1 <- 1-pnorm(0, mean=est.H1, sd=sqrt(sig.H1))
    comp0.H1 <- dnorm(0, mean=0, sd=sqrt(sig.H1/b.H1))
    fit0.H1 <- dnorm(0, mean=est.H1, sd=sqrt(sig.H1))
    
    # calculate BFs under H1
    BFu.H1 <- fit1.H1/comp1.H1
    BFc.H1 <- (fit1.H1/comp1.H1) / ((1-fit1.H1)/(1-comp1.H1))
    BFu0.H1 <- fit0.H1/comp0.H1
    BFs.H1 <- BFu.H1/BFu0.H1
    pmp.a.H1 <- BFu.H1/(BFu.H1 + BFu0.H1)
    
    # return BF01 and BF10
    return(output = list(BF01 = BFs.H0,
                         BF10 = BFs.H1,
                         BF0u = BFu.H0,
                         BF1u = BFu.H1,
                         BF0c = BFc.H0,
                         BF1c = BFc.H1))
    
  } else if (hyp == "h0" | hyp == "H0"){  # calculate BF for H0 only -----------
    
    # fit MLM to dataset under H0
    models.H0 <- lme4::lmer(y.H0 ~ t + treat + t:treat + (1 + t | id), data = dat.H0, REML=F, control = lme4::lmerControl(calc.derivs = F))
    est.H0 <- models.H0@beta[4]  # extract estimate for coefficient of interaction under H0
    sig.H0 <- vcov(models.H0)[4,4]  # extract residual variance under H0  
    
    # calculate Neff and b
    n_eff <- get_neff(model=models.H0)
    b <- fraction/n_eff
    
    # calculate fits and complexities under H0
    comp0.H0 <- dnorm(0, mean=0, sd=sqrt(sig.H0/b))
    fit0.H0 <- dnorm(0, mean=est.H0, sd=sqrt(sig.H0))
    comp1.H0 <- 1-pnorm(0, mean=0, sd=sqrt(sig.H0/b))
    fit1.H0 <- 1-pnorm(0, mean=est.H0, sd=sqrt(sig.H0))
    
    # calculate BFs under H0
    BFu.H0 <- fit0.H0/comp0.H0
    BFc.H0 <- BFu.H0
    BFu1.H0 <- fit1.H0/comp1.H0
    BFs.H0 <- BFu.H0/BFu1.H0
    pmp.a.H0 <- BFu.H0/(BFu.H0 + BFu1.H0)    
    
    # return BF01
    return(list(BF01 = BFs.H0,
                BF0u = BFu.H0,
                BF0c = BFc.H0))
    
  } else if (hyp == "h1" | hyp == "H1"){  # calculate BF for H1 only -----------
    
    # fit MLM to dataset under H1
    models.H1 <- lme4::lmer(y.H1 ~ t + treat + t:treat + (1 + t | id), data = dat.H1, REML=F, control = lme4::lmerControl(calc.derivs = F))
    est.H1 <- models.H1@beta[4]  # extract estimate for coefficient of interaction under H1
    sig.H1 <- vcov(models.H1)[4,4]  # extract residual variance under H0
    
    # calculate Neff and b
    n_eff <- get_neff(model=models.H1)
    b <- fraction/n_eff
    
    # calculate fits and complexities under H1
    comp1.H1 <- 1-pnorm(0, mean=0, sd=sqrt(sig.H1/b))
    fit1.H1 <- 1-pnorm(0, mean=est.H1, sd=sqrt(sig.H1))
    comp0.H1 <- dnorm(0, mean=0, sd=sqrt(sig.H1/b))
    fit0.H1 <- dnorm(0, mean=est.H1, sd=sqrt(sig.H1))
    
    # calculate BFs under H1
    BFu.H1 <- fit1.H1/comp1.H1
    BFc.H1 <- (fit1.H1/comp1.H1) / ((1-fit1.H1)/(1-comp1.H1))
    BFu0.H1 <- fit0.H1/comp0.H1
    BFs.H1 <- BFu.H1/BFu0.H1
    pmp.a.H1 <- BFu.H1/(BFu.H1 + BFu0.H1)
    
    # return BF01 and BF10
    return(list(BF10 = BFs.H1,
                BF1u = BFu.H1,
                BF1c = BFc.H1))
  }
  
}

# END OF FUNCTION --------------------------------------------------------------
