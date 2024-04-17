#================================================================================
# Function for data generation and Bayes Factor calculation for a single data set
#================================================================================

# The function "getbf" uses the following arguments:
# N = the total sample size (number of subjects)
# t.points = position of the measurement occasions in time
# var.u0 = intercept variance
# var.u1 = slope variance
# cov = covariance between intercept and slope variance
# var.e = error variance
# eff.size = effect size defined as beta/sqrt(var.u1), where beta is the coefficient of interaction
# Bfthres = threshold a BF needs to exceed to be considered substantial
# fraction = fraction of information used to specify prior, b = fraction/N
# Neff = if "worst": effective sample size = N, if "best": effective sample size = N*n, 
# where n = number of measurement occasions
# log.grow = indicates whether to use logarithmic (TRUE) or linear growth (FALSE)

# Note that this function is required by the function "getpower".


getbf <- function(N, t.points, var.u0, var.u1, cov, var.e, eff.size, BFthres, fraction, Neff, log.grow){
  
  n <- length(t.points)  # number of measurement occasions
  ifelse(Neff=="worst",  # if Neff="worst", then Neff=N, otherwise Neff=N*n
         b <- fraction/N,  # b fraction to specify prior = fraction / Neff
         b <- fraction/N*n)
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
  
  multinorm <- MASS::mvrnorm(n=2*N, mu=c(0,0), matrix(c(var.u0, cov, cov, var.u1), nrow=2, ncol=2))  # draw random efects
  u0.H0 <- rep(multinorm[1:(nrow(multinorm)/2),1], each=n)  # random intercepts for H0
  u0.H1 <- rep(multinorm[(nrow(multinorm)/2+1):(nrow(multinorm)),1], each=n)  # random intercepts for H1
  u1.H0 <- rep(multinorm[1:(nrow(multinorm)/2),2], each=n)  # random slopes for H0
  u1.H1 <- rep(multinorm[(nrow(multinorm)/2+1):(nrow(multinorm)),2], each=n)  # random slopes for H1
  e.H0 <- rnorm(N*n, 0, sqrt(var.e))  # error variance for H0
  e.H1 <- rnorm(N*n, 0, sqrt(var.e))  # error variance for H1
  y.H0 <- u0.H0 + 0*treat*t + u1.H0*t + e.H0  # create outcome variable y under H0
  y.H1 <- u0.H1 + beta2.H1*treat*t + u1.H1*t + e.H1  # create outcome variable y under H1
  dat.H0 <- data.frame(dat0, y.H0)  # add y under H0 to data frame
  dat.H1 <- data.frame(dat0, y.H1)  # add y under H1 to data frame
  
  models.H0 <- lmer(y.H0 ~ t + t:treat + (t | id), data = dat.H0, control = lmerControl(calc.derivs = F))  # fit MLM model under H0
  est.H0 <- models.H0@beta[3]  # extract estimate for coefficient of interaction under H0
  sig.H0 <- vcov(models.H0)[3,3]  # extract residual variance under H0  
  
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
  
  models.H1 <- lmer(y.H1 ~ t + t:treat + (t | id), data = dat.H1, control = lmerControl(calc.derivs = F))  # fit MLM model under H1
  est.H1 <- models.H1@beta[3]  # extract estimate for coefficient of interaction under H1
  sig.H1 <- vcov(models.H1)[3,3]  # extract residual variance under H0
  
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
  
  # return power for H0 and H1
  return(output = list(BF01 = BFs.H0,
                       BF10 = BFs.H1))
}