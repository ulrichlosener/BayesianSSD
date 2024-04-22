
################################################################################
########## This function generates datasets under both hypotheses, #############
########## fits a multilevel model to each one and calculates BFs  #############
################################################################################

# This function works as follows: under each hypothesis, m datasets are generated.
# Then, a multilevel models is fitted to each dataset

# The function is a vectorized (i.e., faster) version

# The ingredients for this function are as follows:
# gen.H0 is a logical indicating whether data under H0 should be generated, this is 
# not necessary for some results
# m is the number of datasets generated for each of the two hypotheses at every iteration
# t.points captures the number and position of measurement occasions per person
# var.u0 is the intercept variance
# var.u1 is the slope variance
# eff.size is the effect size (sqrt(var.u1) / beta2)
# BFthres is the threshold which a BF needs to exceed in order to be considered of importance
# eta is the desired power level (i.e., the probability of obtaining a BF>BFthres)
# These values are passed on to the data generation function "dat.gen.vec"
# 
# m=1000                # number of datasets per scenario per iteration
# N=30                  # number of individuals
# log=F                 # log linear growth?
# t.points=c(0,1,2,3,4) # vector containing time points
# var.u0=0.0333         # intercept variance
# var.u1=.1             # slope variance
# var.e=.02             # residual variance
# eff.size=.8           # effect size
# cov=0                 # covariance between random effects
# BFthres=3             # threshold which the BF needs to exceed to be considered convincing
# Neff="worst"          # if "worst": Neff=N, if "best", Neff=N*n

dat.gen.vec.hand.b <- function(m=1000, N=30, log=F, t.points=c(0,1,2,3,4), var.u0=0.0333, var.u1=.1, cov=0, var.e=.02, eff.size=.8, BFthres=3, Neff="worst"){
  
  set.seed(123)                                                              # set a seed for reproducibility
  hypotheses <- c("t:treat>0;t:treat=0")                                     # the two competing hypotheses of interest 
  n <- length(t.points)                                                      # number of measurements per person
  ifelse(log==F,                                                             # time variable: linear or log linear
         t <- rep(t.points, N), 
         ifelse(min(t.points)==0,
                t <- rep(log(t.points+1), N),                                # if the first time point is zero, we add 1 to all time points because log(0) is undefined
                t <- rep(log(t.points), N)
         )
  )                                                                          
  id <- rep(seq_len(N), each=n)                                              # create ID variable
  ifelse(Neff=="worst",                                                      # in the worst case: Neff=N, in the best case: Neff=N*n
         b1 <- 1/N,
         b1 <- 1/N*n)                                                 # fraction of the data for prior: b1=1/N, b2=2/N, b3=3/N
  b2 <- 2*b1                                                                  
  b3 <- 3*b1
  treat <- as.numeric(as.character(gl(n=2, k=n, length=N*n, labels=c(0,1)))) # create treatment variable
  dat0 <- data.frame(id, treat, t)                                           # create an empty data frame with all variables except the outcome
  sigma <- matrix(c(var.u0,cov,cov,var.u1),2,2)                              # variance covariance matrix of random effects
  rand.eff <- MASS::mvrnorm(n=N, mu=c(0,0), Sigma=sigma)                     # draw random effects u0 and u1 from multivariate normal distribution
  
  fit0.H0 <- vector("list", m)                                               # empty objects for storing fit, complexity, and the Bayes Factors
  fit1.H0 <- vector("list", m)
  comp0.H0.b1 <- vector("list", m)
  comp0.H0.b2 <- vector("list", m)
  comp0.H0.b3 <- vector("list", m)
  comp1.H0.b1 <- vector("list", m)
  comp1.H0.b2 <- vector("list", m)
  comp1.H0.b3 <- vector("list", m)
  
  fit0.H1 <- vector("list", m)
  fit1.H1 <- vector("list", m)
  comp0.H1.b1 <- vector("list", m)
  comp0.H1.b2 <- vector("list", m)
  comp0.H1.b3 <- vector("list", m)
  comp1.H1.b1 <- vector("list", m)
  comp1.H1.b2 <- vector("list", m)
  comp1.H1.b3 <- vector("list", m)
  
  BFu.H0.b1 <- vector("list", m)
  BFu.H0.b2 <- vector("list", m)
  BFu.H0.b3 <- vector("list", m)
  BFu1.H0.b1 <- vector("list", m)
  BFu1.H0.b2 <- vector("list", m)
  BFu1.H0.b3 <- vector("list", m)
  BFs.H0.b1 <- vector("list", m)
  BFs.H0.b2 <- vector("list", m)
  BFs.H0.b3 <- vector("list", m)
  
  BFu.H1.b1 <- vector("list", m)
  BFu.H1.b2 <- vector("list", m)
  BFu.H1.b3 <- vector("list", m)
  BFu0.H1.b1 <- vector("list", m)
  BFu0.H1.b2 <- vector("list", m)
  BFu0.H1.b3 <- vector("list", m)
  BFs.H1.b1 <- vector("list", m)
  BFs.H1.b2 <- vector("list", m)
  BFs.H1.b3 <- vector("list", m)
  
  beta2.H1 <- eff.size * sqrt(var.u1) # calculate the beta2 based on effect size and slope variance
  
  make.y0 <- function(N){
    cbind(dat0, 0 + rep(rand.eff[,1], each=n) + rep(rand.eff[,2], each=n)*t + rnorm(N*n, 0, sqrt(var.e)))
  }                                                                                                        # make a function for creating outcome variable y under H0
  dat.H0.empty <- rep(list(N),m)                                                                             # make m empty datasets
  dat.H0 <- lapply(dat.H0.empty, make.y0)                                                                    # apply this function to each empty dataset
  dat.H0 <- lapply(dat.H0, setNames, c("id", "treat", "t", "y"))                                             # give names to variables in dataset
  pars.H0 <- list(theta = c(0, sqrt(var.u1), sqrt(var.e)), fixef = c(0, 0, 0))                             # make a list of fixed and random effects under H0
  models.H0 <- lapply(dat.H0, function(x) {
    lmer(y ~ t + t:treat + (t | id), data=x, start = pars.H0, control = lmerControl(calc.derivs = F))})      # fit a multilevel model (mlm) to each dataset
  est.H0 <- as.numeric(lapply(models.H0, function(x) {x@beta[3]}))                                           # extract and store estimates for beta2 from each mlm
  names(est.H0) <- rep("t:treat", m)                                                                         # name the estimates for beta2
  sig.H0 <- lapply(models.H0, function(x) {as.matrix(vcov(x)[3,3])})                                         # extract and store the variance of each estimate
  # finally, calculate and store Bayes Factors for H0 and H1 when H0 is true
  results.H0 <- lapply(1:m, function(i) {
    comp0.H0.b1[[i]] <<- dnorm(0, mean=0, sd=sqrt(sig.H0[[i]]/b1))
    comp0.H0.b2[[i]] <<- dnorm(0, mean=0, sd=sqrt(sig.H0[[i]]/b2))
    comp0.H0.b3[[i]] <<- dnorm(0, mean=0, sd=sqrt(sig.H0[[i]]/b3))
    fit0.H0[[i]] <<- dnorm(0, mean=est.H0[i], sd=sqrt(sig.H0[[i]]))
    comp1.H0.b1[[i]] <<- 1-pnorm(0, mean=0, sd=sqrt(sig.H0[[i]]/b1))
    comp1.H0.b2[[i]] <<- 1-pnorm(0, mean=0, sd=sqrt(sig.H0[[i]]/b2))
    comp1.H0.b3[[i]] <<- 1-pnorm(0, mean=0, sd=sqrt(sig.H0[[i]]/b3))
    fit1.H0[[i]] <<- 1-pnorm(0, mean=est.H0[i], sd=sqrt(sig.H0[[i]]))
    
    BFu.H0.b1[[i]] <<- fit0.H0[[i]]/comp0.H0.b1[[i]]
    BFu.H0.b2[[i]] <<- fit0.H0[[i]]/comp0.H0.b2[[i]]
    BFu.H0.b3[[i]] <<- fit0.H0[[i]]/comp0.H0.b3[[i]]
    BFu1.H0.b1[[i]] <<- fit1.H0[[i]]/comp1.H0.b1[[i]]
    BFu1.H0.b2[[i]] <<- fit1.H0[[i]]/comp1.H0.b2[[i]]
    BFu1.H0.b3[[i]] <<- fit1.H0[[i]]/comp1.H0.b3[[i]]
    
    BFs.H0.b1[[i]] <<- BFu.H0.b1[[i]]/BFu1.H0.b1[[i]]
    BFs.H0.b2[[i]] <<- BFu.H0.b2[[i]]/BFu1.H0.b2[[i]]
    BFs.H0.b3[[i]] <<- BFu.H0.b3[[i]]/BFu1.H0.b3[[i]]
  })
  
  
  # do the same under H1 
  make.y1 <- function(N){
    cbind(dat0, 0 + rep(rand.eff[,1], each=n) + 0*t + beta2.H1*treat*t + rep(rand.eff[,2], each=n)*t + rnorm(N*n, 0, sqrt(var.e)))
  }
  dat.H1.empty <- rep(list(N),m)
  dat.H1 <- lapply(dat.H1.empty, make.y1)
  dat.H1 <- lapply(dat.H1, setNames, c("id", "treat", "t", "y"))
  pars.H1 <- list(theta = c(0, sqrt(var.u1), sqrt(var.e)), fixef = c(0, 0, eff.size * sqrt(var.u1)))
  models.H1 <- lapply(dat.H1, function(x) {
    lmer(y ~ t + t:treat + (t | id), data=x, start = pars.H1, control = lmerControl(calc.derivs = F))})
  est.H1 <- as.numeric(lapply(models.H1, function(x) {x@beta[3]}))
  names(est.H1) <- rep("t:treat", m)
  sig.H1 <- lapply(models.H1, function(x) {as.matrix(vcov(x)[3,3])})
  results.H1 <- lapply(1:m, function(i) {
    comp0.H1.b1[[i]] <<- dnorm(0, mean=0, sd=sqrt(sig.H1[[i]]/b1))
    comp0.H1.b2[[i]] <<- dnorm(0, mean=0, sd=sqrt(sig.H1[[i]]/b2))
    comp0.H1.b3[[i]] <<- dnorm(0, mean=0, sd=sqrt(sig.H1[[i]]/b3))
    fit0.H1[[i]] <<- dnorm(0, mean=est.H1[i], sd=sqrt(sig.H1[[i]]))
    comp1.H1.b1[[i]] <<- 1-pnorm(0, mean=0, sd=sqrt(sig.H1[[i]]/b1))
    comp1.H1.b2[[i]] <<- 1-pnorm(0, mean=0, sd=sqrt(sig.H1[[i]]/b2))
    comp1.H1.b3[[i]] <<- 1-pnorm(0, mean=0, sd=sqrt(sig.H1[[i]]/b3))
    fit1.H1[[i]] <<- 1-pnorm(0, mean=est.H1[i], sd=sqrt(sig.H1[[i]]))
    
    BFu.H1.b1[[i]] <<- fit1.H1[[i]]/comp1.H1.b1[[i]]
    BFu.H1.b2[[i]] <<- fit1.H1[[i]]/comp1.H1.b2[[i]]
    BFu.H1.b3[[i]] <<- fit1.H1[[i]]/comp1.H1.b3[[i]]
    BFu0.H1.b1[[i]] <<- fit0.H1[[i]]/comp0.H1.b1[[i]]
    BFu0.H1.b2[[i]] <<- fit0.H1[[i]]/comp0.H1.b2[[i]]
    BFu0.H1.b3[[i]] <<- fit0.H1[[i]]/comp0.H1.b3[[i]]
    
    BFs.H1.b1[[i]] <<- BFu.H1.b1[[i]]/BFu0.H1.b1[[i]]
    BFs.H1.b2[[i]] <<- BFu.H1.b2[[i]]/BFu0.H1.b2[[i]]
    BFs.H1.b3[[i]] <<- BFu.H1.b3[[i]]/BFu0.H1.b3[[i]]
  })
  
  # generate output with proportions of BFs larger than BFthres 
  output <-  list(Prop_BF0.b1 = length(BFs.H0.b1[BFs.H0.b1>BFthres])/m,
                  Prop_BF0.b2 = length(BFs.H0.b2[BFs.H0.b2>BFthres])/m,
                  Prop_BF0.b3 = length(BFs.H0.b3[BFs.H0.b3>BFthres])/m,
                  Prop_BF1.b1 = length(BFs.H1.b1[BFs.H1.b1>BFthres])/m,
                  Prop_BF1.b2 = length(BFs.H1.b2[BFs.H1.b2>BFthres])/m,
                  Prop_BF1.b3 = length(BFs.H1.b3[BFs.H1.b3>BFthres])/m
  )
  return(output)
}




