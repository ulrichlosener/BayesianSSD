
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
# m=1000
# N=30
# log=F
# t.points=c(1,2,3,4,5,6,7,8,9,10)
# var.u0=0.0333
# var.u1=.1
# var.e=.02
# eff.size=.8
# BFthres=3
# fraction=1
# cov=0

dat.gen.vec.hand <- function(m=1000, N=30, log=F, t.points=c(1,2,3,4,5), var.u0=0.0333, var.u1=.1, var.e=.02, eff.size=.8, BFthres=3, fraction=1, cov=0, Neff="worst"){
  
  set.seed(123) # set a seed for reproducibility
  
  hypotheses <- c("t:treat>0;t:treat=0") # these are the two competing hypotheses of interest. 
  
  n <- length(t.points)                                                      # number of measurements per person
  ifelse(log==F, 
         t <- rep(t.points, N), 
         t <- rep(log(t.points), N))                                         # time variable: linear or loglinear
  id <- rep(seq_len(N), each=n)                                              # create ID variable
  ifelse(Neff=="worst",
        b <- fraction/N,
        b <- fraction/N*n)                                                   # fraction of the data for prior
  treat <- as.numeric(as.character(gl(n=2, k=n, length=N*n, labels=c(0,1)))) # create treatment variable
  dat0 <- data.frame(id, treat, t)                                           # create an empty data frame with all variables except the outcome
  sigma <- matrix(c(var.u0,cov,cov,var.u1),2,2)                              # variance covariance matrix of random effects
  rand.eff <- MASS::mvrnorm(n=N, mu=c(0,0), Sigma=sigma)                     # draw random effects u0 and u1 from multivariate normal distribution
  
  
  # Empty objects for the results to be stored in
  fit0.H0 <- vector("list", m)
  fit1.H0 <- vector("list", m)
  comp0.H0 <- vector("list", m)
  comp1.H0 <- vector("list", m)
  
  fit0.H1 <- vector("list", m)
  fit1.H1 <- vector("list", m)
  comp0.H1 <- vector("list", m)
  comp1.H1 <- vector("list", m)
  
  BFs.H0 <- vector("list", m)
  BFs.H1 <- vector("list", m)
  BFc.H0 <- vector("list", m)
  BFc.H1 <- vector("list", m)
  BFu.H0 <- vector("list", m)
  BFu.H1 <- vector("list", m)
  BFu1.H0 <- vector("list", m)
  BFu0.H1 <- vector("list", m)
  # pmp.c.H0 <- vector("list", m)
  # pmp.c.H1 <- vector("list", m)
  
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
    # finally, calculate and store Bayes Factors (BFs) for H0 
    results.H0 <- lapply(1:m, function(i) {
      comp0.H0[[i]] <<- dnorm(0, mean=0, sd=sqrt(sig.H0[[i]]/b))
      fit0.H0[[i]] <<- dnorm(0, mean=est.H0[i], sd=sqrt(sig.H0[[i]]))
      comp1.H0[[i]] <<- pnorm(0, mean=0, sd=sqrt(sig.H0[[i]]/b))
      fit1.H0[[i]] <<- 1-pnorm(0, mean=est.H0[i], sd=sqrt(sig.H0[[i]]))
      
      BFu.H0[[i]] <<- fit0.H0[[i]]/comp0.H0[[i]]
      BFc.H0[[i]] <<- BFu.H0[[i]]
      BFu1.H0[[i]] <<- fit1.H0[[i]]/comp1.H0[[i]]
      BFs.H0[[i]] <<- BFu.H0[[i]]/BFu1.H0[[i]]
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
    lmer(y ~ t + t:treat + (t | id), data=x)})
  est.H1 <- as.numeric(lapply(models.H1, function(x) {x@beta[3]}))
  names(est.H1) <- rep("t:treat", m)
  sig.H1 <- lapply(models.H1, function(x) {as.matrix(vcov(x)[3,3])})
  results.H1 <- lapply(1:m, function(i) {
    comp0.H1[[i]] <<- dnorm(0, mean=0, sd=sqrt(sig.H1[[i]]/b))
    fit0.H1[[i]] <<- dnorm(0, mean=est.H1[i], sd=sqrt(sig.H1[[i]]))
    comp1.H1[[i]] <<- pnorm(0, mean=0, sd=sqrt(sig.H1[[i]]/b))
    fit1.H1[[i]] <<- 1-pnorm(0, mean=est.H1[i], sd=sqrt(sig.H1[[i]]))
    
    BFu.H1[[i]] <<- fit1.H1[[i]]/comp1.H1[[i]]
    BFc.H1[[i]] <<- (fit1.H1[[i]]/comp1.H1[[i]]) / ((1-fit1.H1[[i]])/(1-comp1.H1[[i]]))
    BFu0.H1[[i]] <<- fit0.H1[[i]]/comp0.H1[[i]]
    BFs.H1[[i]] <<- BFu.H1[[i]]/BFu0.H1[[i]]
  })
  
  # generate output 
  output <-  list(Median_BF0 = median(unlist(BFs.H0)), 
                  Median_BF1 = median(unlist(BFs.H1)), 
                  Prop_BF0 = length(BFs.H0[BFs.H0>BFthres])/m,
                  Prop_BF1 = length(BFs.H1[BFs.H1>BFthres])/m,
                  Prop_BFc0 = length(BFc.H0[BFc.H0>BFthres])/m,
                  Prop_BFc1 = length(BFc.H1[BFc.H1>BFthres])/m,
                  Median_BF_c0 = median(unlist(BFc.H0)),
                  Median_BF_c1 = median(unlist(BFc.H1)),
                  Median_BF_u0 = median(unlist(BFu.H0)),
                  Median_BF_u1 = median(unlist(BFu.H1)),
                  m.est.H0 = mean(est.H0),
                  m.est.H1 = mean(est.H1),
                  m.var.H0 = mean(unlist(sig.H1)),
                  m.var.H1 = mean(unlist(sig.H0)),
                  m.comp0_H0 = mean(unlist(comp0.H0)),
                  m.comp1_H1 = mean(unlist(comp1.H1)),
                  m.fit0_H0 = mean(unlist(fit0.H0)),
                  m.fit1_H1 = mean(unlist(fit1.H1)),
                  Median_BF0u = median(unlist(BFu.H0)),
                  Median_BF1u = median(unlist(BFu.H1))
                  # MeanPMP_c0 = mean(unlist(pmp.c.H0)),
                  # MeanPMP_c1 = mean(unlist(pmp.c.H1))
                  )
  return(output)
}




