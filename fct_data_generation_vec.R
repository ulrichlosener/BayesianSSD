
################################################################################
########## This function generates datasets under both hypotheses, #############
###### fits a multilevel model to each one and calculates BFs  and PMPs ########
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

# var.u0 <- .0333
# var.u1 <- .1
# var.e <- .0262
# t.points <- c(1,2,3,4,5)
# N <- 100
# m <- 100
# eff.size <- .8


dat.gen.vec <- function(gen.H0=T, m=10000, N=72, t.points=c(1,2,3,4,5), var.u0=0, var.u1=.1, var.e=.02, eff.size=.8, BFthres=3, fraction=1){
  library(bain) # Bayesian estimation
  library(lme4) # Multilevel models
  
  set.seed(123) # set a seed for reproducibility
  
  hypotheses <- c("t:treat>0;t:treat=0") # these are the two competing hypotheses of interest. 
  
  n <- length(t.points)                                                      # number of measurements per person
  t <- rep(t.points, N)                                                      # time variable storage 
  id <- rep(seq_len(N), each=n)                                              # create ID variable
  treat <- as.numeric(as.character(gl(n=2, k=n, length=N*n, labels=c(0,1)))) # create treatment variable
  dat0 <- data.frame(id, treat, t)                                           # create an empty data frame with all variables except the outcome

  # Empty objects for the results to be stored in
  BFs.H0 <- vector("list", m)
  BFs.H1 <- vector("list", m)
  BFc.H0 <- vector("list", m)
  BFc.H1 <- vector("list", m)
  BFu.H0 <- vector("list", m)
  BFu.H1 <- vector("list", m)
  pmp.c.H0 <- vector("list", m)
  pmp.c.H1 <- vector("list", m)
  
  beta2.H1 <- eff.size * sqrt(var.u1) # calculate the beta2 based on effect size and slope variance

  gen.H0 <- T # generate data under H0?
  
  if (gen.H0==T){
    make.y0 <- function(N){cbind(dat0, rep(rnorm(N, 0, sqrt(var.u1)), each=n)*t + rnorm(N*n, 0, sqrt(var.e)))} # make a function for creating outcome variable y under H0
    dat.H0.empty <- rep(list(N),m)                                                                             # make m empty datasets
    dat.H0 <- lapply(dat.H0.empty, make.y0)                                                                    # apply this function to each empty dataset
    dat.H0 <- lapply(dat.H0, setNames, c("id", "treat", "t", "y"))                                             # give names to variables in dataset
    pars.H0 <- list(theta = c(0, sqrt(var.u1), sqrt(var.e)), fixef = c(0, 0, 0))                               # make a list of fixed and random effects under H0
    models.H0 <- lapply(dat.H0, function(x) {
      lmer(y ~ t + t:treat + (t | id), data=x, start = pars.H0, control = lmerControl(calc.derivs = F))})      # fit a multilevel model (mlm) to each dataset
    est.H0 <- as.numeric(lapply(models.H0, function(x) {x@beta[3]}))                                           # extract and store estimates for beta2 from each mlm
    names(est.H0) <- rep("t:treat", m)                                                                         # name the estimaes for beta2
    sig.H0 <- lapply(models.H0, function(x) {as.matrix(vcov(x)[3,3])})                                         # extract and store the variance of each estimate
    # finally, calculate and store Bayes Factors (BFs) and posterior model probabilities (PMP) for H0
    results.H0 <- lapply(1:m, function(i) {
      result.H0 <- bain(x = est.H0[i], Sigma = list(sig.H0[[i]]), hypothesis = hypotheses, n=N, group_parameters = 1, joint_parameters = 0, fraction=fraction)
      BFs.H0[[i]] <<- result.H0$BFmatrix[2, 1]
      BFc.H0[[i]] <<- result.H0$fit$BF.c[2] 
      BFu.H0[[i]] <<- result.H0$fit$BF.u[2]
      pmp.c.H0[[i]] <<- result.H0$fit$PMPc[2]
    })
  }
  
  # do the same under H1
  make.y1 <- function(N){cbind(dat0, beta2.H1*treat*t + rep(rnorm(N, 0, sqrt(var.u1)), each=n)*t + rnorm(N*n, 0, sqrt(var.e)))}
  dat.H1.empty <- rep(list(N),m)
  dat.H1 <- lapply(dat.H1.empty, make.y1)
  dat.H1 <- lapply(dat.H1, setNames, c("id", "treat", "t", "y"))
  pars.H1 <- list(theta = c(0, sqrt(var.u1), sqrt(var.e)), fixef = c(0, 0, eff.size * sqrt(var.u1)))
  models.H1 <- lapply(dat.H1, function(x) 
    {lmer(y ~ t + t:treat + (t | id), data=x, start = pars.H1, control = lmerControl(calc.derivs = F))})
  est.H1 <- as.numeric(lapply(models.H1, function(x) {x@beta[3]}))
  names(est.H1) <- rep("t:treat", m)
  sig.H1 <- lapply(models.H1, function(x) {as.matrix(vcov(x)[3,3])})
  results.H1 <- lapply(1:m, function(i) {
    result.H1 <- bain(x = est.H1[i], Sigma = list(sig.H1[[i]]), hypothesis = hypotheses, n=N, group_parameters = 1, joint_parameters = 0, fraction=fraction)
    BFs.H1[[i]] <<- result.H1$BFmatrix[1, 2]
    BFc.H1[[i]] <<- result.H1$fit$BF.c[1]
    BFu.H1[[i]] <<- result.H1$fit$BF.u[1]
    pmp.c.H1[[i]] <<- result.H1$fit$PMPc[1]
    
  })
  
  # generate output if data under H0 was generated
  ifelse(gen.H0==T, output <-  list(Median_BF0 = median(unlist(BFs.H0)), 
                              Median_BF1 = median(unlist(BFs.H1)), 
                              Prop_BF0 = length(BFs.H0[BFs.H0>BFthres])/m,
                              Prop_BF1 = length(BFs.H1[BFs.H1>BFthres])/m,
                              Prop_BFc0 = length(BFc.H0[BFc.H0>BFthres])/m,
                              Prop_BFc1 = length(BFc.H1[BFc.H1>BFthres])/m,
                              Median_BF_c0 = median(unlist(BFc.H0)),
                              Median_BF_c1 = median(unlist(BFc.H1)),
                              Median_BF_u0 = median(unlist(BFu.H0)),
                              Median_BF_u1 = median(unlist(BFu.H1)),
                              MeanPMP_c0 = mean(unlist(pmp.c.H0)),
                              MeanPMP_c1 = mean(unlist(pmp.c.H1))
                    ),
  # generate output when data under H0 was NOT generated
                    output <- list(Median_BF1 = median(unlist(BFs.H1)), 
                              Prop_BF1 = length(BFs.H1[BFs.H1>BFthres])/m,
                              Prop_BFc1 = length(BFc.H1[BFc.H1>BFthres])/m,
                              Median_BF_c1 = median(unlist(BFc.H1)),
                              Median_BF_u1 = median(unlist(BFu.H1)),
                              MeanPMP_c1 = mean(unlist(pmp.c.H1))
                )
        )
  return(output)
}




