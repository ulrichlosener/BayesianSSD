
################################################################################
###### This function executes the SSD using the data generation function #######
################################################################################

# This function works as follows: the function "dat.gen.vec" is executed iteratively
# with a different N at each iteration until the power criterion is met.
# At each iteration, it is evaluated whether the condition is met and if it is not,
# the next iteration starts

# The ingredients for this function are as follows:
# m is the number of datasets generated for each of the two hypotheses at every iteration
# t.points captures the number and position of measurement occasions per person
# var.u0 is the intercept variance
# var.u1 is the slope variance
# eff.size is the effect size (beta2/ sqrt(var.u1))
# BFthres is the threshold which a BF needs to exceed in order to be considered of importance
# eta is the desired power level (i.e., the probability of obtaining a BF>BFthres)
# These values are passed on to the data generation function "dat.gen.vec"

SSD <- function(m=1000, t.points=c(1,2,3,4,5), var.u0=0, var.u1=.1, var.e=.02, eff.size=.8, BFthres=3, eta=.8) {
  
  library(lme4)       # fit multilevel model
  library(bain)       # Bayesian estimation
  
  source("fct_data_generation_vec_hand.R") # get the function for data generation
  start <- Sys.time()                      # measure time it takes to execute function
  
  N <- 30            # initial N
  condition <- FALSE # condition initially false  
  i <- 1             # iteration number
  maxit <- 300
  
  gen.H0 <- F # generate data under H0 as well? (for some plots not necessary)
  
  # make empty objects to store results in 
  medBF.H0 <- vector("list", maxit)
  medBF.H1 <- vector("list", maxit)
  medBFc.H0 <- vector("list", maxit)
  medBFc.H1 <- vector("list", maxit)
  medBFu.H0 <- vector("list", maxit)
  medBFu.H1 <- vector("list", maxit)
  prop.BF.H0 <- vector("list", maxit)
  prop.BFc.H0 <- vector("list", maxit)
  prop.BF.H1 <- vector("list", maxit)
  prop.BFc.H1 <- vector("list", maxit)
  PMP.c0 <- vector("list", maxit)
  PMP.c1 <- vector("list", maxit)
  
  while (condition==F) {  # while the power criterion is not met, do the following

    N <- N+2 # add one person per group in each iteration (i.e., increase N by 2)
    
    # store results of data generation and BFs in an object 
    results <- dat.gen.vec.hand(gen.H0=T, m=m, N=N, t.points=t.points, var.u0=var.u0, var.u1=var.u1, var.e=var.e, eff.size=eff.size, BFthres=BFthres)
    
    if(gen.H0==TRUE){
      medBF.H0[[i]] <- results$Median_BF0
      prop.BF.H0[[i]] <- results$Prop_BF0
      prop.BFc.H0[[i]] <- results$Prop_BFc0
      
      medBFc.H0[[i]] <- results$Median_BF_c0
      medBFu.H0[[i]] <- results$Median_BF_u0
      PMP.c0[[i]] <- results$MeanPMP_c0
    }
    
    medBF.H1[[i]] <- results$Median_BF1
    prop.BF.H1[[i]] <- results$Prop_BF1
    prop.BFc.H1[[i]] <- results$Prop_BFc1
    medBFc.H1[[i]] <- results$Median_BF_c1
    medBFu.H1[[i]] <- results$Median_BF_u1
    PMP.c1[[i]] <- results$MeanPMP_c1
    
    # condition met, i.e., is the proportion of BFs>BFthres at least eta?
    ifelse(prop.BF.H0[[i]]>eta & prop.BF.H1[[i]]>eta, condition <- TRUE, condition <- FALSE)
    print(list(N, prop.BF.H1[[i]]))
    i <- i+1
  }
  
  ifelse(gen.H0==TRUE, 
         measures <- cbind(medBF.H0,medBF.H1, prop.BF.H0, prop.BFc.H0, prop.BF.H1, medBFc.H0, medBFc.H1, medBFu.H0, medBFu.H1, PMP.c0, PMP.c1),
         measures <- cbind(medBF.H1, prop.BF.H1, prop.BFc.H1, medBFc.H1, medBFu.H1, PMP.c1)
  )
  
  return(measures[i-1,])
  
  print(Sys.time() - start)
}
