#===============================================================================
# Function for power determination for a given sample size
#===============================================================================

# Hypotheses to be tested:
# H0: beta2 = 0
# H1: beta2 > 0
# where beta2 is the coefficient of interaction of time and treatment condition (treatment effect)

# The function "getpower" uses the following arguments:
# m = number of datasets created under each hypothesis
# N = the total sample size (number of subjects)
# t.points = position of the measurement occasions in time
# var.u0 = intercept variance
# var.u1 = slope variance
# cov = covariance between intercept and slope variance
# var.e = error variance
# eff.size = effect size defined as beta/sqrt(var.u1), where beta is the coefficient of interaction
# Bfthres = threshold a Bayes Factor needs to exceed to be considered substantial
# fraction = fraction of information used to specify prior, b = fraction/N
# Neff = if "worst": effective sample size = N, if "best": effective sample size = N*n, 
# where n = number of measurement occasions
# log.grow = indicates whether to use logarithmic (TRUE) or linear growth (FALSE)

# Note that this function requires loading the function "getbf" in the global environment.

#-------------------------------------------------------------------------------

getpower <- function(m=1000, N=72, t.points=c(0,1,2,3,4), 
                         var.u0=0.0333, var.u1=.1, var.e=.0262, cov=0, eff.size=.8, 
                         BFthres=3, fraction=1, Neff="worst", log.grow=F, seed=NULL){
  
  source("getbf.R")  # call function to simulate data and calculate BFs
  
  if(!is.null(seed)) {set.seed(seed)}  # set user-specified seed for reproducibility
  
  Ns <- rep(N, m)  # object to use lapply on with first argument for the function (N)
  suppressMessages({
  # for every m, execute function to simulate data and calculate BFs
  bfs <- lapply(Ns, getbf, t.points=t.points, var.u0=var.u0, var.u1=var.u1, cov=cov, var.e=var.e, 
         eff.size=eff.size, fraction=fraction, Neff=Neff, log.grow=log.grow)
  })
  bf0 <- sapply(bfs, "[[", 1) # extract all BF01
  bf1 <- sapply(bfs, "[[", 2) # extract all BF10
  
  # return list with proportion of BFs>threshold, aka power
  return(list(power.H0 = length(bf0[bf0>BFthres])/m,
              power.H1 = length(bf1[bf1>BFthres])/m))
}

#-------------------------------------------------------------------------------
