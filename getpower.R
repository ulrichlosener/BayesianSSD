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
# hyp = for which hypotheses should the power be determined ("H0", "H1", "both" or "h0", "h1", "b")

# Note that this function requires loading the function "getbf" in the global environment.

#-------------------------------------------------------------------------------

getpower <- function(m=1000, N=72, t.points=c(0,1,2,3,4), 
                         var.u0=0.0333, var.u1=.1, var.e=.0262, cov=0, eff.size=.8, 
                         BFthres=3, fraction=1, Neff="worst", log.grow=F, seed=NULL, hyp = "both",
                         test = "alt"){
  
  source("getbf.R")  # call function to simulate data and calculate BFs
  
  if(!is.null(seed)) {set.seed(seed)}  # set user-specified seed for reproducibility
  
  Ns <- rep(N, m)  # object to use lapply on with first argument for the function (N)
  suppressMessages({
  # for every m, execute function to simulate data and calculate BFs
  bfs <- lapply(Ns, getbf, t.points=t.points, var.u0=var.u0, var.u1=var.u1, cov=cov, var.e=var.e, 
         eff.size=eff.size, fraction=fraction, Neff=Neff, log.grow=log.grow, hyp=hyp)
  })
  
  # comparison of H0 and H1 ----------------------------------------------------
  if((hyp == "both" | hyp == "b") & test == "alt"){
    bf0 <- sapply(bfs, "[[", 1) # extract all BF01
    bf1 <- sapply(bfs, "[[", 2) # extract all BF10
    
    # return list with proportion of BFs>threshold, aka power for both hypotheses
    return(list(power.H0 = length(bf0[bf0>=BFthres])/m,
                power.H1 = length(bf1[bf1>=BFthres])/m))
    
  } else if((hyp == "h0" | hyp == "H0") & test == "alt"){
    bf0 <- sapply(bfs, "[[", 1)
      
    # return power for H0 vs H1 only
    return(list(power.H0 = length(bf0[bf0>=BFthres])/m))
    
  } else if((hyp == "h1" | hyp == "1") & test == "alt"){
    bf1 <- sapply(bfs, "[[", 1)
    
    # return power for H1 vs H0 only
    return(list(power.H1 = length(bf1[bf1>=BFthres])/m))
    
    # comparison of H and Hu ---------------------------------------------------
  } else if((hyp == "both" | hyp == "b") & (test == "hu" | test == "Hu")){
    bf0u <- sapply(bfs, "[[", 3) # extract all BF0u
    bf1u <- sapply(bfs, "[[", 4) # extract all BF1u
    
    # return list with proportion of BFs>threshold, aka power for both hypotheses
    return(list(power.H0 = length(bf0u[bf0u>=BFthres])/m,
                power.H1 = length(bf1u[bf1u>=BFthres])/m)) 
    
  } else if((hyp == "h0" | hyp == "H0") & (test == "hu" | test == "Hu")){
    bf0u <- sapply(bfs, "[[", 2)
    
    # return power for H0 vs Hu only
    return(list(power.H0 = length(bf0u[bf0u>=BFthres])/m))
    
  } else if((hyp == "h1" | hyp == "H1") & (test == "hu" | test == "Hu")){
    bf1u <- sapply(bfs, "[[", 2)
    
    # return power for H1 vs Hu only
    return(list(power.H1 = length(bf1u[bf1u>=BFthres])/m))
    
    # comparison of H and Hc ---------------------------------------------------
  } else if((hyp == "both" | hyp == "b") & (test == "hc" | test == "Hc")){
    bf0c <- sapply(bfs, "[[", 5) # extract all BF0u
    bf1c <- sapply(bfs, "[[", 6) # extract all BF1u
    
    # return list with proportion of BFs>threshold, aka power for both hypotheses
    return(list(power.H0 = length(bf0c[bf0c>=BFthres])/m,
                power.H1 = length(bf1c[bf1c>=BFthres])/m)) 
    
  } else if((hyp == "h0" | hyp == "H0") & (test == "hc" | test == "Hc")){
    bf0c <- sapply(bfs, "[[", 3)
    
    # return power for H0 vs Hc only
    return(list(power.H0 = length(bf0c[bf0c>=BFthres])/m))
    
  } else if((hyp == "h1" | hyp == "H1") & (test == "hc" | test == "Hc")){
    bf1c <- sapply(bfs, "[[", 3)
    
    # return power for H1 vs Hc only
    return(list(power.H1 = length(bf1c[bf1c>=BFthres])/m))
  }
}

# END OF FUNCTION --------------------------------------------------------------
