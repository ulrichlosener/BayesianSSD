#===============================================================================
# Function for Sample Size Determination for Longitudinal Trials - Binary Search
#===============================================================================

# Hypotheses to be tested:
# H0: beta2 = 0
# H1: beta2 > 0
# where beta2 is the coefficient of interaction of time and treatment condition (treatment effect)

# The function "BayeSSD" uses the following arguments:
# eta = the desired power level for each hypothesis
# m = number of datasets created under each hypothesis
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
# sensitivity = should a sensitivity analysis be carried out for different values (1,2,3) for "fraction"?
# hyp = for which hypotheses should the SSD be carried out ("H0", "H1", "both" or "h0", "h1", "b")

# Note 1: this function requires loading the function "getpower" in the global environment.
# Note 2: the packages "MASS" and "lme4" need to be installed (not attached) in order to run this function

#-------------------------------------------------------------------------------

BayeSSD <- function(eta=.8, m=1000, log.grow=F, t.points=c(0,1,2,3,4), 
                    var.u0=0.0333, var.u1=.1, cov=0, var.e=.02, 
                    eff.size=.8, BFthres=3, Neff="worst",  fraction=1,
                    sensitivity=F, seed=NULL, hyp = "both", test="alt") {
  
  # error and warning messages in case of incorrect input
  if(eta<0 | eta>1) {stop("'eta' (the desired power level) must be between 0 and 1.")}
  if(m%%1!=0 | m<1) {stop("'m' must be a positive integer.")}
  if(!is.logical(log.grow)) {stop("'log.grow' must be either TRUE or FALSE.")}
  if(is.logical(sensitivity)==F) {stop("'sensitivity' must be either TRUE or FALSE.")}
  if(any(t.points<0)) {stop("all time points must be positive.")}
  if(var.u0<0 | var.u1<0 | cov<0 | var.e<0) {stop("all variance components must be positive.")}
  if(BFthres<0) {stop("'BFthres' must be positive.")}
  if(fraction%%1!=0 | fraction<1) {stop("'fraction' must be a positive integer, b=fraction/N.")}
  if(m<1000) {warning("Results with less than 1000 generated datasets per iteration can be unreliable and result in power < eta.")}
  if(m>4999) {print("Depending on your machine and available memory, it might take some time to run this function with m > 5000. Enjoy a lunch or a coffee break while waiting for the results!")}
  if(hyp != "both" & hyp != "h1" & hyp != "h0" & hyp != "b" & hyp != "H0" & hyp != "H1") {
    stop("Value for 'hyp' must be either 'both'/'b', 'h0'/'H0', or 'h1'/'H1'." )
  }
  if(test != "alt" & test != "hu" & test != "Hu" & test != "hc" & test != "Hc") {
    stop("Value for 'test' must be either 'alt', 'hu'/'Hu', or 'hc'/'Hc'." )
  }
  if((test=="hu" | test=="Hu") & BFthres>2 & (hyp!="h0") | hyp!="H0") 
    {stop("The value for 'BFthres' is too high in this configuration, BF1u cannot exceed 2.")}
  
  start <- Sys.time()   # measure time it takes to execute function
  source("getbf.R")     # call the function for data generation
  source("getpower.R")  # call the function for power analysis
  if(!is.null(seed)) {set.seed(seed)}  # set user-specified seed for reproducibility
  
  prop.BF.H0 <- list()  # empty objects to store results
  prop.BF.H1 <- list()
  N <- list()
  
  Nmin <- 30            # (initial) minimal sample size 
  Nmax <- 1000          # (initial) maximum sample size
  condition <- FALSE    # condition initially FALSE until power criterion is reached
  j <- 1                # iteration counter
  
  suppressMessages({    # suppress warning messages about singular fit of some MLM models
    
    # Without sensitivity analysis ---------------------------------------------
    if(sensitivity==F){ 
      
      while(condition==F){
        
        N[j] <- round((Nmin + Nmax)/2, digits = 0)  # current N is the mid point between Nmin and Nmax
        # generate data and store BFs
        results <- getpower(m=m, N=unlist(N[j]), log.grow=log.grow, fraction=fraction, 
                            t.points=t.points, var.u0=var.u0, var.u1=var.u1, 
                            cov=cov, var.e=var.e, eff.size=eff.size, 
                            BFthres=BFthres, Neff=Neff, hyp=hyp, test=test)
        
        if(hyp == "both" | hyp == "b"){
          prop.BF.H0[j] <- results$power.H0
          prop.BF.H1[j] <- results$power.H1
          
          # if power>eta in both scenarios, Nmax becomes the current N; if power<eta in both scenarios, Nmin becomes the current N
          ifelse(prop.BF.H0[j]>=eta && prop.BF.H1[j]>=eta, 
                 Nmax <- unlist(N[j]),
                 Nmin <- unlist(N[j])
          )
        } else if(hyp == "h0" | hyp == "H0"){
          prop.BF.H0[j] <- results$power.H0
          
          # if power>eta under H0, Nmax becomes the current N; if power<eta under H0, Nmin becomes the current N
          ifelse(prop.BF.H0[j]>=eta, 
                 Nmax <- unlist(N[j]),
                 Nmin <- unlist(N[j])
          )
        } else if(hyp == "h1" | hyp == "H1"){
          prop.BF.H1[j] <- results$power.H1
          
          # if power>eta under H1, Nmax becomes the current N; if power<eta under H1, Nmin becomes the current N
          ifelse(prop.BF.H1[j]>=eta, 
                 Nmax <- unlist(N[j]),
                 Nmin <- unlist(N[j])
          )
        }
        
        # print text with info about iteration number and current N
        cat("Iteration number", j, "\n", "Power evaluation for a total sample size of N =", unlist(N[j]), "\n")
        
        # if N increases by only 1, condition is met and the algorithm stops
        if(N[j]==Nmin+1 | Nmax==Nmin) {condition=TRUE}
        
        # increase iteration by 1
        j <- j+1
      }
      
      # print results and return power levels
      if(hyp == "both" | hyp == "b"){
        cat("\n", "The recommended sample size to achieve a power of at least", eta, "using b =", fraction, "/ N is N =", unlist(N[j-1]), "\n", 
            "Power for H0:", "P(BF01 >", BFthres, "| H0) =", unlist(prop.BF.H0[j-1]), "\n",
            "Power for H1:", "P(BF10 >", BFthres, "| H1) =", unlist(prop.BF.H1[j-1]), "\n", "\n")
        
        # return if function is assigned to an object but not print
        return(invisible(list(Recommended_N = unlist(N[j-1]),
                         Power_H0 = unlist(prop.BF.H0[j-1]),
                         Power_H1 = unlist(prop.BF.H1[j-1]),
                         b_fraction = paste(fraction, "/ N")
                         )
              )
        )
      } else if(hyp == "h0" | hyp == "H0"){
        
        cat("\n", "The recommended sample size to achieve a power of at least", eta, "using b =", fraction, "/ N is N =", unlist(N[j-1]), "\n", 
            "Power for H0:", "P(BF01 >", BFthres, "| H0) =", unlist(prop.BF.H0[j-1]), "\n", "\n")
        
        return(invisible(list(Recommended_N = unlist(N[j-1]),
                        Power_H0 = unlist(prop.BF.H0[j-1]),
                        b_fraction = paste(fraction, "/ N")
                        )
               )
        )
      } else if(hyp == "h1" | hyp == "H1"){
        
        cat("\n", "The recommended sample size to achieve a power of at least", eta, "using b =", fraction, "/ N is N =", unlist(N[j-1]), "\n", 
            "Power for H1:", "P(BF10 >", BFthres, "| H1) =", unlist(prop.BF.H1[j-1]), "\n", "\n")
      
        return(invisible(list(Recommended_N = unlist(N[j-1]),
                         Power_H1 = unlist(prop.BF.H1[j-1]),
                         b_fraction = paste(fraction, "/ N")
                         )
               )
        )
      }
      
      # print total runtime
      print(Sys.time() - start)
      
    # With sensitivity analysis ------------------------------------------------
    } else if(sensitivity==T){ 
      
      for (i in 1:3) {
        
        N <- list()
        Nmin <- 30
        Nmax <- 1000
        condition <- FALSE
        j <- 1
        
        while(condition==F){
          
          N[j] <- round((Nmin + Nmax)/2, digits = 0)
          results <- getpower(m=m, N=unlist(N[j]), fraction=i, log.grow=log.grow, t.points=t.points, 
                              var.u0=var.u0, var.u1=var.u1, cov=cov, var.e=var.e, 
                              eff.size=eff.size, BFthres=BFthres, Neff=Neff, hyp=hyp)
          
          if(hyp == "both" | hyp == "b"){
            prop.BF.H0[j] <- results$power.H0
            prop.BF.H1[j] <- results$power.H1
            
            # if power>eta in both scenarios, Nmax becomes the current N; if power<eta in both scenarios, Nmin becomes the current N
            ifelse(prop.BF.H0[j]>=eta && prop.BF.H1[j]>=eta, 
                   Nmax <- unlist(N[j]),
                   Nmin <- unlist(N[j])
            )
          } else if(hyp == "h0" | hyp == "H0"){
            prop.BF.H0[j] <- results$power.H0
            
            # same as above but for H0 only
            ifelse(prop.BF.H0[j]>=eta, 
                   Nmax <- unlist(N[j]),
                   Nmin <- unlist(N[j])
            )
          } else if(hyp == "h1" | hyp == "H1"){
            prop.BF.H1[j] <- results$power.H1
            
            # same as above but for H1 only
            ifelse(prop.BF.H1[j]>=eta, 
                   Nmax <- unlist(N[j]),
                   Nmin <- unlist(N[j])
            )
          }
          
          cat("Iteration number", j, "\n", "Power evaluation for a total sample size of N =", unlist(N[j]), "with b =", i, "* b_min", "\n") 
          
          if(N[j]==Nmin+1 | Nmax==Nmin) {condition=TRUE}
          
          j <- j+1
        }
        
        # print results
        if(hyp == "both" | hyp == "b"){
          cat("\n", "The recommended sample size to achieve a power of at least", eta, "using b =", i, "/ N is N =", unlist(N[j-1]), "\n", 
              "Power for H0:", "P(BF01 >", BFthres, "| H0) =", unlist(prop.BF.H0[j-1]), "\n",
              "Power for H1:", "P(BF10 >", BFthres, "| H1) =", unlist(prop.BF.H1[j-1]), "\n", "\n")
          
        } else if(hyp == "h0" | hyp == "H0"){
          
          cat("\n", "The recommended sample size to achieve a power of at least", eta, "using b =", i, "/ N is N =", unlist(N[j-1]), "\n", 
              "Power for H0:", "P(BF01 >", BFthres, "| H0) =", unlist(prop.BF.H0[j-1]), "\n", "\n")
          
        } else if(hyp == "h1" | hyp == "H1"){
          
          cat("\n", "The recommended sample size to achieve a power of at least", eta, "using b =", i, "/ N is N =", unlist(N[j-1]), "\n", 
              "Power for H1:", "P(BF10 >", BFthres, "| H1) =", unlist(prop.BF.H1[j-1]), "\n", "\n")
        }
      }
      
      print(Sys.time() - start)
    }
  })
}

# END OF FUNCTION --------------------------------------------------------------
