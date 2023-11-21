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

bayesianSSD <- function(m=100, t.points=c(1,2,3,4,5), var.u0=0, var.u1=.1, var.e=.02, eff.size=.8, BFthres=3, eta=.8) {
  
  suppressMessages({
  
    library(lme4)       # fit multilevel model
    library(bain)       # Bayesian estimation
    
    source("fct_data_generation_vec_hand.R") # get the function for data generation
    start <- Sys.time()             # measure time it takes to execute function
    
    n.steps <- 20
    
    prop.BF.H0 <- vector("list", n.steps)
    prop.BF.H1 <- vector("list", n.steps)
    
    set.seed(123)      # for reproducibility
    
    N <- numeric(n.steps)
    Nmin <- 30
    Nmax <- 1000
    condition <- F
    j <- 1
    
    
    while(condition==F){
      
      N[j] <- round((Nmin + Nmax)/2, digits = 0)
      
      results <- dat.gen.vec.hand(gen.H0=T, m=m, N=N[j], t.points=t.points, var.u0=var.u0, var.u1=var.u1, var.e=var.e, eff.size=eff.size, BFthres=BFthres, fraction=1)
      prop.BF.H0[j] <- results$Prop_BF0
      prop.BF.H1[j] <- results$Prop_BF1
      
      print(N[j]) # also print explaining text 
      
      if(N[j]==Nmin+1) {condition=T}
      
      ifelse(prop.BF.H0[[j]]>eta & prop.BF.H1[[j]]>eta, Nmax <- N[j], Nmin <- N[j])
      
      j <- j+1
    }

    print(Sys.time() - start) 
    cat("The recommended sample size for fraction=1 is", N[j-1], "\n")
    
    N2 <- numeric(n.steps)
    Nmin2 <- 30
    Nmax2 <- 1000
    condition2 <- F
    k <- 1
    
    while(condition2==F){
      
      ifelse(k==1, N2[k] <- N[j-1], N2[k] <- ceiling((Nmin2 + Nmax2)/2))
  
      results <- dat.gen.vec.hand(gen.H0=T, m=m, N=N2[k], t.points=t.points, var.u0=var.u0, var.u1=var.u1, var.e=var.e, eff.size=eff.size, BFthres=BFthres, fraction=2)
      prop.BF.H0[k] <- results$Prop_BF0
      prop.BF.H1[k] <- results$Prop_BF1
      
      print(N2[k])
      
      if(N2[k]==Nmin2+1) {condition2 <- T}
      
      ifelse(prop.BF.H0[[k]]>eta & prop.BF.H1[[k]]>eta, Nmax2 <- N2[k], Nmin2 <- N2[k])
      
      k <- k+1
    }
    print(Sys.time() - start) 
    cat("The recommended sample size for fraction=2 is", N2[k-1], "\n")
    
    N3 <- numeric(n.steps)
    
    Nmin3 <- 30
    Nmax3 <- 1000
    
    condition3 <- F
    l <- 1
    
    while(condition3==F){
      
      ifelse(l==1, N3[l] <- N[j-1], N3[l] <- ceiling((Nmin3 + Nmax3)/2))
      
      results <- dat.gen.vec.hand(gen.H0=T, m=m, N=N3[l], t.points=t.points, var.u0=var.u0, var.u1=var.u1, var.e=var.e, eff.size=eff.size, BFthres=BFthres, fraction=3)
      prop.BF.H0[l] <- results$Prop_BF0
      prop.BF.H1[l] <- results$Prop_BF1
      
      print(N3[l])
      
      if(N3[l]==Nmin3+1) {condition3 <- T}
      
      ifelse(prop.BF.H0[[l]]>eta & prop.BF.H1[[l]]>eta, Nmax3 <- N3[l], Nmin3 <- N3[l])
      
      l <- l+1
    }
    print(Sys.time() - start) 
    cat("The recommended sample size for fraction=3 is", N3[l-1], "\n")
    
  })
}
  

