################################################################################
###################### ALGORITHM 2 -- Binary Search ############################
################################################################################

library(lme4)        # fit multilevel model

SSD2 <- function(eta=.8, m=1000, N=30, log=F, t.points=c(0,1,2,3,4), var.u0=0.0333, var.u1=.1, cov=0, var.e=.02, eff.size=.8, BFthres=3, fraction=1, cov=0, Neff="worst") {
  
  suppressMessages({
  
    library(lme4)                            # to fit multilevel model to simulated data
    
    source("fct_data_generation_vec_hand.R") # call the function for data generation
    start <- Sys.time()                      # measure time it takes to execute function
    
    n.steps <- 20
    
    prop.BF.H0 <- vector("list", n.steps)
    prop.BF.H1 <- vector("list", n.steps)
    
    set.seed(123)                            # for reproducibility
    
    N <- numeric(n.steps)
    Nmin <- 30
    Nmax <- 1000
    condition <- FALSE
    condition.b1 <- FALSE
    condition.b2 <- FALSE
    condition.b3 <- FALSE
    j <- 1
    
    while(condition==F){
      
      N[j] <- round((Nmin + Nmax)/2, digits = 0)
      results <- dat.gen.vec.hand.b(m=m, N=N[j], log=log, t.points=t.points, var.u0=var.u0, var.u1=var.u1, cov=cov, var.e=var.e, eff.size=eff.size, BFthres=BFthres, Neff=Neff)
      prop.BF.H0.b1[j] <- results$Prop_BF0.b1
      prop.BF.H0.b2[j] <- results$Prop_BF0.b2
      prop.BF.H0.b3[j] <- results$Prop_BF0.b3
      prop.BF.H1.b1[j] <- results$Prop_BF1.b1
      prop.BF.H1.b2[j] <- results$Prop_BF1.b2
      prop.BF.H1.b3[j] <- results$Prop_BF1.b3
      
      if(prop.BF.H0.b1[j]>eta && prop.BF.H1.b1[j]>eta){
        N.b1 <- N[j]
        condition.b1 <- TRUE
      }
      
      if(prop.BF.H0.b2[j]>eta && prop.BF.H1.b2[j]>eta){
        N.b2 <- N[j]
        condition.b2 <- TRUE
      }      
      
      if(prop.BF.H0.b3[j]>eta && prop.BF.H1.b3[j]>eta){
        N.b3 <- N[j]
        condition.b3 <- TRUE
      }
      
      if(N[j]==Nmin+1) {condition=T}
      
      ifelse(prop.BF.H0.b1[[j]]>eta & prop.BF.H1.b1[[j]]>eta, Nmax <- N[j], Nmin <- N[j])
      
      ifelse(condition.b1==TRUE && condition.b2==TRUE && condition.b3==TRUE,
             Nmax <- N[j], 
             Nmin <- N[j])

      
      print(N[j]) # also print explaining text 
      j <- j+1
    }

    print(Sys.time() - start) 
    cat("The recommended sample size for b=1/N is", N.b1, "\n", 
        "The recommended sample size for b=2/N is", N.b2, "\n",
        "The recommended sample size for b=3/N is", N.b3, "\n")
    

  })
}
  

