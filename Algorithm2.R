################################################################################
###################### ALGORITHM 2 -- Binary Search ############################
################################################################################

SSD2 <- function(eta=.8, m=1000, log=F, t.points=c(0,1,2,3,4), 
                 var.u0=0.0333, var.u1=.1, cov=0, var.e=.02, 
                 eff.size=.8, BFthres=3, Neff="worst",  fraction=1,
                 sensitivity=F) {
  
  if(eta<0 | eta>1) {stop("'eta' (the desired power level) must be between 0 and 1")}
  if(m%%1!=0 | m<1) {stop("'m' must be a positive integer")}
  if(is.logical(log)==F) {stop("'log' must be either TRUE or FALSE")}
  if(is.logical(sensitivity)==F) {stop("'sensitivity' must be either TRUE or FALSE")}
  if(any(t.points<0)) {stop("all time points must be positive")}
  if(var.u0<0 | var.u1<0 | cov<0 | var.e<0) {stop("all variance components must be between 0 and 1")}
  if(BFthres<0) {stop("'BFthres' must be positive")}
  if(fraction%%1!=0 | fraction<1) {stop("'fraction' must be a positive integer, b=fraction/N")}
  
    
  source("fct_data_generation_hand.R") # call the function for data generation
  start <- Sys.time()                  # measure time it takes to execute function
  
  prop.BF.H0.b1 <- list()
  prop.BF.H0.b2 <- list()
  prop.BF.H0.b3 <- list()
  prop.BF.H1.b1 <- list()
  prop.BF.H1.b2 <- list()
  prop.BF.H1.b3 <- list()
  
  N <- list()
  Nmin <- 30
  Nmax <- 1000
  condition <- FALSE
  j <- 1
  
  suppressMessages({
    
# Without sensitivity analysis -------------------------------------------------
    if(sensitivity==F){ 
        
       while(condition==F){
          
          N[j] <- round((Nmin + Nmax)/2, digits = 0)
          results <- dat.gen.hand(m=m, N=unlist(N[j]), log=log, fraction=fraction, t.points=t.points, var.u0=var.u0, var.u1=var.u1, cov=cov, var.e=var.e, eff.size=eff.size, BFthres=BFthres, Neff=Neff)
          prop.BF.H0[j] <- results$Prop_BF0
          prop.BF.H1[j] <- results$Prop_BF1
          
          ifelse(prop.BF.H0[j]>eta && prop.BF.H1[j]>eta, 
            Nmax <- unlist(N[j]),
            Nmin <- unlist(N[j])
          )
          
          cat("Iteration number", j, "\n", "Power evaluation for a total sample size of N =", unlist(N[j]), "\n")
          
          if(N[j]==Nmin+1 | N[j]==Nmin+2) {condition=TRUE}
          
          j <- j+1
        }
        
        cat("\n", "The recommended sample size to achieve a power of at least", eta, "using b =", b,  "is N =", unlist(N[j-1]), "\n", 
              "Power for H0:", "P(BF01 >", BFthres, "| H0) =", unlist(prop.BF.H0[j-1]), "\n",
              "Power for H1:", "P(BF10 >", BFthres, "| H1) =", unlist(prop.BF.H1[j-1]), "\n", "\n")
          
        print(Sys.time() - start)
        
# With sensitivity analysis ----------------------------------------------------
    } else if(sensitivity==T){ 
      
      for (i in 1:3) {
        
        N <- list()
        Nmin <- 30
        Nmax <- 1000
        condition <- FALSE
        j <- 1
        
        while(condition==F){
          
          N[j] <- round((Nmin + Nmax)/2, digits = 0)
          results <- dat.gen.hand(m=m, N=unlist(N[j]), fraction=i, log=log, t.points=t.points, var.u0=var.u0, var.u1=var.u1, cov=cov, var.e=var.e, eff.size=eff.size, BFthres=BFthres, Neff=Neff)
          prop.BF.H0[j] <- results$Prop_BF0
          prop.BF.H1[j] <- results$Prop_BF1
          
          ifelse(prop.BF.H0[j]>eta && prop.BF.H1[j]>eta, 
                 Nmax <- unlist(N[j]),
                 Nmin <- unlist(N[j]))
          
          cat("Iteration number", j, "\n", "Power evaluation for a total sample size of N =", unlist(N[j]), "with b =", i, "\n") 
          
          if(N[j]==Nmin+1 | N[j]==Nmin+2) {condition=TRUE}
          
          j <- j+1
        }
        
        cat("\n", "The recommended sample size to achieve a power of at least", eta, "using b =", i, "is N =", unlist(N[j-1]), "\n", 
            "Power for H0:", "P(BF01 >", BFthres, "|H0) =", unlist(prop.BF.H0[j-1]), "\n",
            "Power for H1:", "P(BF10 >", BFthres, "|H1) =", unlist(prop.BF.H1[j-1]), "\n", "\n")
        
      }
      print(Sys.time() - start)
    }
    
  })
}






# OLD VERSION

# SSD2.b <- function(eta=.8, m=1000, log=F, t.points=c(0,1,2,3,4), var.u0=0.0333, var.u1=.1, cov=0, var.e=.02, eff.size=.8, BFthres=3, Neff="worst") {
#   
#   if(eta<0 | eta>1) {stop("eta must be between 0 and 1")}
#   if(m%%0==F | m<1) {stop("m must be a positive integer")}
#   if(log!=T | log!=F) {stop("log must be either TRUE or FALSE")}
#   if(any(t.points<0)) {stop("all time points must be positive")}
#   if(var.u0<0 | var.u1<1 | cov<0 | var.e<0) {stop("all variance components must be between 0 and 1")}
#   if(BFthres<0) {stop("BFthres must be positive")}
#   
#   suppressMessages({
#   
#     source("fct_data_generation_hand.R") # call the function for data generation
#     start <- Sys.time()                  # measure time it takes to execute function
#     
#     n.steps <- 20
#     
#     prop.BF.H0.b1 <- list()
#     prop.BF.H0.b2 <- list()
#     prop.BF.H0.b3 <- list()
#     prop.BF.H1.b1 <- list()
#     prop.BF.H1.b2 <- list()
#     prop.BF.H1.b3 <- list()
#     
#     N <- numeric(n.steps)
#     Nmin <- 30
#     Nmax <- 1000
#     condition <- FALSE
#     condition.b1 <- FALSE
#     condition.b2 <- FALSE
#     condition.b3 <- FALSE
#     j <- 1
#     
#     while(condition==F){
#       
#       N[j] <- round((Nmin + Nmax)/2, digits = 0)
#       results <- dat.gen.hand.b(m=m, N=N[j], log=log, t.points=t.points, var.u0=var.u0, var.u1=var.u1, cov=cov, var.e=var.e, eff.size=eff.size, BFthres=BFthres, Neff=Neff)
#       prop.BF.H0.b1[j] <- results$Prop_BF0.b1
#       prop.BF.H0.b2[j] <- results$Prop_BF0.b2
#       prop.BF.H0.b3[j] <- results$Prop_BF0.b3
#       prop.BF.H1.b1[j] <- results$Prop_BF1.b1
#       prop.BF.H1.b2[j] <- results$Prop_BF1.b2
#       prop.BF.H1.b3[j] <- results$Prop_BF1.b3
#       
#       if(prop.BF.H0.b1[j]>eta && prop.BF.H1.b1[j]>eta){
#         N.b1 <- N[j]
#         condition.b1 <- TRUE
#       } else {
#         condition.b1 <- FALSE
#       }
#       
#       if(prop.BF.H0.b2[j]>eta && prop.BF.H1.b2[j]>eta){
#         N.b2 <- N[j]
#         condition.b2 <- TRUE
#       } else {
#         condition.b2 <- FALSE
#       }     
#       
#       if(prop.BF.H0.b3[j]>eta && prop.BF.H1.b3[j]>eta){
#         N.b3 <- N[j]
#         condition.b3 <- TRUE
#       } else {
#         condition.b3 <- FALSE
#       }
#       
#       if(N[j]==Nmin+1) {condition=TRUE}
#       
#       ifelse(condition.b1==TRUE && condition.b2==TRUE && condition.b3==TRUE,
#              Nmax <- N[j], 
#              Nmin <- N[j])
# 
#       
#       print(N[j]) # also print explaining text 
#       j <- j+1
#     }
# 
#     print(Sys.time() - start)
#     
#     if(Neff=="worst"){ 
#            cat("The recommended sample size for b=1/N is", N.b1, "\n", 
#                "The recommended sample size for b=2/N is", N.b2, "\n",
#                "The recommended sample size for b=3/N is", N.b3, "\n")
#     } else {
#            cat("The recommended sample size for b=1/N*n is", N.b1, "\n", 
#                "The recommended sample size for b=2/N*n is", N.b2, "\n",
#                "The recommended sample size for b=3/N*n is", N.b3, "\n")           
#     }
# 
#   })
# }
#   
# 
