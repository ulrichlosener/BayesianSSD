################################################################################
###################### Power determination for fixed N #########################
################################################################################

PD <- function(m=1000, N=72, t.points=c(0,1,2,3,4), 
               var.u0=0.0333, var.u1=.1, var.e=.0262, cov=0, eff.size=.8, 
               BFthres=3, fraction=1, Neff="worst", log=F,
               dropout=F, omega=.5, gamma=1){
  
  if(m%%1!=0 | m<1) {stop("'m' must be a positive integer")}
  if(is.logical(log)==F) {stop("'log' must be either TRUE or FALSE")}
  if(is.logical(sensitivity)==F) {stop("'sensitivity' must be either TRUE or FALSE")}
  if(any(t.points<0)) {stop("all time points must be positive")}
  if(var.u0<0 | var.u1<0 | cov<0 | var.e<0) {stop("all variance components must be between 0 and 1")}
  if(BFthres<0) {stop("'BFthres' must be positive")}
  if(fraction%%1!=0 | fraction<1) {stop("'fraction' must be a positive integer, b=fraction/N")}
  
  source("fct_data_generation_hand_mis.R") # call the function for data generation
  
  suppressMessages({
  suppressWarnings({
    res <- dat.gen.hand.mis(m=m, N=N, t.points=t.points, 
                            var.u0=var.u0, var.u1=var.u1, var.e=var.e, cov=cov, 
                            eff.size=eff.size, BFthres=BFthres, fraction=fraction, 
                            Neff=Neff, log=log, dropout=dropout, omega=omega, gamma=gamma)
    
    cat("\n", "The statistical power achieved by N =", N, "with b =", fraction, "/", N, "is as follows: \n", 
        "Power for H0:", "P(BF01 >", BFthres, "| H0) =", unlist(res$Prop_BF0), "\n",
        "Power for H1:", "P(BF10 >", BFthres, "| H1) =", unlist(res$Prop_BF1), "\n", "\n")

  })
  })
}
