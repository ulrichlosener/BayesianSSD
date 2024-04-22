# 
# var.u0 <- .0333
# var.u1 <- .1
# var.e <- .0262
# t.points <- c(0,1,2,3,4)
# N <- 100
# m <- 100
# eff.size <- .8
# fraction <- 1
# cov <- 0

dat.gen.hand <- function(m=1000, N=72, t.points=c(0,1,2,3,4), 
                         var.u0=0.0333, var.u1=.1, var.e=.0262, cov=0, eff.size=.8, 
                         BFthres=3, fraction=1, Neff="worst", log.grow=F){
  n <- length(t.points)
  ifelse(Neff=="worst",
         b <- fraction/N,
         b <- fraction/N*n)
  ifelse(log.grow==F, 
         t <- rep(t.points, N), 
         ifelse(min(t.points)==0,
                t <- rep(log(t.points+1), N),                                # if the first timepoint is zero, we add 1 to all timepoints because log(0) is undefined
                t <- rep(log(t.points), N)
         )
  )
  id <- rep(seq_len(N), each=n)                                              # create ID variable
  treat <- as.numeric(as.character(gl(n=2, k=n, length=N*n, labels=c(0,1)))) # create treatment variable
  dat0 <- data.frame(id, treat, t)                                           # combine into data frame
  
  
  models.H0 <- list()
  models.H1 <- list()
  BFs.H0 <- list()
  BFs.H1 <- list()
  BFc.H0 <- list()
  BFc.H1 <- list()
  BFu.H0 <- list()
  BFu.H1 <- list()
  pmp.a.H0 <- list()
  pmp.a.H1 <- list()
  
  beta2.H1 <- eff.size * sqrt(var.u1)
  beta2.H0 <- 0
  
  for (i in 1:m) {
    multinorm <- MASS::mvrnorm(n=2*N, mu=c(0,0), matrix(c(var.u0, 0, 0, var.u1), nrow=2, ncol=2)) # draw individual deviation from treatment intercept and slope from a multivariate normal distribution with mean 0.
    u0.H0 <- rep(multinorm[1:(nrow(multinorm)/2),1], each=n)
    u0.H1 <- rep(multinorm[(nrow(multinorm)/2+1):(nrow(multinorm)),1], each=n)
    
    u1.H0 <- rep(multinorm[1:(nrow(multinorm)/2),2], each=n)
    u1.H1 <- rep(multinorm[(nrow(multinorm)/2+1):(nrow(multinorm)),2], each=n)
    e.H0 <- rnorm(N*n, 0, sqrt(var.e))
    e.H1 <- rnorm(N*n, 0, sqrt(var.e))
    
    y.H0 <- u0.H0 + beta2.H0*treat*t + u1.H0*t + e.H0 # data-generating mechanism under H0
    y.H1 <- u0.H1 + beta2.H1*treat*t + u1.H1*t + e.H1 # data-generating mechanism under H1
    dat.H0 <- data.frame(dat0, y.H0)
    dat.H1 <- data.frame(dat0, y.H1)
    
    models.H0[[i]] <- lmer(y.H0 ~ t + t:treat + (t | id), data = dat.H0, control = lmerControl(calc.derivs = F))
    est.H0 <- models.H0[[i]]@beta[3]
    names(est.H0) <- c("t:treat")
    sig.H0 <- vcov(models.H0[[i]])[3,3]
    
    # calculate fits and complexities under H0
    comp0.H0 <- dnorm(0, mean=0, sd=sqrt(sig.H0/b))
    fit0.H0 <- dnorm(0, mean=est.H0, sd=sqrt(sig.H0))
    comp1.H0 <- 1-pnorm(0, mean=0, sd=sqrt(sig.H0/b))
    fit1.H0 <- 1-pnorm(0, mean=est.H0, sd=sqrt(sig.H0))
    
    # calculate BFs under H0
    BFu.H0[[i]] <- fit0.H0/comp0.H0
    BFc.H0[[i]] <- BFu.H0[[i]]
    BFu1.H0 <- fit1.H0/comp1.H0
    BFs.H0[[i]] <- BFu.H0[[i]]/BFu1.H0
    pmp.a.H0[[i]] <- BFu.H0[[i]]/(BFu.H0[[i]] + BFu1.H0)
    
    models.H1[[i]] <- lmer(y.H1 ~ t + t:treat + (t | id), data = dat.H1, control = lmerControl(calc.derivs = F))
    est.H1 <- models.H1[[i]]@beta[3]
    names(est.H1) <- c("t:treat")
    sig.H1 <- vcov(models.H1[[i]])[3,3]
    
    # calculate fits and complexities under H1
    comp1.H1 <- 1-pnorm(0, mean=0, sd=sqrt(sig.H1/b))
    fit1.H1 <- 1-pnorm(0, mean=est.H1, sd=sqrt(sig.H1))
    comp0.H1 <- dnorm(0, mean=0, sd=sqrt(sig.H1/b))
    fit0.H1 <- dnorm(0, mean=est.H1, sd=sqrt(sig.H1))
    
    # calculate BFs under H1
    BFu.H1[[i]] <- fit1.H1/comp1.H1
    BFc.H1[[i]] <- (fit1.H1/comp1.H1) / ((1-fit1.H1)/(1-comp1.H1))
    BFu0.H1 <- fit0.H1/comp0.H1
    BFs.H1[[i]] <- BFu.H1[[i]]/BFu0.H1
    pmp.a.H1[[i]] <- BFu.H1[[i]]/(BFu.H1[[i]] + BFu0.H1)
    
  }
  return(output = list(Median_BF0 = median(unlist(BFs.H0)), 
                       Median_BF1 = median(unlist(BFs.H1)), 
                       Prop_BF0 = length(BFs.H0[BFs.H0>BFthres])/m,
                       Prop_BF1 = length(BFs.H1[BFs.H1>BFthres])/m,
                       Median_BF_c0 = median(unlist(BFc.H0)),
                       Median_BF_c1 = median(unlist(BFc.H1)),
                       Median_BF_u0 = median(unlist(BFu.H0)),
                       Median_BF_u1 = median(unlist(BFu.H1)),
                       MeanPMP_a0 = mean(unlist(pmp.a.H0)),
                       MeanPMP_a1 = mean(unlist(pmp.a.H1)),
                       BFs.H0 = BFs.H0,
                       BFs.H1 = BFs.H1)
        )
  
}
