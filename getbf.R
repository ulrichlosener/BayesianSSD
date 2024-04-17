

getbf <- function(N, var.u0, var.u1, var.e, eff.size, fraction, BFthres, t.points, Neff, log.grow){
  
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
  beta2.H1 <- eff.size * sqrt(var.u1)
  
  multinorm <- MASS::mvrnorm(n=2*N, mu=c(0,0), matrix(c(var.u0, 0, 0, var.u1), nrow=2, ncol=2))
  u0.H0 <- rep(multinorm[1:(nrow(multinorm)/2),1], each=n)
  u0.H1 <- rep(multinorm[(nrow(multinorm)/2+1):(nrow(multinorm)),1], each=n)
  u1.H0 <- rep(multinorm[1:(nrow(multinorm)/2),2], each=n)
  u1.H1 <- rep(multinorm[(nrow(multinorm)/2+1):(nrow(multinorm)),2], each=n)
  e.H0 <- rnorm(N*n, 0, sqrt(var.e))
  e.H1 <- rnorm(N*n, 0, sqrt(var.e))
  y.H0 <- u0.H0 + 0*treat*t + u1.H0*t + e.H0
  y.H1 <- u0.H1 + beta2.H1*treat*t + u1.H1*t + e.H1
  dat.H0 <- data.frame(dat0, y.H0)
  dat.H1 <- data.frame(dat0, y.H1)
  
  models.H0 <- lmer(y.H0 ~ t + t:treat + (t | id), data = dat.H0, control = lmerControl(calc.derivs = F))
  est.H0 <- models.H0@beta[3]
  sig.H0 <- vcov(models.H0)[3,3]
  
  # calculate fits and complexities under H0
  comp0.H0 <- dnorm(0, mean=0, sd=sqrt(sig.H0/b))
  fit0.H0 <- dnorm(0, mean=est.H0, sd=sqrt(sig.H0))
  comp1.H0 <- 1-pnorm(0, mean=0, sd=sqrt(sig.H0/b))
  fit1.H0 <- 1-pnorm(0, mean=est.H0, sd=sqrt(sig.H0))
  
  # calculate BFs under H0
  BFu.H0 <- fit0.H0/comp0.H0
  BFc.H0 <- BFu.H0
  BFu1.H0 <- fit1.H0/comp1.H0
  BFs.H0 <- BFu.H0/BFu1.H0
  pmp.a.H0 <- BFu.H0/(BFu.H0 + BFu1.H0)
  
  models.H1 <- lmer(y.H1 ~ t + t:treat + (t | id), data = dat.H1, control = lmerControl(calc.derivs = F))
  est.H1 <- models.H1@beta[3]
  names(est.H1) <- c("t:treat")
  sig.H1 <- vcov(models.H1)[3,3]
  
  # calculate fits and complexities under H1
  comp1.H1 <- 1-pnorm(0, mean=0, sd=sqrt(sig.H1/b))
  fit1.H1 <- 1-pnorm(0, mean=est.H1, sd=sqrt(sig.H1))
  comp0.H1 <- dnorm(0, mean=0, sd=sqrt(sig.H1/b))
  fit0.H1 <- dnorm(0, mean=est.H1, sd=sqrt(sig.H1))
  
  # calculate BFs under H1
  BFu.H1 <- fit1.H1/comp1.H1
  BFc.H1 <- (fit1.H1/comp1.H1) / ((1-fit1.H1)/(1-comp1.H1))
  BFu0.H1 <- fit0.H1/comp0.H1
  BFs.H1 <- BFu.H1/BFu0.H1
  pmp.a.H1 <- BFu.H1/(BFu.H1 + BFu0.H1)
  
  return(output = list(BF01 = BFs.H0,
                       BF10 = BFs.H1))
}

