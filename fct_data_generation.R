
# var.u0 <- .0333
# var.u1 <- .1
# var.e <- .0262
# t.points <- c(1,2,3,4,5)
# N <- 100
# m <- 100
# eff.size <- .8


dat.gen <- function(m=100, N=72, t.points=c(1,2,3,4,5), var.u0=0, var.u1=.1, var.e=.02, eff.size=.8, BFthres=3){
  library(bain)        # Bayesian estimation
  
  set.seed(123)
  
  hypotheses <- c("t:treat>0;t:treat=0")
  
  n <- length(t.points)
  t <- rep(t.points, N)                                                      # time variable storage 
  id <- rep(seq_len(N), each=n)                                              # create ID variable
  treat <- as.numeric(as.character(gl(n=2, k=n, length=N*n, labels=c(0,1)))) # create treatment variable
  dat0 <- data.frame(id, treat, t)                                           # combine into data frame
  
  models.H0 <- vector("list", m)
  models.H1 <- vector("list", m)
  BFs.H0 <- vector("list", m)
  BFs.H1 <- vector("list", m)
  BFc.H0 <- vector("list", m)
  BFc.H1 <- vector("list", m)
  BFu.H0 <- vector("list", m)
  BFu.H1 <- vector("list", m)
  pmp.c.H0 <- vector("list", m)
  pmp.c.H1 <- vector("list", m)
  
  beta2.H1 <- eff.size * sqrt(var.u1)
  beta2.H0 <- 0
  
  for (i in 1:m) {
    multinorm <- mvrnorm(n=N, mu=c(0,0), matrix(c(var.u0, 0, 0, var.u1), nrow=2, ncol=2)) # draw individual deviation from treatment intercept and slope from a multivariate normal distribution with mean 0.
    u0 <- rep(multinorm[,1], each=n)
    u1 <- rep(multinorm[,2], each=n)
    e <- rnorm(N*n, 0, sqrt(var.e))
    
    yH0 <- u0 + beta2.H0*treat*t + u1*t + e # data-generating mechanism under H0
    yH1 <- u0 + beta2.H1*treat*t + u1*t + e # data-generating mechanism under H1
    dat.H0 <- data.frame(dat0, yH0)
    dat.H1 <- data.frame(dat0, yH1)
    
    models.H0[[i]] <- lmer(yH0 ~ t + t:treat + (t | id), data = dat.H0, control = lmerControl(calc.derivs = F))
    est.H0 <- models.H0[[i]]@beta[3]
    names(est.H0) <- c("t:treat")
    sig.H0 <- list(as.matrix(vcov(models.H0[[i]])[3,3]))
    result.H0 <- bain(est.H0, hypothesis = hypotheses, n = N, Sigma = sig.H0, group_parameters = 1, joint_parameters = 0)
    BFs.H0[[i]] <- result.H0$BFmatrix[2,1]
    
    models.H1[[i]] <- lmer(yH1 ~ t + t:treat + (t | id), data = dat.H1, control = lmerControl(calc.derivs = F))
    est.H1 <- models.H1[[i]]@beta[3]
    names(est.H1) <- c("t:treat")
    sig.H1 <- list(as.matrix(vcov(models.H1[[i]])[3,3]))
    result.H1 <- bain(est.H1, hypothesis = hypotheses, n = N, Sigma = sig.H1, group_parameters = 1, joint_parameters = 0) 
    BFs.H1[[i]] <- result.H1$BFmatrix[1,2]
  }
  return(output = list(Median_BF0 = median(unlist(BFs.H0)), 
                       Median_BF1 = median(unlist(BFs.H1)), 
                       Prop_BF0 = length(BFs.H0[BFs.H0>BFthres])/m,
                       Prop_BF1 = length(BFs.H1[BFs.H1>BFthres])/m,
                       Median_BF_c0 = median(unlist(BFc.H0)),
                       Median_BF_c1 = median(unlist(BFc.H1)),
                       Median_BF_u0 = median(unlist(BFu.H0)),
                       Median_BF_u1 = median(unlist(BFu.H1)),
                       MeanPMP_c0 = mean(unlist(pmp.c.H0)),
                       MeanPMP_c1 = mean(unlist(pmp.c.H1))
                       )
          )
}
