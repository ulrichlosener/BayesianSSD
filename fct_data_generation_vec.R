### Vectorized dat generation

# 
# var.u0 <- .0333
# var.u1 <- .1
# var.e <- .0262
# t.points <- c(1,2,3,4,5)
# N <- 100
# m <- 100
# eff.size <- .8


dat.gen.vec <- function(m=100, N=72, t.points=c(1,2,3,4,5), var.u0=0, var.u1=.1, var.e=.02, eff.size=.8, BFthres=3){
  library(bain)        # Bayesian estimation
  
  set.seed(123)
  
  hypotheses <- c("t:treat>0;t:treat=0")
  
  n <- length(t.points)
  t <- rep(t.points, N)                                                      # time variable storage 
  id <- rep(seq_len(N), each=n)                                              # create ID variable
  treat <- as.numeric(as.character(gl(n=2, k=n, length=N*n, labels=c(0,1)))) # create treatment variable
  dat0 <- data.frame(id, treat, t)

  BFs.H0 <- vector("list", m)
  BFs.H1 <- vector("list", m)
  BFc.H0 <- vector("list", m)
  BFc.H1 <- vector("list", m)
  BFu.H0 <- vector("list", m)
  BFu.H1 <- vector("list", m)
  pmp.c.H0 <- vector("list", m)
  pmp.c.H1 <- vector("list", m)
  
  beta2.H1 <- eff.size * sqrt(var.u1)
  
  make.y0 <- function(N){cbind(dat0, rep(rnorm(N, 0, sqrt(var.u1)), each=n)*t + rnorm(N*n, 0, sqrt(var.e)))}
  make.y1 <- function(N){beta2.H1*treat*t + rep(rnorm(N, 0, sqrt(var.u1)), each=n)*t + rnorm(N*n, 0, sqrt(var.e))}
  
  dat.H0.empty <- rep(list(N),m)
  dat.H1.empty <- rep(list(N),m)
  
  dat.H0 <- lapply(dat.H0.empty, function(N){cbind(dat0, rep(rnorm(N, 0, sqrt(var.u1)), each=n)*t + rnorm(N*n, 0, sqrt(var.e)))})
  dat.H1 <- lapply(dat.H1.empty, function(N){cbind(dat0, beta2.H1*treat*t + rep(rnorm(N, 0, sqrt(var.u1)), each=n)*t + rnorm(N*n, 0, sqrt(var.e)))})
  
  dat.H0 <- lapply(dat.H0, setNames, c("id", "treat", "t", "y"))
  dat.H1 <- lapply(dat.H1, setNames, c("id", "treat", "t", "y"))
  
  pars.H0 <- list(theta = c(0, sqrt(var.u1), sqrt(var.e)), fixef = c(0, 0, 0))
  pars.H1 <- list(theta = c(0, sqrt(var.u1), sqrt(var.e)), fixef = c(0, 0, eff.size * sqrt(var.u1)))
  
  models.H0 <- lapply(dat.H0, function(x) {lmer(y ~ t + t:treat + (t | id), data=x, start = pars.H0, control = lmerControl(calc.derivs = F))})
  models.H1 <- lapply(dat.H1, function(x) {lmer(y ~ t + t:treat + (t | id), data=x, start = pars.H1, control = lmerControl(calc.derivs = F))})
  
  est.H0 <- as.numeric(lapply(models.H0, function(x) {x@beta[3]}))
  est.H1 <- as.numeric(lapply(models.H1, function(x) {x@beta[3]}))
  
  names(est.H0) <- rep("t:treat", m)
  names(est.H1) <- rep("t:treat", m)
  
  sig.H0 <- lapply(models.H0, function(x) {as.matrix(vcov(x)[3,3])})
  sig.H1 <- lapply(models.H1, function(x) {as.matrix(vcov(x)[3,3])})
  
  results.H0 <- lapply(1:m, function(i) {
    result.H0 <- bain(x = est.H0[i], Sigma = list(sig.H0[[i]]), hypothesis = hypotheses, n=N, group_parameters = 1, joint_parameters = 0)
    BFs.H0[[i]] <<- result.H0$BFmatrix[2, 1]
    BFc.H0[[i]] <<- result.H0$fit$BF.c 
    BFu.H0[[i]] <<- result.H0$fit$BF.u
    pmp.c.H0[[i]] <<- result.H0$fit$PMPc[2]
  })
  
  results.H1 <- lapply(1:m, function(i) {
    result.H1 <- bain(x = est.H1[i], Sigma = list(sig.H1[[i]]), hypothesis = hypotheses, n=N, group_parameters = 1, joint_parameters = 0)
    BFs.H1[[i]] <<- result.H1$BFmatrix[1, 2]
    BFc.H1[[i]] <<- result.H1$fit$BF.c 
    BFu.H1[[i]] <<- result.H1$fit$BF.u 
    pmp.c.H1[[i]] <<- result.H1$fit$PMPc[1]
    
  })

  return(output = list(Median_BF0 = median(unlist(BFs.H0)), 
                Median_BF1 = median(unlist(BFs.H1)), 
                Prop_BF0 = length(BFs.H0[BFs.H0>BFthres])/m,
                Prop_BF1 = length(BFs.H1[BFs.H1>BFthres])/m,
                Median_BF_c0 = median(unlist(BFc.H0)),
                Median_BF_c1 = median(unlist(BFc.H1)),
                Median_BF_u0 = median(unlist(BFu.H0)),
                Median_BF_u1 = median(unlist(BFu.H1)),
                PMP_c0 = mean(unlist(pmp.c.H0)),
                PMP_c1 = mean(unlist(pmp.c.H1))
                )
        )
}




