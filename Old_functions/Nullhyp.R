library(MASS)        # multinorm - already included in lme4?
library(lme4)


n.steps <- 49  # number of different sample sizes to evaluate
m <- 100      # number of datasets
eff.size <- .8 # effect size (Cohen's d): .2 (small), .5 (medium), .8 (large)

t.points <- c(0,1,2,3,4) # time of measurements
n <- length(t.points)    # number of measurements per subject

sigmasq.u0 <- 1   # variance of individual deviation from treatment intercept 
sigmasq.u1 <- 1   # variance of individual deviation from treatment slope
sigma.u0.u1 <- 0  # covariance between sigmasq.u0 and sigmasq.u1. If positive, then individuals with higher(lower) initial values tend to have a higher (lower) rate of change over time.
sigmasq.e <- 1    # error variance


ss.seq <- rep(NA, n.steps)
models.0 <- rep(list(rep(vector("list", m))),n.steps)
models.1 <- rep(list(rep(vector("list", m))),n.steps)

fit1 <- rep(NA, n.steps)
fit2 <- rep(NA, n.steps)

bf.12 <- rep(NA, m)
bf.21 <- rep(NA, m)

pmp1 <- rep(NA, m)
pmp2 <- rep(NA, m)

medbf.12 <- rep(NA, n.steps)
medbf.21 <- rep(NA, n.steps)
varbf.12 <- rep(NA, n.steps)
varbf.21 <- rep(NA, n.steps)
madbf.12 <- rep(NA, n.steps)
madbf.21 <- rep(NA, n.steps)
meanpmp1 <- rep(NA, n.steps)
meanpmp2 <- rep(NA, n.steps)

for(j in 1:n.steps){
  
  N <- j+1 #*4+16      # number of subjects
  ss.seq[j] <- N
  
  # create data vectors
  y <- rep(NA, N)    # data storage
  t <- rep(t.points, N)
  id <- rep(seq_len(N), each=n)
  treat <- as.numeric(as.character(gl(n=2, k=n, length=N*n, labels=c(0,1))))
  dat <- data.frame(id, treat, t)
  
  # create population parameters for data generation
  beta0 <- 0 # average y at t0 
  beta1 <- 0 # average increase for x=0
  beta2.0 <- 0 # average difference in slopes between conditions under H0
  beta2.1 <- eff.size * sqrt(sigmasq.u1) # average difference in slopes between conditions under H1
  
  for (i in 1:m) {
    #seeds[[i]] <- .Random.seed
    multinorm <- mvrnorm(n=N, mu=c(0,0), matrix(c(sigmasq.u0, sigma.u0.u1, sigma.u0.u1, sigmasq.u1), nrow=2, ncol=2)) # draw individual deviation from treatment intercept and slope from a multivariate normal distribution with mean 0.
    u0 <- rep(multinorm[,1], each=n)
    u1 <- rep(multinorm[,2], each=n)
    e <- rnorm(N*n, 0, sqrt(sigmasq.e))
    
    y.0 <- beta0 + u0 + beta1*t + beta2.0*treat*t + u1*t + e # data-generating mechanism under H0
    #y.1 <- beta0 + u0 + beta1*t + beta2H1*treat*t + u1*t + e # data-generating mechanism under H1
    dat.0 <- data.frame(dat, y.0)
    #dat.1 <- data.frame(dat, yH1)
    
    #inter <- lme(y ~ t + t:treat, random =~ t | id, data = dat, control = lmeControl(opt="optim"))
    models.0[[j]][[i]] <- lmer(y.0 ~ t + t:treat + (t | id), data = dat.0, control = lmerControl(calc.derivs = F))
    #models.1[[j]][[i]] <- lmer(yH1 ~ t + t:treat + (t | id), data = datH1, control = lmerControl(calc.derivs = F))
    
    #if (isSingular(models[[j]][[i]])) {next} # uncomment if singular models are not to be included
    
    est.0 <- models.0[[j]][[i]]@beta[3]
    #est.1 <- modelsH1[[j]][[i]]@beta[3]
    
    sig.0 <- list(as.matrix(vcov(models.0[[j]][[i]])[3,3]))
    #sig.1 <- list(as.matrix(vcov(modelsH1[[j]][[i]])[3,3]))
    
    comp1 <- .5
    comp2 <- .5
    
    fit1 <- pnorm(0, mean=est.0, sd=sqrt(as.numeric(sig.0)))
    fit2 <- 1-fit1
    
    bf.12[i] <- (fit1/comp1) / (fit2/comp2)
    bf.21[i] <- 1/bf.12[i]
    
    pmp1[i] <- bf.12[i]/(bf.12[i]+bf.21[i])
    pmp2[i] <- 1-pmp1[i]
  }
  
  medbf.12[j] <- median(bf.12)
  medbf.21[j] <- median(bf.21)
  meanpmp1[j] <- mean(pmp1)
  meanpmp2[j] <- mean(pmp2)
  varbf.12[j] <- var(bf.12)
  varbf.21[j] <- var(bf.21)
  madbf.12[j] <- median(abs(bf.12-medbf.12[j]))
  madbf.21[j] <- median(abs(bf.21-medbf.21[j]))
  
}

#pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/Plots/Nullhyp_N.pdf",   
#    width = 6, height = 6) 
par(mfrow=c(3,2))
plot(x=ss.seq, y=medbf.12[1:length(ss.seq)], type="l", xlab="N", ylab="median BF12")
plot(x=ss.seq, y=medbf.21[1:length(ss.seq)], type="l", xlab="N", ylab="median BF01")
plot(x=ss.seq, y=meanpmp1[1:length(ss.seq)], type="l", xlab="N", ylab="mean PMP1")
plot(x=ss.seq, y=meanpmp2[1:length(ss.seq)], type="l", xlab="N", ylab="mean PMP2")    
plot(x=ss.seq, y=madbf.12[1:length(ss.seq)], type="l", xlab="N", ylab="variance of BF12")    
plot(x=ss.seq, y=madbf.21[1:length(ss.seq)], type="l", xlab="N", ylab="variance of BF21")    

#dev.off()


