#===============================================================================
# Simulation study table 2 and table 3 
#===============================================================================

# table 2: frequency and duration, fix N=100 -----------------------------------
# comparison against alternative hypothesis
dur <- 1:8  # duration
freq <- 1:6  # frequency
dat <- cbind(expand.grid(dur, freq), rep(NA, 48), rep(NA, 48))  # space to store results
names(dat) <- c("dur", "freq", "power.H0", "power.H1")  # name columns

set.seed(123)  # set seed for reproducibility outside the loop
for(i in 2:nrow(dat)){
  start <- Sys.time()
  a <- getpower(m=10000, N=100, t.points = seq(from=0, to=dat[i,"dur"], by=1/dat[i,"freq"]), var.u1 = .001)
  dat[i,"power.H0"] <- a$power.H0  # in each iteration, save power for H0
  dat[i,"power.H1"] <- a$power.H1  # in each iteration, save power for H1
  print(Sys.time()-start)
  print(i)
}

# comparison against Hc
dur <- 1:8  # duration
freq <- 1:6  # frequency
dat1c <- cbind(expand.grid(dur, freq), rep(NA, 48), rep(NA, 48))  # space to store results
names(dat1c) <- c("dur", "freq", "power.H0", "power.H1")  # name columns

set.seed(123)  # set seed for reproducibility outside the loop
for(i in 2:nrow(dat1c)){
  start <- Sys.time()
  a <- getpower(m=10000, N=100, t.points = seq(from=0, to=dat1c[i,"dur"], by=1/dat1c[i,"freq"]), var.u1 = .001, test = "hc")
  dat1c[i,"power.H0"] <- a$power.H0  # in each iteration, save power for H0
  dat1c[i,"power.H1"] <- a$power.H1  # in each iteration, save power for H1
  print(Sys.time()-start)
  print(i)
}

# comparison against Hu
dur <- 1:8  # duration
freq <- 1:6  # frequency
dat1u <- cbind(expand.grid(dur, freq), rep(NA, 48), rep(NA, 48))  # space to store results
names(dat1u) <- c("dur", "freq", "power.H0", "power.H1")  # name columns

set.seed(123)  # set seed for reproducibility outside the loop
for(i in 2:nrow(dat1u)){
  start <- Sys.time()
  a <- getpower(BFthres=1.8, m=10000, N=100, t.points = seq(from=0, to=dat1u[i,"dur"], by=1/dat1u[i,"freq"]), var.u1 = .001, test = "hu")
  dat1u[i,"power.H0"] <- a$power.H0  # in each iteration, save power for H0
  dat1u[i,"power.H1"] <- a$power.H1  # in each iteration, save power for H1
  print(Sys.time()-start)
  print(i)
}

# table 3: N and duration, fix f=1 ---------------------------------------------

# comparison against H1
dur2 <- 2:8  # duration
Ns <- seq(20, 200, by=20)  # sample sizes
dat2 <- cbind(expand.grid(dur2, Ns), rep(NA, 70), rep(NA, 70))  # space to store results
names(dat2) <- c("dur", "N", "power.H0", "power.H1")  # name columns

set.seed(123)  # set seed for reproducibility outside the loop
for(i in 1:nrow(dat2)){
  start <- Sys.time()
  b <- getpower(m=10000, N = dat2[i,"N"], t.points = seq(0, dat2[i,"dur"]), var.u1 = .001)
  dat2[i,"power.H0"] <- b$power.H0  # in each iteration, save power for H0
  dat2[i,"power.H1"] <- b$power.H1  # in each iteration, save power for H1
  print(Sys.time()-start)
  print(i)
}

# comparison against Hc
dur2 <- 2:8  # duration
Ns <- seq(20, 200, by=20)  # sample sizes
dat2c <- cbind(expand.grid(dur2, Ns), rep(NA, 70), rep(NA, 70))  # space to store results
names(dat2c) <- c("dur", "N", "power.H0", "power.H1")  # name columns

set.seed(123)  # set seed for reproducibility outside the loop
for(i in 1:nrow(dat2c)){
  start <- Sys.time()
  b <- getpower(m=10000, N = dat2c[i,"N"], t.points = seq(0, dat2c[i,"dur"]), var.u1 = .001, test = "hc")
  dat2c[i,"power.H0"] <- b$power.H0  # in each iteration, save power for H0
  dat2c[i,"power.H1"] <- b$power.H1  # in each iteration, save power for H1
  print(Sys.time()-start)
  print(i)
}

# comparison against Hu
dur2 <- 2:8  # duration
Ns <- seq(20, 200, by=20)  # sample sizes
dat2u <- cbind(expand.grid(dur2, Ns), rep(NA, 70), rep(NA, 70))  # space to store results
names(dat2u) <- c("dur", "N", "power.H0", "power.H1")  # name columns

set.seed(123)  # set seed for reproducibility outside the loop
for(i in 1:nrow(dat2u)){
  start <- Sys.time()
  b <- getpower(m=10000, N = dat2u[i,"N"], t.points = seq(0, dat2u[i,"dur"]), var.u1 = .001, test = "hu")
  dat2u[i,"power.H0"] <- b$power.H0  # in each iteration, save power for H0
  dat2u[i,"power.H1"] <- b$power.H1  # in each iteration, save power for H1
  print(Sys.time()-start)
  print(i)
}
