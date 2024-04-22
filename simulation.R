### RAUDENBUSH TABLES ###   

# After importing RData file

for(i in 1:5){
  assign(paste0("duration_", i+1), res[[i]])
  assign(paste0("prop_BF0_d", i+1), rep(NA, 6))
  assign(paste0("prop_BF1_d", i+1), rep(NA, 6))
}

for(i in 1:6){
  prop_BF0_d2[i] <- duration_2[[i+1]][3]
  prop_BF1_d2[i] <- duration_2[[i+1]][4]
  
  prop_BF0_d3[i] <- duration_3[[i+1]][3]
  prop_BF1_d3[i] <- duration_3[[i+1]][4]
  
  prop_BF0_d4[i] <- duration_4[[i+1]][3]
  prop_BF1_d4[i] <- duration_4[[i+1]][4]
  
  prop_BF0_d5[i] <- duration_5[[i+1]][3]
  prop_BF1_d5[i] <- duration_5[[i+1]][4]
  
  prop_BF0_d6[i] <- duration_6[[i+1]][3]
  prop_BF1_d6[i] <- duration_6[[i+1]][4]
}

## Part I: Neff = N

# Table 1: Effect of study duration and frequency of observation


d <- rep(NA, 7)
f <- rep(NA, 7)
t_m <- vector("list", 7)
res.tab1 <- vector("list", 7)

for(i in 1:7){
d[i] <- 4                                   # study duration
ifelse(i==1, f[i] <- 0.5, f[i] <- i-1)      # frequency
t_m[[i]] <- seq(from=0, to=d[i], by=1/f[i])
}


# increasing frequency while holding duration constant

set.seed(123)
for(i in 5:7){
  suppressMessages({
    start <- Sys.time()
    res.tab1[[i]] <- dat.gen.hand(m=10000, N=100, t.points = t_m[[i]], var.u1 = .001, Neff = "worst")
    print(i)
    print(Sys.time()-start)
  })
}



tab1_d3 <- res.tab1

# save(tab1_d3, file = "tab1_d3.RData")

# plot results
med_BF0_d3 <- rep(NA, 6)
med_BF1_d3 <- rep(NA, 6)
prop_BF0_d3 <- rep(NA, 6)
prop_BF1_d3 <- rep(NA, 6)


for(i in 1:6){
  med_BF0_d3[i] <- tab1_d3[[i+1]][1]
  med_BF1_d3[i] <- tab1_d3[[i+1]][2]
  prop_BF0_d3[i] <- tab1_d3[[i+1]][3]
  prop_BF1_d3[i] <- tab1_d3[[i+1]][4]
}




# Table 2: study duration and sample size
  # holding N constant while increasing D

library(lme4)

res.tab2 <- vector("list", 10)
N <- seq(20, 200, by=20)



for(i in 1:10){
  start <- Sys.time()
  res.tab2[[i]] <- dat.gen.vec.hand(m=1000, N=N[i], t.points = c(0,1,2), var.u1 = .001, Neff = "worst")
  print(i)
  print(Sys.time()-start)
}

  tab2_N140 <- res.tab2

prop_BF0_N140 <- rep(NA, 5)
prop_BF1_N140 <- rep(NA, 5)
med_BF0_N140 <- rep(NA, 5)
med_BF1_N140 <- rep(NA, 5)

for(i in 1:7){
  med_BF0_N140[i] <- tab2_N140[[i]][1]
  med_BF1_N140[i] <- tab2_N140[[i]][2]
  prop_BF0_N140[i] <- tab2_N140[[i]][3]
  prop_BF1_N140[i] <- tab2_N140[[i]][4]
}



# holding study duration constant
res.tab3 <- vector("list", 8)

for(i in 1:8){
  res.tab3[[i]] <- dat.gen.vec.hand.no.t(m=1000, N=100, t.points = seq(0,i+1), var.u1 = .001, Neff = "worst")
  print(i)
}

res.b <- res.tab3

prop_BF0.b <- rep(NA, 10)
prop_BF1.b <- rep(NA, 10)
med_BF0.b <- rep(NA, 10)
med_BF1.b <- rep(NA, 10)
m_est0.b <- rep(NA, 10)
m_est1.b <- rep(NA, 10)
m_var0.b <- rep(NA, 10)
m_var1.b <- rep(NA, 10)
comp_0.H0.b <- rep(NA, 10)
comp_1.H1.b <- rep(NA, 10)
fit_0.H0.b <- rep(NA, 10)
fit_1.H1.b <- rep(NA, 10)
med_BFu0.b <- rep(NA, 10) 
med_BFu1.b <- rep(NA, 10) 

for(i in 1:10){
  med_BF0.b[i] <- res.b[[i]][1]
  med_BF1.b[i] <- res.b[[i]][2]
  prop_BF0.b[i] <- res.b[[i]][3]
  prop_BF1.b[i] <- res.b[[i]][4]
  m_est0.b[i] <- res.b[[i]][11]
  m_est1.b[i] <- res.b[[i]][12]
  m_var0.b[i] <- res.b[[i]][13]
  m_var1.b[i] <- res.b[[i]][14]
  comp_0.H0.b[i] <- res.b[[i]][15]
  comp_1.H1.b[i] <- res.b[[i]][16]
  fit_0.H0.b[i] <- res.b[[i]][17]
  fit_1.H1.b[i] <- res.b[[i]][18]
  med_BFu0.b[i] <- res.b[[i]][19]
  med_BFu1.b[i] <- res.b[[i]][20]
}

plot(unlist(prop_BF0), type="l", ylim=c(-1,1))



plot(seq(3,12), prop_BF0_50c, type="l")
plot(seq(3,12), med_BF0_50c, type="l")

plot(seq(3,12), prop_BF0_150c, type="l")
plot(seq(3,12), med_BF0_150c, type="l")


plot(seq(3,12), m_est0_50c, type="l")
plot(seq(3,12), m_var0_50c, type="l")

par(mfrow=c(1,2))
plot(seq(3,12), comp_0.H0c, type="l")
plot(seq(3,12), fit_0.H0c, type="l")

par(mfrow=c(1,2))
plot(seq(3,12), med_BFu0_150c, type="l")
plot(seq(3,12), med_BF1_150c, type="l")

plot(seq(3,12), fit_0.H0c, type="l")
plot(seq(3,12), fit_1.H1c, type="l")


ratio135 <- vector("list", 10)

for(i in 1:10){
  ratio135[[i]] <- m_est0_135c[[i]]^2/m_var0_135c[[i]]*2
}

par(mfrow=c(1,2))
plot(seq(3,12), ratio135, type="l")
plot(seq(3,12), ratio50, type="l")





#pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.50.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,13),med_BF0_50, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),med_BF1_50, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),prop_BF0_50, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,13),prop_BF1_50, type="l", xlab = "number of measurements", ylab = "power")
#dev.off()

#pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.100.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,13),med_BF0_100, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),med_BF1_100, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),prop_BF0_100, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,13),prop_BF1_100, type="l", xlab = "number of measurements", ylab = "power")
#dev.off()

# from here on, trend reverses for H0
#pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.150.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,13),med_BF0_150, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),med_BF1_150, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),prop_BF0_150, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,13),prop_BF1_150, type="l", xlab = "number of measurements", ylab = "power")
#dev.off()

#pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.200.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,12),med_BF0_200, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,12),med_BF1_200, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,12),prop_BF0_200, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,12),prop_BF1_200, type="l", xlab = "number of measurements", ylab = "power")
#dev.off()

# m=1000
#pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.125.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,13),med_BF0_125, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),med_BF1_125, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),prop_BF0_125, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,13),prop_BF1_125, type="l", xlab = "number of measurements", ylab = "power")
#dev.off()

#pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.135.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,13),med_BF0_135, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),med_BF1_135, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),prop_BF0_135, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,13),prop_BF1_135, type="l", xlab = "number of measurements", ylab = "power")
#dev.off()

# pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.140.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,13),med_BF0_140, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),med_BF1_140, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),prop_BF0_140, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,13),prop_BF1_140, type="l", xlab = "number of measurements", ylab = "power")
# dev.off()



## Part II: Neff = N*n


# Table 3: frequency of observation and sample size
library(lme4)



res.tab3b <- vector("list", 11)

for(i in 1:11){
  res.tab3b[[i]] <- dat.gen.vec.hand(m=1000, N=25, t.points = c(1:(i+2)), var.u1 = .001, Neff = "best")
}


res25 <- res.tab3b

prop_BF0_25 <- rep(NA, 11)
prop_BF1_25 <- rep(NA, 11)
med_BF0_25 <- rep(NA, 11)
med_BF1_25 <- rep(NA, 11)

for(i in 1:11){
  med_BF0_25[i] <- res25[[i]][1]
  med_BF1_25[i] <- res25[[i]][2]
  prop_BF0_25[i] <- res25[[i]][3]
  prop_BF1_25[i] <- res25[[i]][4]
}

# pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.50b.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,13),med_BF0_50, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),med_BF1_50, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),prop_BF0_50, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,13),prop_BF1_50, type="l", xlab = "number of measurements", ylab = "power")
# dev.off()

# pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.100b.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,13),med_BF0_100, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),med_BF1_100, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),prop_BF0_100, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,13),prop_BF1_100, type="l", xlab = "number of measurements", ylab = "power")
# dev.off()

# pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.150b.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,13),med_BF0_150, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),med_BF1_150, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),prop_BF0_150, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,13),prop_BF1_150, type="l", xlab = "number of measurements", ylab = "power")
# dev.off()

# pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.200b.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,13),med_BF0_200, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),med_BF1_200, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),prop_BF0_200, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,13),prop_BF1_200, type="l", xlab = "number of measurements", ylab = "power")
# dev.off()

# pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.250b.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,13),med_BF0_250, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),med_BF1_250, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),prop_BF0_250, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,13),prop_BF1_250, type="l", xlab = "number of measurements", ylab = "power")
# dev.off()


# pdf(file = "C://Users/losen002/OneDrive - Universiteit Utrecht/Desktop/PhD/BayesianSSD/meas.occ.25b.pdf", width = 6, height = 6)
par(mfrow = c(2, 2))
plot(seq(3,13),med_BF0_25, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),med_BF1_25, type="l", xlab = "number of measurements", ylab = "median BF")
plot(seq(3,13),prop_BF0_25, type="l", xlab = "number of measurements", ylab = "power")
plot(seq(3,13),prop_BF1_25, type="l", xlab = "number of measurements", ylab = "power")
# dev.off()



# table 3: frequency of observation and sample size
d <- 4
f <- c(0.5, seq(1:6))
t_m <- vector("list", 7)

res.tab3 <- vector("list", 7)

for(i in 1:7){
  t_m[[i]] <- seq(from=0, to=d, by=1/f[i])
}



for(i in 1:7){
  start <- Sys.time()
  res.tab2[[i]] <- dat.gen.vec.hand(m=10000, N=20, t.points = t_m[[i]], var.u1 = .001, Neff = "worst")
  print(i)
  print(Sys.time()-start)
}



