### RAUDENBUSH TABLES ###   


## Part I: Neff = N


# Table 3: frequency of observation and sample size
library(lme4)



res.tab3 <- vector("list", 11)

for(i in 1:11){
  res.tab3[[i]] <- dat.gen.vec.hand(m=1000, N=150, t.points = c(1:(i+2)), var.u1 = .001, Neff = "worst")
}

res150 <- res.tab3

prop_BF0_150 <- rep(NA, 11)
prop_BF1_150 <- rep(NA, 11)
med_BF0_150 <- rep(NA, 11)
med_BF1_150 <- rep(NA, 11)
m_est0_150 <- m_est1_150 <- rep(NA, 11)
m_var0_150 <- rep(NA, 11)
m_var1_150 <- rep(NA, 11)
comp_0.H0 <- rep(NA, 11)
comp_1.H1 <- rep(NA, 11)
fit_0.H0 <- rep(NA, 11)
fit_1.H1 <- rep(NA, 11)
med_BFu0_150 <- rep(NA, 11) 
med_BFu1_150 <- rep(NA, 11) 

for(i in 1:11){
  med_BF0_150[i] <- res150[[i]][1]
  med_BF1_150[i] <- res150[[i]][2]
  prop_BF0_150[i] <- res150[[i]][3]
  prop_BF1_150[i] <- res150[[i]][4]
  m_est0_150[i] <- res150[[i]][11]
  m_est1_150[i] <- res150[[i]][12]
  m_var0_150[i] <- res150[[i]][13]
  m_var1_150[i] <- res150[[i]][14]
  comp_0.H0[i] <- res150[[i]][15]
  comp_1.H1[i] <- res150[[i]][16]
  fit_0.H0[i] <- res150[[i]][17]
  fit_1.H1[i] <- res150[[i]][18]
  med_BFu0_150[i] <- res150[[i]][19]
  med_BFu1_150[i] <- res150[[i]][20]
}



# holding study duration constant
res.tab3c <- vector("list", 10)
meas.occs <- vector("list", 10)

for(i in 5:10){
  meas.occs[[i]] <-  seq(1,3,by=1/i)
  res.tab3c[[i]] <- dat.gen.vec.hand(m=1000, N=150, t.points = seq(1,5,by=1/i), var.u1 = .001, Neff = "worst")
  print(i)
}

res150c <- res.tab3c

prop_BF0_150c <- rep(NA, 10)
prop_BF1_150c <- rep(NA, 10)
med_BF0_150c <- rep(NA, 10)
med_BF1_150c <- rep(NA, 10)
m_est0_150c <- m_est1_150c <- rep(NA, 10)
m_var0_150c <- rep(NA, 10)
m_var1_150c <- rep(NA, 10)
comp_0.H0c <- rep(NA, 10)
comp_1.H1c <- rep(NA, 10)
fit_0.H0c <- rep(NA, 10)
fit_1.H1c <- rep(NA, 10)
med_BFu0_150c <- rep(NA, 10) 
med_BFu1_150c <- rep(NA, 10) 

for(i in 1:10){
  med_BF0_150c[i] <- res150c[[i]][1]
  med_BF1_150c[i] <- res150c[[i]][2]
  prop_BF0_150c[i] <- res150c[[i]][3]
  prop_BF1_150c[i] <- res150c[[i]][4]
  m_est0_150c[i] <- res150c[[i]][11]
  m_est1_150c[i] <- res150c[[i]][12]
  m_var0_150c[i] <- res150c[[i]][13]
  m_var1_150c[i] <- res150c[[i]][14]
  comp_0.H0c[i] <- res150c[[i]][15]
  comp_1.H1c[i] <- res150c[[i]][16]
  fit_0.H0c[i] <- res150c[[i]][17]
  fit_1.H1c[i] <- res150c[[i]][18]
  med_BFu0_150c[i] <- res150c[[i]][19]
  med_BFu1_150c[i] <- res150c[[i]][20]
}



  

plot(seq(3,12), prop_BF0_150c, type="l")
plot(seq(3,12), m_est0_150c, type="l")
plot(seq(3,12), m_var0_150c, type="l")

par(mfrow=c(1,2))
plot(seq(3,12), comp_0.H0c, type="l")
plot(seq(3,12), fit_0.H0c, type="l")

par(mfrow=c(1,2))
plot(seq(3,12), med_BFu0_150c, type="l")
plot(seq(3,12), med_BF1_150c, type="l")

plot(seq(3,12), fit_0.H0c, type="l")
plot(seq(3,12), fit_1.H1c, type="l")


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






