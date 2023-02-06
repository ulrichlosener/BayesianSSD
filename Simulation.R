library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(conflicted)
library(dplyr)
library(faux)
conflict_prefer("lmer", "lmerTest")

set.seed(123)

N <- 100
n <- 5
mat <- matrix(NA, nrow=N, ncol=n)
#colnames(mat) <- paste("t", 1:n, sep="")
colnames(mat) <- seq(1:n)
rownames(mat) <- seq(1:N)
id <- 1:N
betw.var <- c(sample(0:1, N, replace = T))
dat <- cbind(mat, id, betw.var)

# correlation matrix between tpoints with decreasing r for more distant tpoints
cmat <- abs(row(diag(n)) - col(diag(n)))
sequ <- seq(.4, 0, -.1)
for (i in 1:n) {
  cmat[cmat==i] <- sequ[i]
}
diag(cmat) <- 1

# generate data
dat0 <- rnorm_multi(n=N/2, vars=n, mu=seq(1,n), sd=1, r=cmat)
dat0$gender <- 0
dat1 <- rnorm_multi(n=N/2, vars=n, mu=seq(1, (n/2), length=n), sd=1, r=cmat)
dat1$gender <- 1
dat <- rbind(dat0, dat1)
colnames(dat)[1:n] <- seq(1,n)
dat$id <- 1:N

# wide to long
dat2 <- as.data.frame(dat) %>% 
  pivot_longer(cols = c("1", "2", "3", "4", "5"), names_to = "time", values_to = "y") 
dat2$time <- as.numeric(dat2$time)
dat2$gender <- as.factor(dat2$gender)
dat2

# descriptives
tapply(dat2$y, dat2$gender, summary)


# plot
# base
p <- ggplot(data = dat2, aes(x = time, y = y, group = id, color = gender)) 
# different lines for gender
p + geom_point() + 
  geom_smooth(group = 0, method = "lm", data = subset(dat2, dat2$gender == 0)) +
  geom_smooth(group = 0, method = "lm", data = subset(dat2, dat2$gender == 1))
# different plots for gender
p + geom_point() + stat_smooth(group = 0, method = "lm") + 
  facet_grid(dat2$gender) 


# fit models
# intercept-only
int.only <- lmer(y ~ 1 + (1 | id), data = dat2)
summary(int.only) 
logLik(int.only) 
performance::icc(int.only) # .33
p + geom_line()  + geom_smooth(aes(group = 1), method = "lm", formula = y ~ 1)

# level 1 predictor: time
lvl1 <- lmer(y ~ time + (1 | id), data = dat2)
summary(lvl1)
logLik(lvl1) # better fit compared to int.only
performance::icc(lvl1)
p + geom_line() + geom_smooth(aes(group = 1), method = "lm")

# add level 2  predictor: gender
lvl2 <- lmer(y ~ time + gender + (1 | id), data = dat2)
summary(lvl2)
logLik(lvl2) # better fit compared to lvl1 predictor only
performance::icc(lvl2) 
p + geom_jitter() + stat_smooth(group = 0, method = "lm", formula = y ~ time)
p + geom_jitter() + 
  stat_smooth(group = 0, method = "lm", fullrange = T, data = subset(dat2, dat2$gender == 0)) +
  stat_smooth(group = 0, method = "lm", fullrange = T, data = subset(dat2, dat2$gender == 1))

# add random slopes for gender
rand.slop <- lmer(y ~ time + (1 + gender | id), data = dat2)
summary(rand.slop)
logLik(rand.slop) 
performance::icc(rand.slop) 







