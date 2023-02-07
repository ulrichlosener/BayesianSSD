library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(conflicted)
library(dplyr)
library(faux)
conflict_prefer("lmer", "lmerTest")

set.seed(123)

N <- 100 # number of subjects
n <- 5   # number of timepoints

# correlation matrix between tpoints with decreasing r for more distant tpoints
cmat <- abs(row(diag(n)) - col(diag(n)))
x <- 0:5
sequ <- .2 * exp(-.8*x)
#plot(x, sequ, type = "l", xlim = c(0,5), ylim = c(0,1))
for (i in 1:n) {
  cmat[cmat==i] <- sequ[i]
}
diag(cmat) <- 1





# generate data
dat0 <- rnorm_multi(n=N/2, vars=n, mu=seq(1,n), sd=2, r=cmat)
dat0$gender <- 0
dat1 <- rnorm_multi(n=N/2, vars=n, mu=seq(1, (n/4), length=n), sd=2, r=cmat)
dat1$gender <- 1
dat.wide <- rbind(dat0, dat1)
colnames(dat.wide)[1:n] <- seq(1,n)
dat.wide$id <- 1:N

# wide to long
dat <- as.data.frame(dat.wide) %>% 
  pivot_longer(cols = c("1", "2", "3", "4", "5"), names_to = "time", values_to = "y") 
dat$time <- as.numeric(dat$time)
dat$gender <- as.factor(dat$gender)
dat

# descriptives
tapply(dat$y, dat$gender, summary)

# plots

# explore individual trajectories
ggplot(data = dat, aes(x = time, y = y, group = id, color = gender)) + geom_line()

# base
p <- ggplot(data = dat, aes(x = time, y = y, group = id, color = gender)) + geom_jitter(width = .1, height = 0)
# one line for everyone
ggplot(data = dat, aes(x = time, y = y), color = gender) + 
  geom_jitter(aes(color = gender), width = .1, height = 0) + 
  geom_smooth(method = "lm")
# different lines for gender
p  + geom_smooth(group = 0, method = "lm", data = subset(dat, dat$gender == 0)) +
     geom_smooth(group = 0, method = "lm", data = subset(dat, dat$gender == 1))



# fit models

# intercept-only
int.only <- lmer(y ~ 1 + (1 | id), data = dat)
summary(int.only) 
logLik(int.only) 
performance::icc(int.only) # .121 for N=100

p + geom_smooth(aes(group = 1), method = "lm", formula = y ~ 1)

# level 1 predictor: time
lvl1 <- lmer(y ~ time + (1 | id), data = dat)
summary(lvl1)
logLik(lvl1)          # better fit compared to int.only
anova(int.only, lvl1) # likelihood ratio test significant
performance::icc(lvl1)
p + geom_smooth(aes(group = 1), method = "lm")

# add level 2  predictor: gender
lvl2 <- lmer(y ~ time + gender + (1 | id), data = dat)
summary(lvl2)
logLik(lvl2)      # better fit compared to lvl1 predictor only
anova(lvl1, lvl2) # likelihood ratio test significant
performance::icc(lvl2) 

p + stat_smooth(group = 1, method = "lm", formula = y ~ x)
p + stat_smooth(group = 0, method = "lm", data = subset(dat, dat$gender == 0)) +
    stat_smooth(group = 0, method = "lm", data = subset(dat, dat$gender == 1))

# add random slopes for gender
rand.slop <- lmer(y ~ time + gender + (1 + gender | id), data = dat)
summary(rand.slop)
logLik(rand.slop) 
performance::icc(rand.slop)

# add interaction effect between time and gender
inter <- lmer(y ~ time + gender + time:gender + (1 | id), data = dat)
summary(inter)
logLik(inter) 
performance::icc(inter)








dat$time <- as.factor(dat$time)
dat$time <- as.numeric(dat$time)

summary(lm(y ~ time-1, dat))



# variance diagnostics
for(i in 1:n){
print(var(dat$y[dat$time==i & dat$gender==1]))
}



