################################################################################
# Code for Figure 1 and 2 in the paper #########################################
################################################################################

library(ggplot2)
library(latex2exp)
library(tidyr)

# Figure 1 ---------------------------------------------------------------------
# a) Priors 

# set parameters
beta <- seq(-5, 5, length.out=1000)
N <- 100
b1 <- 1/N
b2 <- 2/N
b3 <- 3/N
brob <- max(2/N, 1/sqrt(N))
sigsq <- .1
sd_prior_b1 <- sqrt(sigsq/b1)
sd_prior_b2 <- sqrt(sigsq/b2)
sd_prior_b3 <- sqrt(sigsq/b3)
sd_prior_brob <- sqrt(sigsq/brob)

# create data
density1 <- dnorm(beta, mean=0, sd=sd_prior_b1)
density2 <- dnorm(beta, mean=0, sd=sd_prior_b2)
density3 <- dnorm(beta, mean=0, sd=sd_prior_b3)
density.b.rob <- dnorm(beta, mean=0, sd=sd_prior_brob)

dat_prior <- as.data.frame(cbind(beta, density1, density2, density3, density.b.rob)) # data frame in wide format with different densities

# plot

ggplot(dat_prior, aes(x = beta)) +
  geom_line(aes(y = density1, linetype = "b min"), lwd = 0.8) +
  geom_line(aes(y = density.b.rob, linetype = "b robust"), lwd = 0.8) +
  geom_area(aes(y = density1, x = ifelse(beta > 0, beta, 0)), fill = "#333333", alpha = 0.5) +
  geom_area(aes(y = density.b.rob, x = ifelse(beta > 0, beta, 0)), fill = "#999999", alpha = 0.5) +
  scale_linetype_manual(name = "b fraction", 
                        values = c("b min" = "dashed", "b robust" = "solid"),
                        labels = c("b min", "b robust")) +
  theme_classic() +
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        legend.position = c(0.8, 0.8), legend.key.width = unit(1, "cm")) +
  geom_segment(aes(x = -5, xend = 0, y = max(density1), yend = max(density1)), linetype = "dotted") +
  geom_segment(aes(x = -5, xend = 0, y = max(density.b.rob), yend = max(density.b.rob)), linetype = "dotted") +
  xlab(bquote(beta[2])) +
  scale_y_continuous(limits = c(0, 0.41), 
                     breaks = c(0, round(max(density1), 2), round(max(density.b.rob), 2))) +
  ylab("prior density") 

# b) Posterior
# create data
beta.post <- seq(-1,2, length.out=1000)
post.density <- dnorm(beta.post, mean=.5, sd=sqrt(.1))
dat_posterior <- as.data.frame(cbind(post.density, beta.post))
annotation <- data.frame(pos.x=.78, pos.y=.2, label="0.94")

# plot
ggplot(data=dat_posterior, aes(x=beta.post, y=post.density)) +
  geom_line(lwd=.8) +
  theme_classic() +
  scale_x_continuous(limits = c(-1, 2), breaks = c(-1, 0, 0.5, 1, 2)) +
  scale_y_continuous(limits = c(0, 1.3), breaks = c(0, 0.36, 0.5, 1)) +
  xlab(bquote(beta[2])) +
  ylab("posterior density") +
  geom_segment(aes(x = -1, xend = 0, y = 0.36, yend = 0.36), linetype="dotted") +
  geom_area(aes(y=post.density, x=ifelse(beta.post>0, beta.post, 0)), fill="#999999", alpha=.5) +
  geom_label(data=annotation, aes(x=pos.x, y=pos.y, label=label), size=3)

# Figure 2 ---------------------------------------------------------------------
# create data, calculate BFs
res <- list() # storage for results
set.seed(123) # for reproducibility
for(i in 1:100000){
  res[[i]] <- getbf(N=72, t.points=c(0,1,2,3,4), var.u0=0.0333, var.u1=.1, var.e=.0262, 
                    cov=0, eff.size=.8, fraction=1, Neff="worst", log.grow=F, hyp = "both")
  if (i%%100 == 0) {print(i)} # print every 100th iteration to track progress
}

# a) Sampling distribution of BFs under H0
BF0 <- sapply(res, "[[", 1) # extract BF01s

# plot
p0 <- ggplot() + 
  geom_density(aes(x=BF0), linewidth=1, outline.type = "full", bounds=c(0, 50)) + 
  labs(y="density", x="BF01") +
  theme_classic() 

dpb0 <- ggplot_build(p0)

p0 + geom_area(data=data.frame(x=dpb0$data[[1]]$x[dpb0$data[[1]]$x>3], y=dpb0$data[[1]]$y[dpb0$data[[1]]$x>3]), aes(x=x, y=y), fill="grey") +
  annotate("text", x=c(1.5,10), y=c(.006,.006), label=c(TeX("$BF_{thresh} = 3$"), TeX("$\\eta \\geq .80$")), angle=c(90,0)) +
  geom_vline(xintercept=3, linetype=2, linewidth=1) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_blank())

# b) Sampling distribution of BFs under H1
BF1 <- sapply(res, "[[", 2) # extract BF10s

# plot
p1 <- ggplot() + 
  geom_density(aes(x=BF1), linewidth=1, outline.type = "full", bounds=c(0, 50)) + 
  labs(y="density", x="BF10") +
  theme_classic() 

dpb1 <- ggplot_build(p1)

p1 + geom_area(data=data.frame(x=dpb1$data[[1]]$x[dpb1$data[[1]]$x>3], y=dpb1$data[[1]]$y[dpb1$data[[1]]$x>3]), aes(x=x, y=y), fill="grey") +
  annotate("text", x=c(1.5,10), y=c(.01,.01), label=c(TeX("$BF_{thresh} = 3$"), TeX("$\\eta \\geq .80$")), angle=c(90,0)) +
  geom_vline(xintercept=3, linetype=2, linewidth=1) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_blank())

# END OF FILE ------------------------------------------------------------------
