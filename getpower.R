# 
# var.u0 <- .0333
# var.u1 <- .1
# var.e <- .0262
# t.points <- c(0,1,2,3,4)
# N <- 100
# m <- 100
# eff.size <- .8
# fraction <- 1
# BFthres <- 3

getpower <- function(m=1000, N=72, t.points=c(0,1,2,3,4), 
                         var.u0=0.0333, var.u1=.1, var.e=.0262, cov=0, eff.size=.8, 
                         BFthres=3, fraction=1, Neff="worst", log.grow=F){
  
  source("getbf.R") #  call function to simulate data and calculate BFs
  
  Ns <- rep(N, m)  #  object to use lapply on with first argument for the function (N)
  suppressMessages({
  #  for every m, execute function to simulate data and calculate BFs
  bfs <- lapply(Ns, getbf, var.u0=var.u0, var.u1=var.u1, var.e=var.e, 
         eff.size=eff.size, fraction=fraction, BFthres=BFthres, t.points=t.points, Neff=Neff, log.grow=log.grow)
  })
  bf0 <- sapply(bfs, "[[", 1) #  extract all BF01
  bf1 <- sapply(bfs, "[[", 2) #  extract all BF10
  
  #  return list with proportion of BFs>threshold, aka power
  return(list(power.H0 <- length(bf0[bf0>BFthres])/m,
              power.H1 <- length(bf1[bf1>BFthres])/m))
}


bm <- bench::mark(
  dat.gen.hand(m=10000),
  getpower(m=10000),
  check=F
)
