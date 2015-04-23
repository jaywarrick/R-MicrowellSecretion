rm(list=ls())
n <- 300
beta <- 4
bg <- 0# rnorm(n, 0.05, 0.01)
alpha <- 3
adjust <- 1
beadsAvg <- 10
beads <- rpois(n,beadsAvg)
virus1 <- rpois(beads,1*beta/beadsAvg) + bg*beads
virus3 <- rpois(beads,3*beta/beadsAvg) + bg*beads
virus5 <- rpois(beads,5*beta/beadsAvg) + bg*beads
virus10 <- rpois(beads,10*beta/beadsAvg) + bg*beads

plot(density(beads))

plot(density(virus10/beads*alpha), ylim=c(0,3), xlim=c(0,4))
lines(density(virus5/beads*alpha, adjust=adjust))
lines(density(virus3/beads*alpha, adjust=adjust))
lines(density(virus1/beads*alpha, adjust=adjust))
lines(density(beads*alpha, adjust=adjust))

plot(density(virus10), ylim=c(0,3), xlim=c(0,4))
lines(density(virus5, adjust=adjust))
lines(density(virus3, adjust=adjust))
lines(density(virus1, adjust=adjust))
lines(density(beads*alpha, adjust=adjust))


duh <- hist(rpois(100000000,1/5000))
duh
