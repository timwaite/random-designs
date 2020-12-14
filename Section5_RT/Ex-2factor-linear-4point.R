#
#  Supporting material for "Minimax efficient random experimental designs with application to 
#	model-robust design for prediction"
# by Waite, T.W and Woods, D.C. (2020)


# Section 5: Model-robust designs
## Section 5.3.2, first example: 4 points, 2 factors, first-order model

source("obj-funs.R")
source("coord-descent.R")
source("optimize-RTD.R")

# set up model
f=function(x) { c(1,x[1],x[2]) }
A=diag(c(4,4/3,4/3))

# try several tau2s 

answers <- list()
tau2s <- c(0.01,0.05,0.1, 0.25)

# loop over tau2 values
# find optimal strategies
# -- NB takes several minutes
for(k in 1:length(tau2s) ) {
  cat("\n\ntau2 = ",tau2s[k], "\n")
  answers[[k]] <- optimize.RTD(f=f, A=A, n=4, q=2, tau2=tau2s[k], sigma2.UB=1, random.starts=10, verbose=T)
}

# plot the strategies
par(mfrow=c(2,2))
for (k in 1:4) {
  xi.bar <- answers[[k]]$xi.bar
  delta <- answers[[k]]$delta
  plot(xi.bar, xlim=c(-1,1), ylim=c(-1,1), pch=19, xlab=expression(x[1]), ylab=expression(x[2]), main=substitute(tau^2 /bar(sigma)^2==val , list(val=tau2s[k])))
  for (l in 1:4) {
    polygon(matrix(xi.bar[l,], ncol=2,nrow=4, byrow=T) + delta/2 *matrix(c(-1,-1, -1,1, 1,1, 1,-1),ncol=2, byrow=T) )
  }
}

save(tau2s,answers, file="out-2factor.Rdata")


