#
#  Supporting material for "Minimax efficient random experimental designs with application to 
#	model-robust design for prediction" by T.W. Waite and D.C. Woods (2019)

# Section 4: Model-robust designs
## Section 4.3.2, first example: 3 points, 1 factors, quadratic model

source("obj-funs.R")
source("coord-descent.R")
source("optimize-RTD.R")

# set up model
f=function(x) { c(1,x,x^2) }
A=diag(c(4,4/3,4/5)); A[1,3] <- 4/3 -> A[3,1]

# try several values of tau^2

# compute optimal strategies
# -- N.B. takes several minutes 

tau2s <- c(0.01, 0.05, 0.10, 0.25, 0.5)
answers <- list()

for (k in 1:length(tau2s)) {
  cat("tau2 = ",tau2s[k],"\n")
  answers[[k]] <- optimize.RTD(f=f, A=A, n=3, q=1, tau2=tau2s[k], sigma2.UB=1, random.starts=5, verbose=T)
}

# plot the strategies
plot(NULL, xlim=c(-1,1), ylim=c(1,5) ,ylab=expression(tau^2/bar(sigma)^2), xlab=expression(x),yaxt='n')
axis(2, at=1:5, labels=tau2s[1:5])
for( k in 1:length(tau2s)) {
  xi.bar <- answers[[k]]$xi.bar
  delta <- answers[[k]]$delta
  points(x=xi.bar, y=rep(k,3), pch=19)
  for( l in 1:3) {
    lines(x=xi.bar[l]+delta/2* c(-1,1), y=rep(k,2) )
  }
}


save(tau2s,answers,file="out-1factor.Rdata")

