#
#  Supporting material for "Minimax efficient random experimental designs with application to 
#	model-robust design for prediction" by T.W. Waite and D.C. Woods (2019)

# Section 4: Model-robust designs
## Section 4.3.2, first example: 3 points, 1 factors, quadratic model

source("obj-funs.R")
source("coord-descent.R")



# set up model
f=function(x) { c(1,x,x^2) }
A=diag(c(4,4/3,4/5)); A[1,3] <- 4/3 -> A[3,1]

# set up discretized translation set (for objective function evaluation)
library(lhs)
Tmc <- 2*randomLHS(n=50, k=1)-1

# set up initial mean-design
xi <- 0.8*(randomLHS(n=3, k=1)-0.5)

tau2=2*1*0.01^2

ans <- coord.descent(c(xi,1 ), Psi.approx.wrap, lower=rep(c(-1, 0 ), c(3, 1)), upper=rep(c(1,1),c(3,1)),  
                     n=3, q=1, Tmc=Tmc, f=f, A=A, sigma2.UB=1, Tmax=Tmc, tau2=tau2, tol=1e-4 )

plot(ans$xcurr[1:3], rep(0,3))


# try several values of tau^2

# compute optimal strategies
tau2s <- 2*1* c(0.01, 0.05, 0.10, 0.25, 0.5, 1)^2
answers <- list()
rawdesigns <- list()
for (k in 1:length(tau2s)) {
  answers[[k]] <- coord.descent(c(xi,1 ), Psi.approx.wrap, lower=rep(c(-1, 0 ), c(3, 1)), upper=rep(c(1,1),c(3,1)),  
                                n=3, q=1, Tmc=Tmc, f=f, A=A, sigma2.UB=1, Tmax=Tmc, tau2=tau2s[k], tol=1e-4 )
  rawdesigns[[k]] <- design.raw(answers[[k]]$xcurr, n=3, q=1, transform.delta=T)
}


# plot the strategies
plot(NULL, xlim=c(-1,1), ylim=c(1,5) ,ylab=expression(tau^2/bar(sigma)^2), xlab=expression(x),yaxt='n')
axis(2, at=1:5, labels=tau2s[1:5])
for( k in 1:length(tau2s)) {
  xi.bar <- rawdesigns[[k]]$xi.bar
  delta <- rawdesigns[[k]]$delta
  points(x=xi.bar, y=rep(k,3), pch=19)
  for( l in 1:3) {
    lines(x=xi.bar[l]+delta/2* c(-1,1), y=rep(k,2) )
  }
}



