#
#  Supporting material for "Minimax efficient random experimental designs with application to 
#	model-robust design for prediction" by T.W. Waite and D.C. Woods (2019)

# Section 4: Model-robust designs
## Section 4.3.2, first example: 4 points, 2 factors, first-order model


source("obj-funs.R")
source("coord-descent.R")

# set up model
f=function(x) { c(1,x[1],x[2]) }
A=diag(c(4,4/3,4/3))

# set up discretized translation set (for objective function evaluation)
library(lhs)
Tmc <- 2*randomLHS(n=100, k=2)-1


# try several tau2s 

answers <- list()
rawdesigns <- list()

tau2s <- c(0.01,0.05,0.1, 0.25)
xi <- 0.8*(randomLHS(n=4, k=2)-0.5)
xi <- matrix(runif(8, min=-0.8, max=0.8),ncol=2,nrow=4)

# find optimal strategies

# loop over tau2 values
for(k in 1:length(tau2s) ) {
  tmp <- NULL
  answers[[k]] <- NULL
  
  # multiple random starts - keep best so far
  for (m in 1:20) {
    xi <- 0.8*(randomLHS(n=4, k=2)-0.5)
    tmp <- coord.descent(c(xi,1 ), Psi.approx.wrap, lower=rep(c(-1, 0 ), c(8, 1)), upper=rep(c(1,1),c(8,1)),  
                               n=4, q=2, Tmc=Tmc, f=f, A=A, sigma2.UB=1, Tmax=Tmc, tau2=tau2s[k], tol=1e-6 )
    if (length(answers)<k) { 
      answers[[k]] <- tmp 
    } else if (min(tmp$ofvals)<min(answers[[k]]$ofvals)) { 
      answers[[k]] <- tmp
    }
  }
  rawdesigns[[k]] <- design.raw(answers[[k]]$xcurr, n=4, q=2)
}

# plot the strategies
par(mfrow=c(2,2))
for (k in 1:4) {
  xi.bar <- rawdesigns[[k]]$xi.bar
  delta <- rawdesigns[[k]]$delta
  plot(xi.bar, xlim=c(-1,1), ylim=c(-1,1), pch=19, xlab=expression(x[1]), ylab=expression(x[2]), main=substitute(tau^2 /bar(sigma)^2==val , list(val=tau2s[k])))
  for (l in 1:4) {
    polygon(matrix(xi.bar[l,], ncol=2,nrow=4, byrow=T) + delta/2 *matrix(c(-1,-1, -1,1, 1,1, 1,-1),ncol=2, byrow=T) )
  }
}

# save(tau2s,answers,rawdesigns, file="out.Rdata")

