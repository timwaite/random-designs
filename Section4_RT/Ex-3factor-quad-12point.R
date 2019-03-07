#
#  Supporting material for "Minimax efficient random experimental designs with application to 
#	model-robust design for prediction" by T.W. Waite and D.C. Woods (2019)

# Section 4: Model-robust designs
## Section 4.3.3, 12 points, 3 factors, quadratic model

source("obj-funs.R")
source("coord-descent.R")
source("optimize-RTD.R")

f <- function(x) { c(1,x[1],x[2],x[3], x[2]*x[3], x[1]*x[3], x[1]*x[2], x[1]^2, x[2]^2, x[3]^2)}

#
# Caclulate A matrix
#
A <- matrix(0, ncol=10, nrow=10)

indxs <- matrix( c(0,0,0,
                   1,0,0,
                   0,1,0,
                   0,0,1,
                   0,1,1,
                   1,0,1,
                   1,1,0,
                   2,0,0,
                   0,2,0,
                   0,0,2), ncol=3, nrow=10, byrow=T)

INT <- function(idx) { 2/(idx+1) * ( (idx+1) %% 2 == 1) }

for (i in 1:10) {
  for (j in 1:10) {
    A[i,j] <- prod( sapply(indxs[i,] + indxs[j,], INT) )
  }
}

library(MASS); fractions(A)

# 
# Find optimal strategy 
# - it is time consuming to perform thorough search
# - the below took around 70 min on a 2.5Ghz Macbook Pro

tau2=2^3*1*0.05^2
system.time( ans <- optimize.RTD(f=f, A=A, n=12, q=3, tau2=tau2, sigma2.UB=1, random.starts=20, verbose=T, pass.max=30, tol=1e-5) )
#save(tau2,ans,"3factor-out.Rdata")

#
# Heuristic strategies
# 
# It is quicker to obtain an efficient heuristic strategy based on a V-optimal deterministic design
#


####
#
# V-optimal deterministic design 
#
#

library(lhs)
init <- c(2*randomLHS(n=12,k=3)-1)
ansV <- coord.descent(init, Vobfun, f=f, n=12, q=3, A=A)   # best value 3.9178, worth trying several repeats
des <- data.frame(matrix(ansV$xcurr,nrow=12,ncol=3))
des <- round(des,digits=3) 
des <- des[do.call(order,des),]


#
# Set up heuristic xi.bar as described in text 
# - N.B. - first check that rows 6 & 7 of des are replicates
#

clip <- function(x, delta=0.1) {
  x <- as.matrix(x)
  x<-  pmin(  x, 1-delta/2 )
  x <- pmax( x, -1+delta/2)
  return(x)
}

# 
heur.des <- function(delta) {
  xi <- as.matrix(des)
  xi <- clip(xi, delta)  
  xi[6,] <- xi[6,] + c(-delta/2,0,0)
  xi[7,] <- xi[7, ] + c(delta/2,0,0)
  return(xi)
}


Tmc <- 2*randomLHS(n=100, k=3)-1 
Tmax <- rbind( 2*randomLHS(50,3)-1, as.matrix( expand.grid(c(-1,1), c(-1,1),c(-1,1)) ) )
# NB inclusion of the corner points in Tmax stabilizes the results as the worst t is usually at one of these points

tau2=2^3*1*0.05^2 # sd of psi at randomly selected point is 5% that of the random error sd

#
# Risk bound for heuristic strategy, as a function of delta
# - optimize delta and plot the risk bound
#

PsiH <- function(delta) { Psi.approx( xi.bar=heur.des(delta), delta, Tmc, f, A, sigma2.UB=1, Tmax, tau2=tau2 ) } 
optH <- optimize(Vectorize(PsiH), lower=0.05, upper=0.5)  
optH
plot(Vectorize(PsiH), from=0.01,to=0.5,type="l", ylim=c(0,5*optH$objective),xlim=c(0,0.5), ylab=expression(hat(Psi)(bar(xi)[delta],delta)), xlab=expression(delta))

#
# Use this as the initialization for a co-ordinate descent
#

delta <- optH$minimum
xi.bar <- heur.des(delta)
init <- c( xi.bar/(1-delta/2), delta/min(dist(xi.bar,method="maximum")) )
design.raw(init, n=12,q=3)
ans2 <- coord.descent( init, objfun=Psi.approx.wrap, lower=rep(c(-1,0), c(36,1)), upper=rep(1,37),  n=12, q=3, Tmc=Tmc, f=f, A=A, sigma2.UB=1, Tmax=Tmax, tau2=tau2  )


#
# how much bigger is the RMISPE due to presence of discrepancy?
#
sqrt( optH$objective/ min( ans$ofvals) )


#
# Bound on loss distribution survivor function, for optimized heuristic strategy
#

par(mfrow=c(1,1))
plot(NULL,xlim=c(0,25), ylim=c(0,1), xlab=expression(u/bar(sigma)^2), ylab="",lwd=2)
lines(seq(from=0,to=25,len=2), c(1,1), col="grey", lwd=2)
bound <- function(u) { min(1, optH$objective/u ) }
us <- seq(from=0,to=25,len=100)
lines(us, sapply(us,bound), lwd=2)


#
# V-efficiency distribution for optimized heuristic strategy 
#

Veffs <- NULL
for (l in 1:1000) {
  xi <-  heur.des(optH$minimum) + runif(2, max=optH$minimum/2)
  Veffs[l] <- tr(solve(M(heur.des(0),f),A))/tr(solve(M(xi,f), A))
}
summary(Veffs)

par(mfrow=c(1,1))
plot(Vectorize(PsiH), from=0.01,to=0.5,type="l", ylim=c(0,5*optH$objective),xlim=c(0,0.5), ylab=expression(hat(Psi)(bar(xi)[delta],delta)), xlab=expression(delta))
abline(v=optH$minimum,lty=2)



boxplot(Veffs, xlab="V-efficiency distribution", ylim=c(0,1))



