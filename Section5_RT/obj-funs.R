#
#  Supporting material for "Minimax efficient random experimental designs, with application to 
#	model-robust design for prediction" by T.W. Waite and D.C. Woods (2019)

#
#  Section 5 - model-robust designs
#   - functions to compute the approximate risk bound hat(Psi) described in Section 4.3.1
#	- also includes functions to implement the design reparameterization used to convert to box constraints


Fmat <- function(xi,f) {
  t(apply(xi,1,f))
}

M <- function(xi,f) {
  Fm <- Fmat(xi,f)
  t(Fm)%*%Fm 
}
tr <- function(X) { sum(diag(X)) }

Vobfun <- function(x,f, n,q, A) {
  xi <- matrix(x, nrow=n, ncol=q) 
  M1 <- M(xi,f)
  tr(solve(M1, A))
}

d <- function(xi.bar, t) { 
  xi.bar + matrix(t, nrow=nrow(xi.bar), ncol=ncol(xi.bar), byrow=T) 
}

MIV.approx <- function( xi.bar, delta, Tmc, f, A, sigma2.UB) {		
  tmp <- 0 
  for (i in 1:nrow(Tmc)) {
    xi <- d(xi.bar, Tmc[i,] *delta/2 )
    tmp <- tmp + tr(solve(M(xi,f),A)) /nrow(Tmc)
  }
  sigma2.UB * tmp
}

K <- function(xi,f,A) {
  F1 <- Fmat(xi,f)
  M1 <- t(F1)%*%F1
  MiFt <- solve( M1, t(F1) )
  t(MiFt)%*%A%*%MiFt
}

biassq.approx <- function( t, xi.bar, delta, tau2, f, A) {
  q <- ncol(xi.bar)
  xi <- d(xi.bar, t*delta/2)
  K1 <- K(xi, f, A)
  ans <- tau2 + tau2 *max(eigen(K1,symmetric=T)$values) / (delta^q) 
  ans
}

Psi.approx <- function( xi.bar, delta,  Tmc, f, A, sigma2.UB, Tmax, tau2) {
  
  MIV1 <- MIV.approx(xi.bar, delta, Tmc, f, A, sigma2.UB)
  
  biassq1 <- rep(0,nrow(Tmax))
  for (i in 1:nrow(Tmax)) {
    biassq1[i] <- biassq.approx( Tmax[i,], xi.bar, delta, tau2, f, A)
  }
  
  MIV1 + max(biassq1)
}

design.raw <- function(x, n, q, transform.delta=T ) {
  N <- length(x)
  if (N!= n*q+1) stop("Incorrect dimensions")
  ctilde <- matrix(x[-N], nrow=n, ncol=q, byrow=F)
  dtilde <- x[N]
  dmin <- min( dist(ctilde,method="maximum") )
  delta <- 2* dtilde*dmin/( 2 + dmin* dtilde)
  xi.bar <- ctilde * (1-delta/2)
  list(xi.bar = xi.bar, delta=delta)
}

Psi.approx.wrap <- function(x, n, q, Tmc, f, A, sigma2.UB, Tmax, tau2, transform.delta=T) {
  des1 <- design.raw(x,n,q, transform.delta)
  Psi.approx(des1$xi.bar, des1$delta, Tmc, f, A, sigma2.UB, Tmax, tau2)
}

# design.rawOLD <- function(x, n, q, transform.delta=T ) {
#   N <- length(x)
#   if (N!= n*q+1) stop("Incorrect dimensions")
#   xi.bar <- matrix(x[-N], nrow=n, ncol=q, byrow=F)
#   delta.max <- min( dist(xi.bar,method="maximum"), 2*abs(xi.bar-1), 2*abs(xi.bar+1)  )
#   if (transform.delta) {	delta <- x[N] * delta.max } else {delta <- x[N]}
#   list(xi.bar = xi.bar, delta=delta, delta.max=delta.max)
# }
# 
# Psi.approx.wrapOLD <- function(x, n, q, Tmc, f, A, sigma2.UB, Tmax, tau2, transform.delta=T) {
#   des1 <- design.raw(x,n,q, transform.delta)
#   Psi.approx(des1$xi.bar, des1$delta, Tmc, f, A, sigma2.UB, Tmax, tau2)
# }
 