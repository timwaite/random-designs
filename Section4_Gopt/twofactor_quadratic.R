#
#  Supporting material for "Minimax efficient random experimental designs with application to 
#	model-robust design for prediction" 
# by Waite, T.W. and Woods, D.C. (2020)

#  Section 4: G-optimal random design examples
#  Two factor quadratic model - Section 4.1.4

library(gmp)  # used to perform exact inversion of (fractional) information matrices -- reduces numerical error

des_sp <- as.matrix(expand.grid(c(-1,0,1),c(-1,0,1)))

designs <- t(combn(1:9, 6))
designs <- designs[-c( 1, 20, 24, 29, 43, 63, 69, 84),]  #exclude designs that later turn out to be singular

# information matrix for a design defined by indices
M <- function(design,des_sp,frac=T) {
	Xmat <- des_sp[as.vector(design),]
	Fmat <- as.matrix( cbind(1, Xmat, Xmat[,1]*Xmat[,2], Xmat[,1]^2, Xmat[,2]^2) )
	if (frac) Fmat <- as.bigq(Fmat)
	colnames(Fmat) <- NULL
	return(t(Fmat)%*%Fmat)
}

M(designs[1,],des_sp)

EL <- function(design,x,des_sp) {
	fx<-c(1,x,x[1]*x[2],x[1]^2, x[2]^2)
	tmp <- try(solve(M(design,des_sp)),silent=TRUE)
	if(class(tmp)!="try-error") return(	t(fx)%*%tmp%*%fx )
	else return(Inf)
}

EL(designs[76,], c(0,0), des_sp)
EL(designs[8,], c(-1,-1), des_sp)


Lmat <- matrix(nrow=nrow(designs),ncol=nrow(des_sp))

for(i in 1:nrow(designs)) {
	for(j in 1:nrow(des_sp)) {
		Lmat[i,j] <- as.numeric(EL(designs[i,],des_sp[j,],des_sp))
	}
}


min(apply(Lmat,1,max))
bestdes <- which.min(apply(Lmat,1,max))
# optimal deterministic design - max exptected loss = 2.75

library(lpSolveAPI)

A <- t(Lmat)
my.lp2 <- make.lp( nrow=nrow(A), ncol=ncol(A) )
lp.control(my.lp2,sense="max")   # make it a maximization problem
for(i in 1:ncol(A)) set.column(my.lp2,i,A[,i])  
set.objfn(my.lp2, rep(1,ncol(A)))
set.constr.type(my.lp2,rep("<=",nrow(A)))
set.rhs(my.lp2, rep(1, nrow(A)))
solve(my.lp2)
get.objective(my.lp2)
optq <- get.variables(my.lp2)
optq <- optq/sum(optq)


max(optq%*%Lmat)

#
# plot the random design strategy
# - this gives Fig. 2 (left panel) from the paper
#

which(apply(Lmat[optq>0,],1,max) == 2.75)

sptdes <- which(optq>0)
par(mfrow=c(3,3))
for(i in 1:length(sptdes)) {
	
	
	main <- substitute(pi(xi[i])==main,list(i=i,main=round(optq[sptdes[i]],4)))
	#main <- paste(round(optq[sptdes[i]],4))
	
	
	plot( des_sp[designs[sptdes[i],],] ,xlim=c(-1,1),ylim=c(-1,1),
		main=main,pch=19, xlab=expression(x[1]),ylab=expression(x[2]))
}
plot( des_sp[designs[bestdes,],] ,xlim=c(-1,1),ylim=c(-1,1), pch=19, main="best deterministic",xlab=expression(x[1]),ylab=expression(x[2]))


Lmat[optq>0,]



# model: full quadratic in 2 factors: p=6, n=6 i.e. smallest possible number of runs
# 
# optimal deterministic design - max expected loss = 2.75
# optimal randomized design contains 8 support designs 
# max expected loss is 1.55


##### 
#
# What about the ROAD? (Randomly-rounded optimal approximate design)
#
# First need to find the D/G-optimal approximate design.
#

X <- as.matrix(expand.grid(c(-1,0,1), c(-1,0,1)))
f <- function(x) { c(1,x[1],x[2],x[1]*x[2],x[1]^2,x[2]^2) }

M2 <- function(w) {
  Fmat <- t(apply(X,1,f) ) 
  out <- t(Fmat)%*%diag(w)%*%Fmat
  return(out)
}
Psi <- function(w) { 
  lambda <- eigen(M2(w),symmetric=T)$values
  sum(log(lambda))
}
Phi <- function(x,w=rep(1/9,9)) {
  fx <- f(x)
  t(fx)%*%solve(M2(w), fx) 
}

Shift <- function(wk, k=1, Des=rep(1/9,9)) {
  DesNew <- Des
  DesNew[k] <- wk
  DesNew[-k] <- (1-wk)*DesNew[-k]/sum(DesNew[-k])
  DesNew
}

Psi2 <- function(wk, k=1, Des=rep(1/9,9)) {
  -Psi(w = Shift(wk, k, Des))
}

f(c(0.5,1))
f(X[4,])
psi(rep(1/9,9))
apply(X, 1, Phi, w=rep(1/9,9))

Des <- rep(1/9,9)
for(i in 1:100) {
  for (k in 1:9) {
    opt <- optimize(Psi2, interval=c(0,1), k=k, Des=Des)
    wk <- opt$minimum
    Des <- Shift(wk, k, Des)
    print(opt$objective)
  }
}
Des

round(Des,digits=4)

#
# Kiefer roundings of the G(=D)-optimal approximate design 
#  - also 'randomized roundings'
#
#####
# indices of the design points in des_sp
# corner points 1,3,7,9
# center 5
# edge mid-points 2,4,6,8

# the Kiefer roundings contain all 4 corner points and the center point, plus one edge mid-point
# there are 4 different rounded designs corresponding to the choice of mid-point

# 1,2,3,5,7,9 <-- design 14
# 1,3,4,5,7,9 <-- design 36
# 1,3,5,6,7,9 <-- design 42
# 1,3,5,7,8,9 <-- design 44

kieferq <- rep(0,nrow(designs))
kieferq[c(14,36,42,44)] <- 1/4


# Efficiencies
# fixed Kiefer rounding of approximate design
delta <- function(i) { out <- rep(0,length=nrow(designs)); out[i] <-1; return(out) } 
max(optq%*%Lmat)/max(delta(14)%*%Lmat)
# randomized Kiefer rounding of approximate design
max(optq%*%Lmat)/max(kieferq%*%Lmat)



#####
#
#  Bounds on loss distribution survivor function
#
#####

loss.surv <- function(l,design.prob=optq,x.ind=1,s2=1) {
	ans <- 0
	for (i in 1:nrow(designs)) {
			x <- des_sp[x.ind,]
			f <- c(1,x,x[1]*x[2], x[1]^2, x[2]^2)
			#qform <- t(f)%*% solve(M(designs[i,],des_sp,frac=F)) %*% f   # 1x1 BigRational matrix
			predsd <- s2* (1+ t(f)%*% solve(M(designs[i,],des_sp,frac=F)) %*% f)
				# do not use exact fractional computation, slows things down too much here
			term <- pnorm(sqrt(l), mean=0, sd=predsd ) - pnorm(-sqrt(l),mean=0,sd=predsd)
			ans <- ans + design.prob[i] * (1-term)
	}
	return(ans)
}
loss.surv(1,optq,1,1)

wc.loss.surv <- function(l,design.prob=optq,s2=1) {
	vals<- rep(NA,nrow(des_sp))
	for (i in 1:nrow(des_sp)) {
		vals[i] <- loss.surv(l,design.prob,i,s2)
	}
	max(vals)
 }
wc.loss.surv(1.3) 
 
# plotting grid
lgrid <- exp(seq(from=log(0.01),to=log(100),length=100))

#loss distribution for optimal random strategy 
formals(wc.loss.surv)$design.prob = optq
opt.rds.loss.surv <- sapply(lgrid, Vectorize(wc.loss.surv))

delta8 <- delta(8)

# loss distribution for minimax deterministic
formals(wc.loss.surv)$design.prob= delta(8)
mm.det.loss.surv <- sapply(lgrid, Vectorize(wc.loss.surv))

# loss distribution for fixed rounded approximate design
formals(wc.loss.surv)$design.prob = delta(14)
fixed.kiefer.loss.surv <- sapply(lgrid, Vectorize(wc.loss.surv))

# loss distribution for randomly-rounded approximate design
formals(wc.loss.surv)$design.prob = kieferq
random.kiefer.loss.surv <- sapply(lgrid, Vectorize(wc.loss.surv))

#
# plot the bounds for the survivor function or optimal RDS and minimax deterministic
# -  this is Fig.2 (Right panel) from the paper
#
par(mfrow=c(1,1))
plot(sqrt(lgrid), opt.rds.loss.surv, type="l", main="", ylab="", xlab=expression(u/bar(sigma)^2) ,lwd=4) 
lines(sqrt(lgrid), mm.det.loss.surv, type="l",col="grey",lwd=4)
lines(sqrt(lgrid), fixed.kiefer.loss.surv, type="l", col=1,lty=2, lwd=2)
lines(sqrt(lgrid), random.kiefer.loss.surv, type="l",col="lightgrey",lwd=4)
