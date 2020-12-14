#
#  Supporting material for "Minimax efficient random experimental designs with application to 
#	model-robust design for prediction" 
# by Waite, T.W. and Woods, D.C. (2020)

#  Section 4: G-optimal random design examples
#  One factor quadratic model - Section 4.1.3


designs <- expand.grid(c(-1,0,1),c(-1,0,1),c(-1,0,1),c(-1,0,1))		
# n=4 - to obtain results for n=4, uncomment this line and comment out the line below 
#designs <- expand.grid(c(-1,0,1),c(-1,0,1),c(-1,0,1),c(-1,0,1),c(-1,0,1))	# n=5 - to obtain results for n=5, uncomment this line 

# restrict to designs with at least 3 distinct support points, since the quadratic model has 3 parameters
designs <- designs[sapply(apply(designs,1,unique),length)>=3,]	 

# information matrix
M <- function(xi) {
	X <- cbind(1,xi,xi^2)
	t(X)%*%X
}
M(t(designs[2,]))

# conditional variance of prediction at x, given the design
EL <- function(xi,x) {
	tmp <- try(solve(M(xi),c(1,x,x^2)),silent=TRUE)
	if(class(tmp)!="try-error") return(	t(c(1,x,x^2))%*%tmp)
	else return(Inf)
}

EL(t(designs[3,]),1)


# calculate (transposed) modified risk matrix
Lmat <- matrix(nrow=nrow(designs),ncol=3)
for(i in 1:nrow(designs)){
	for(j in 1:3) {
		Lmat[i,j] <- EL(t(designs[i,]),c(-1,0,1)[j])	
	}	
}

# inspect the designs and the risk profiles
cbind(designs,Lmat, apply(Lmat,1,max))
min(apply(Lmat,1,max))
designs[which.min(apply(Lmat,1,max)),]


L <- t(Lmat)
library(lpSolveAPI)

# set up linear programming problem -- equation (3) from Section 2.1
my.lp2 <- make.lp(nrow=nrow(L), ncol=ncol(L))

# constraints
for(i in 1:ncol(L)) set.column(my.lp2,i,L[,i])
set.constr.type(my.lp2, rep("<=", nrow(L)))
set.rhs(my.lp2, rep(1,nrow(L)))

# set up objective function
lp.control(my.lp2,sense="max")   # make it a maximization problem
set.objfn(my.lp2, rep(1,ncol(L)))

# solve
solve(my.lp2)
pi_dagger <- get.variables(my.lp2)
pi_star <- pi_dagger/sum(pi_dagger)

designs[ pi_star>0, ]

max( L%*%pi_star )

# n=4
# the optimal randomized design assigns weight (1/3) to (1,0,-1,-1), (1,0,0,-1) and (1,1,0,-1)
# the effective number of replications at -1,0,1 is 4/3
# equivalently, always use (-1,0,1) plus an extra point chosen uniformly at random from {-1,0,1}
# the max expected loss for the deterministic design is 1, for the randomizaed design it is 0.8333

# n=5
# the optimal randomized design contains two replicates of two points from (-1,0,1) and one replicate of the remaining point. The single-rep point is chosen uniformly at random
# the max expected loss for the deterministic design is 1, for the randomized design it is 2/3



