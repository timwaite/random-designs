#
#  Supporting material for "Minimax efficient random experimental designs with application to 
#	model-robust design for prediction" by T.W. Waite and D.C. Woods (2019)
#
# 
#  Section 3: G-optimal random design examples
#  One factor first-order model - Section 3.1.1

# completely enumerate all possible three point designs over X= {-1,0,1}

designs <- expand.grid(c(-1,0,1),c(-1,0,1),c(-1,0,1))
designs <- designs[-c(1,14,27),]   # exclude the singular designs (i.e. only one support point)

# function to calculate the information matrix for the first-order model
M <- function(xi) {
	X <- cbind(1,xi)
	t(X)%*%X
}
M(t(designs[1,]))

# conditional variance of prediction at x, given the design
EL <- function(xi,x) {
	tmp <- try(solve(M(xi),c(1,x)),silent=TRUE)
	if(class(tmp)!="try-error") return(	t(c(1,x))%*%tmp)
	else return(Inf)
}

EL(t(designs[1,]),1)


# calculate the (transpose of the) modified risk matrix \tilde{R}

Lmat <- matrix(nrow=nrow(designs),ncol=3)
for(i in 1:nrow(designs)){
	for(j in 1:3) {
		Lmat[i,j] <- EL(t(designs[i,]),c(-1,0,1)[j])	
	}	
}

# present the designs, risk, and max-risk
cbind(designs,Lmat, apply(Lmat,1,max))


L <- t(Lmat)		# set up the L matrix / modified risk matrix in the correct transposition 

# set up and solve the linear programming problem from (3) in Section 2.1

library(lpSolveAPI)

# set up a linear programming problem with the correct numbers of decision variables and inequality constraints
my.lp2 <- make.lp( nrow=nrow(L), ncol=ncol(L) )  

# set up the linear inequality constraints:
for (i in 1:ncol(L)) set.column(my.lp2,i,L[,i])	# set coefficients of linear constraints
set.constr.type(my.lp2, rep("<=",nrow(L))) 		# set direction of constraints
set.rhs(my.lp2, rep(1,nrow(L)))					# set the bounds
# NB lpSolve also include non-negativity constraints by default.

# set up the linear objective function.
lp.control(my.lp2,sense="max")   # make it a maximization problem
set.objfn(my.lp2, rep(1,ncol(L)))   

# solve the linear program and extract the optimal values of the decision variables
solve(my.lp2)
get.objective(my.lp2)
( pi_dagger <- get.variables(my.lp2) )

pi_star <- pi_dagger/sum(pi_dagger)	
# normalize the decision vector to find the optimal mixed strategy

max(pi_star%*%Lmat)		# max risk of optimal strategy


designs[ pi_star>0,]



