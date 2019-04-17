#
#  Supporting material for "Minimax efficient random experimental designs, with application to 
#	model-robust design for prediction" by T.W. Waite and D.C. Woods (2019)
#
#	Computations for Section 3 - Potential outcomes designs 
#   - three factors, 20 runs, full quadratic model, factors take levels -1,0,1
#   - use (discrete) co-ordinate exchange to optimize designs
#   - in code, define designs as a matrix

des <- matrix(sample(c(-1,0,1), 20*3, replace=T), nrow=20, ncol=3)

# define quadratic regression model in three factors
f <- function(x) { c(1, x[1], x[2], x[3], 
                     x[1]*x[2], x[2]*x[3], x[1]*x[3], 
                     x[1]^2, x[2]^2, x[3]^2  ) }

# A_S objective function

Psi.crrod <- function(design) {
	Xmat <- design
	Fmat <- t(apply(Xmat,1,f))
	p <- ncol(Fmat) -1 
	n <- nrow(Xmat)
	colnames(Fmat) <- NULL
	M <- t(Fmat)%*%Fmat
	Minv <- try(solve(M), silent=T)
	if (class(Minv)!="try-error") {
		S <- t(diag(1, p +1 , p+1)[,-1])
		return( (n/(n-1)) * sum(diag(Minv%*%t(S)%*%S) ))
		} else return(Inf)	
}

#debugonce(AS_objfun)
Psi.crrod(des)
 
# maximum risk for a deterministic design
# computed using class Theta2; otherwise this is a lower bound for the maximum 

Psi.det <- function(design) {
	Xmat <- design
	Fmat <- t(apply(Xmat,1,f))
	p <- ncol(Fmat) -1 
	n <- nrow(Xmat)
	colnames(Fmat) <- NULL
	M <- t(Fmat)%*%Fmat
	Minv <- try(solve(M), silent=T)
	
		if (class(Minv)!="try-error") {
		S <- t(diag(1, p +1 , p+1)[,-1])
		Amat <- (t(S)%*%S %*% Minv)
		return( n * max(eigen(Amat)$values)  )		
		} else return(Inf)	
	
}

#debugonce(Psi.det)
Psi.det(des)
 
# discrete co-ordinate exchange algorithm 

coord.exchange <- function(init_des, objfun = AS_objfun, tol=1e-6, verbose=F) {
  stop <- F
  cur_des <- init_des
  of_best <- objfun(init_des) # need to catch failed evaluations
  itno <- 0
  
  # keep making passes until I say stop
  while (!stop) {
    stop <- T   # stop after this pass unless it makes an improvement
    itno <- itno +1
    for (i in 1:nrow(cur_des)) {  # iterate over rows
      for (j in 1:ncol(cur_des)) {  # iterate over columns
        for (x in c(-1,0,1)) { # evaluate all possible exchanges
          tmp_des <- cur_des;
          tmp_des[i,j] <- x;
          ofval <- objfun(tmp_des); # need to catch failed evaluations
          if (ofval > -Inf) { # do not save changes that lead to problems with infinities
            if ( ofval < of_best - tol ) {   
              # if the proposed exchange is beneficial, then keep it and do at least one more pass
              cur_des <- tmp_des
              of_best <- ofval;
              stop <- FALSE;  
            }
          }
         }
        # #commented code: continuous version
        # of1d <- Vectorize( function(x) { tmp_des <- cur_des; tmp_des[i,j] <- x; objfun(tmp_des) } )
        # a <- optimize(of1d, interval=c(-1,1))
        # if ( a$objective / of_best   <  1 - tol) { 
          # cur_des[i,j] <- a$minimum
          # of_best <- a$objective;
          # stop <- FALSE;
        # }
      }
    }
    if (verbose) cat("Pass " , itno,  " Obj fun ", of_best, "\n")
  }
  
  return(cur_des)
}

	
coord.exchange.multistart <- function( random.starts = 1000, objfun = AS_objfun) {
	loc.opt <- list()
	of_best <- Inf	
	of_vals <- rep(NA,random.starts)
	for (k in 1:random.starts) {
		# revise: take more care to generate nonsingular starting design
		init  <- matrix(sample(c(-1,0,1), 20*3, replace=T), nrow=20, ncol=3)
		this.start.des <- coord.exchange(init, objfun=objfun)
		ofval <- objfun(this.start.des)
		of_vals[k] <- ofval
		if (ofval < of_best) {	
			 best.des <- this.start.des
			 of_best <- ofval
		}
	}
	return(list(best.design=best.des,objfun.vals=of_vals))
}

crrod <- coord.exchange.multistart(random.starts=1000, objfun=Psi.crrod)
plot(crrod$objfun.vals)


tr <- function(A) sum(diag(A))

# Max risk for minimax random design strategy

( Psi.crrod(crrod$best.design) )


# The max risk efficiency of the unrandomized  AS-optimal design
#  compared to the randomized AS-optimal design is  1.7%
Psi.crrod(crrod$best.design) / Psi.det(crrod$best.design)

# It is not enough to randomize the run order of a poor choice of treatments.
# To demonstrate this, the code below finds a CRROD which has an MR-efficiency of 
# less than 1% relative to the best deterministic design.
#
stop <- F
tmp <- Psi.crrod(crrod$best.design)
while (!stop) {
  bad.des <- matrix( sample(c(-1,0,1), 3*20, replace=T), ncol=3, nrow=20)
  ofval <- Psi.crrod( bad.des )
   if ( ofval > 100*tmp  & ofval < Inf) {
    stop <- TRUE
  }
}
# efficiency of CRROD relative to best deterministic is 24.9%
Psi.crrod(crrod$best.des) / Psi.crrod( bad.des ); 
 
# order the design points lexicographically for pretty printing
library(kdtools)
 

library( xtable )
xtable( cbind( lex_sort(crrod$best.design), NA, 
		lex_sort(bad.des)), digits=0) 
#save(deterministic,crrod, bad.des , file="CRROD_new.Rdata")

# survivor function upper bounds (Figure 1 in the paper, left panel)

# unrandomized AS-optimal

lmax.ASunran <- Psi.det(crrod$best.design)
survbound.ASunran <- function(u) { return( u <= lmax.ASunran ) }
plot( survbound.ASunran, from=0, to= 10 ,ylab="", ylim=c(0,1) ,xlab=expression(u/sigma^2), n=1000, lwd=4, col="grey")

# randomized AS-optimal 

tmp2 <- Psi.crrod(crrod$best.design)
n <- nrow(crrod$best.design)

survbound.mMcrrod <- Vectorize( function(u) { min(1,  tmp2/u * ( u <= lmax.ASunran ) ) } )
curve( survbound.mMcrrod , from=1e-6, to=10, n=1000, add=TRUE ,lwd=4, lty=2)








