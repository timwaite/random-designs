#
#  Supporting material for "Minimax efficient random experimental designs, with application to 
#	model-robust design for prediction" by T.W. Waite and D.C. Woods (2019)

# Section 4: Main algorithm 
#
# Algorithm to find minimax hypercuboidal random translation design strategy for an approximate linear model
# Used in Section 4 of the paper
#
# The algorithm uses multiple random initializations of a co-ordinate descent algorithm to find an approximate minimax strategy. 
# The objective function is evaluated as described in Section 4.3.1, using Monte Carlo integration for the first term 
#  of equation (14) in the paper. A discretized translation set to perform the maximization for the third term. 
# A transformation is used to convert the problem into a box-constrained optimization. 
# The discretized translation set is iteratatively refined using an approach similar to Pronzato and Walter (1987)
#
# Inputs:
# - f is the vector of regressor functions - given a vector x in the design space, f returns (f_0(x), ... f_p(x))
# - A is the integral of f(x)f^T(x) over [-1,1]^q (needs to be specified by the user)
# - n is the number of design points
# - q is the number of factors (each assumed to take values in [-1,1])
# - tau2 is the maximum L2-norm of the discrepancy function 
# - sigma2.UB is the upper bound on the random error variance
# - Tmc is the Latin Hypercube Sample of possible translations (for the Monte Carlo integration) - defaults should be ok
# - Tmax.init is the initial discretization of the translation (for maximization part) - defaults should be ok 
# - tol is a tolerance for the optimization
# - random.starts is the number of random initializations used in the co-ordinate descent algorithm
# - pass.max is the maximum number of passes using cyclic co-ordinate descent
# - verbose = T gives further progress messages and plots


library(lhs)

optimize.RTD <- function(f, A, n, q, tau2, sigma2.UB=1, Tmc=(2*randomLHS(n=100,k=q)-1), Tmax.init=as.matrix(expand.grid(rep(list(c(-1,1)),q))), tol=1e-6, random.starts=20, pass.max=100, verbose=F) {
  
  Tmax <- Tmax.init
  
  stop <- F
  while (!stop) {
    
    opt.des <- NULL
    
    initial <- list()
    for( m in 1:random.starts) {
        initial[[m]] <- c(2*randomLHS(n,q)-1,1)
    }
 
    for (m in 1:length(initial)) {
      cat("Optimizing design strategy using co-ordinate descent, random initialization ", m, "/", random.starts, "\n")
      init.des <- initial[[m]]
      # ideally do multiple random starts here and save the best 
      tmp <- coord.descent(init.des, Psi.approx.wrap, lower=rep(c(-1, 0 ), c(n*q, 1)), upper=rep(c(1,1),c(n*q,1)),  tol=tol, pass.max=pass.max, verbose=verbose,
                               n=n, q=q, Tmc=Tmc, f=f, A=A, sigma2.UB=sigma2.UB, Tmax=Tmax, tau2=tau2)
 
      if (is.null(opt.des)) { 
        opt.des <- tmp 
        } else if (min(tmp$ofvals)< min(opt.des$ofvals)) {
        opt.des <- tmp
        }
    }
    # check if the discretization of the translation set needs refining
    cat("Checking if refinement needed in discretization of translation set...\n")
    desraw <- design.raw(opt.des$xcurr, n=n, q=q) 
    worst.t <- coord.descent(rep(0,q), function(x) { -biassq.approx(x, xi.bar=desraw$xi.bar, delta=desraw$delta, tau2=tau2, f=f, A=A)} )
    ofnew1 <- MIV.approx(xi.bar=desraw$xi.bar, delta=desraw$delta, Tmc=Tmc, f=f, A=A, sigma2.UB=sigma2.UB)
    ofnew2 <- biassq.approx( t=worst.t$xcurr, xi.bar=desraw$xi.bar, delta=desraw$delta, tau2=tau2, f=f, A=A)
    if ( ofnew1+ofnew2 > min(opt.des$ofvals) + tol) { 
        Tmax <- rbind(Tmax, worst.t$xcurr)
        cat(" -- Refining discretization.\n")
    } else {
        stop <- T
        cat(" -- No further refinement needed.\n")
      }
   }
   return(list(xi.bar=desraw$xi.bar, delta=desraw$delta, ofval=min(opt.des$ofvals)))
}


