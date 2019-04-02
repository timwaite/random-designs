#
#  Supporting material for "Minimax efficient random experimental designs, with application to 
#	model-robust design for prediction" by T.W. Waite and D.C. Woods (2019)

# Co-ordinate descent algorithm used for Section 5 - model-robust designs

coord.descent <- function(par, objfun, lower=rep(-1,length(par)), upper=rep(1,length(par)), tol=0.0001, pass.max=100, verbose=T,...) {
  
  best <- Inf
  stop <- F
  
  matrix <- F; N <- length(par)  
  
  xcurr <- par
  xnew <- xcurr
  
  ofcurr <- objfun(xcurr, ...)
  ofvals <- ofcurr
  
  pass.ctr <- 1
  
  while (!stop & pass.ctr<=pass.max) {
    if (verbose) { cat("Pass ", pass.ctr, ", current objective ", ofcurr, "\n") }
    for (i in 1:N) {
      objfuni <- function(xch) { x <- xnew; x[i] <- xch; objfun(x, ...) }
      opt <- optimize(Vectorize(objfuni), lower=lower[i], upper=upper[i])
      if (verbose) { 
        #xs <- seq(from=0.95*opt$minimum,to=1.05*opt$minimum,len=100); 
        #xs <- seq(from=lower[i],to=upper[i],len=100); 
        #plot(xs, sapply(xs,objfuni),type="l", main=i, xlim=c(lower[i],upper[i]), ylim=c(0*opt$objective, 2*opt$objective)); 
        #abline(v=xnew[i]); 
        #abline(v=opt$minimum,lty=2)   
      }
      if (opt$objective < objfun(xnew, ...)) { 
        xnew[i] <- opt$minimum
        # print(c(pass.ctr, i , xnew[i], opt$objective))
      }
    }
    
    ofnew <- objfun(xnew, ...)
    if ( ofcurr-ofnew < tol) { stop <- T }
    xcurr <- xnew 
    ofcurr <- ofnew
    ofvals <- c(ofvals, objfun(xcurr,...))
    pass.ctr <- pass.ctr + 1
  }
  return(list(xcurr=xcurr,ofvals=ofvals))
}
