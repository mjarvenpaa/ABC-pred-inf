ifelsem <- function(test, yes, no) {
  # As ifelse but test must be scalar logical while the yes and no can be pretty
  # much anything, e.g. vectors that are returned as such based on the test. 
  # This function is meant for a one-line shorthand for the code line
  # if(test){return(yes)}else{return(no)}
  
  if (length(test)!=1) {stop('Incorrect input.')}
  if (test) {
    return(yes)
  }
  return(no)
}

rmunif <- function(n=1, a=0, b=1) {
  # Generates n random vectors from a distribution whose components are 
  # independent uniform distributed with bounds in vectors a and b.
  
  if (!is.vector(a)||!is.vector(b)||length(a)!=length(b)||any(a>b)||!all(is.finite(c(a,b)))) {
    stop('Incorrect input.')
  }
  return(drop(a + (b-a)*matrix(runif(n = n*length(a)),length(a))))
}

logpmunif <- function(x, a=-Inf, b=Inf) {
  # Returns the log-density value (up to normalization!) at x of a density whose 
  # components are independent uniform distributed with bounds in vectors a and 
  # b. That is, returns 0 if x is inside the bounds and -Inf otherwise. 
  # No input checking is made for speed.
  
  if (all(a<=x&x<=b)) { return(0) }
  return(-Inf)
}

gaussian.mixture.quantile <- function(q, w, ms, sigmas) {
  # Finds (approximate) q-quantile of Gaussian mixture with weights 'w', 
  # means 'ms' and standard deviations 'sigmas'. 
  
  gm.cdf.eval <- function(x) sum(w*pnorm(x,ms,sigmas))
  f <- function(x) gm.cdf.eval(x)-q
  ival <- c(min(ms-5*sigmas),max(ms+5*sigmas))
  return(uniroot(f,ival)$root)
}

density.bound <- function(x, lb=-Inf, ub=Inf, ...) {
  # Kernel density estimation that takes into account strict bounds of samples 
  # such as positivity (in a somewhat ad-hoc fashion). Intended only for 
  # visualization purposes. Partially based on the code given in
  # stats.stackexchange.com/questions/65866/good-methods-for-density-plots-of-non-negative-variables-in-r
  
  if (lb>=ub) {
    stop('Incorrect bounds.')
  } else if (is.finite(lb)&&is.finite(ub)) {
    stop('Currently only one bound condition is supported.')
  } else if (lb==-Inf&&ub==Inf) {
    #stop('No bounds.')
    return(density(x, ...))
  } else if (any(x<lb)||any(x>ub)) {
    stop('Samples do not obey the bounds.')
  }
  # transform data do that it is positive and truncated at 0
  if (is.finite(ub)) {
    x <- ub - x
  } else {
    x <- x - lb
  }
  h <- density(x, kernel="gaussian", ...)$bw
  w <- 1/pnorm(0, mean=x, sd=h, lower.tail=FALSE)
  suppressWarnings({d <- density(x, bw=h, kernel="gaussian", weights=w/length(x))})
  if (min(d$x)<0) {
    d$y[d$x<0] <- 0
    m <- max(d$x[d$x<0])
    d$y <- d$y[d$x>=m] # remove all but one negative grid points
    d$x <- d$x[d$x>=m]
  }
  # transform data back to original domain
  if (is.finite(ub)) {
    d$y <- rev(d$y)
    d$x <- -rev(d$x) + ub
  } else {
    d$x <- d$x + lb
  }
  return(d)
}

pred.post.stats <- function(samples, q, sm=F, x=NA) {
  # Computes mean, median and q-quantiles of some predictive distribution based
  # on given MCMC or ABC samples. 
  
  qnames <- c('mean','med','u1','l1','u2','l2')
  qs <- c(NA, 0.5, 1-q[1]/2, q[1]/2)
  if (length(q) >= 2) {
    qs <- c(qs, 1-q[2]/2, q[2]/2) # extra quantiles
  }
  pred <- list()
  for (i in 1:length(qs)) {
    if (i==1) {
      vals <- rowMeans(samples)
    } else {
      vals <- rowQuantiles(samples, probs = qs[i]) # 'matrixStats'
    }
    # smooth the values (could also use quantile regression)
    if (sm) {
      #vals <- predict(loess(vals~x)) # would require adjusting the span using e.g. CV
      vals <- smooth.spline(x, vals)$y # smoothness determined using CV
    }
    pred[[qnames[i]]] <- vals
  }
  return(pred)
}

