mg1.mcmc <- function(params0, y, opt, n.pred=0) {
  # MCMC sampler for M/G/1 based on the arXiv paper by A. Y. Shestopaloff and 
  # R. M. Neal (2014).
  
  do.rate.upd <- !is.null(opt[["c.rate"]]) && opt$c.rate != 1
  check.constr <- FALSE # for debugging
  
  ns <- opt$n.samples
  nb <- opt$burnin
  nf <- opt$n.final
  n <- length(y)
  if (n*ns > 10^8) {
    stop('Too much memory needed.') # all sampled v currently saved to memory
  }
  x <- cumsum(y)
  C.L <- t(chol(opt$C))
  
  # initialize parameters
  thetas <- matrix(NA,3,ns)
  vs <- matrix(NA,n,ns)
  thetas[,1] <- params0$theta
  vs[,1] <- params0$v
  
  # Run sampler
  for (j in 2:ns) {
    ## 1/3: sample v's
    ##################
    rj <- runif(n)
    for (i in 1:n) {
      # i=1 case:
      if (i == 1) {
        a <- max(0,x[1] - thetas[2,j-1])
        b <- min(vs[2,j-1], x[1] - thetas[1,j-1])
        vs[1,j] <- runifbounds(a,b,rj[i],i,j)
      }
      
      # 1 < i < n case:
      if (1 < i && i < n) {
        b <- min(vs[i+1,j-1], x[i] - thetas[1,j-1]) # same for both cases
        if (y[i] > thetas[2,j-1]) {
          a <- x[i] - thetas[2,j-1]
        } else {
          a <- vs[i-1,j] # note that we have already sampled this
        }
        vs[i,j] <- runifbounds(a,b,rj[i],i,j)
      }
      
      # i=n case:
      if (i == n) {
        U <- x[n] - thetas[1,j-1] # same for both cases
        if (y[n] > thetas[2,j-1]) {
          L <- x[n] - thetas[2,j-1]
        } else {
          L <- vs[n-1,j] # note that we have already sampled this
        }
        z <- rj[i]
        # NOTE: division by (-thetas[3,j-1]) is missing in Eq. 12 of the 2014 arXiv paper
        vs[n,j] <- log((1-z)*exp(-thetas[3,j-1]*L) + z*exp(-thetas[3,j-1]*U)) / (-thetas[3,j-1])
      }
    }
    
    if (check.constr && !is.all.constraint.ok(thetas[,j-1], y, x, vs[,j])) {
      print('1: constraints not ok!')
    }
    
    ## 2/3 sample theta using M-H
    #############################
    # TODO: It would be possible to do many M-H updates in a row with little 
    # additional computational cost but this is not implemented below.
    ## UPDATE ALL THREE THETA-PARAMETERS USING M-H
    if (0) {
      theta.prop <- thetas[,j-1] + C.L%*%rnorm(3)
      if (!is.theta.constraint.ok(theta.prop, y, x, vs[,j])) {
        # reject
        thetas[,j] <- thetas[,j-1]
      } else {
        # M-H check
        logpdf_old <- n*log(thetas[3,j-1]) - thetas[3,j-1]*vs[n,j] - n*log(thetas[2,j-1]-thetas[1,j-1])
        logpdf_prop <- n*log(theta.prop[3]) - theta.prop[3]*vs[n,j] - n*log(theta.prop[2]-theta.prop[1])
        alpha <- logpdf_prop - logpdf_old
        rmh <- log(runif(1))
        if (!is.finite(logpdf_old)) {
          print(logpdf_old)
          stop('logpdf_old not finite!') # should not happen
        }  
        if (!is.finite(logpdf_prop) || alpha < rmh) {
          # reject
          thetas[,j] <- thetas[,j-1]
        } else {
          # accept
          thetas[,j] <- theta.prop
        }
      }
    } else {
      # USE ETA PARAMETRIZATION (as in the 2014 arXiv paper):
      eta.prop <- to.eta(thetas[,j-1]) + C.L%*%rnorm(3)
      theta.prop <- to.theta(eta.prop)
      if (!is.theta.constraint.ok(theta.prop, y, x, vs[,j])) {
        # reject
        thetas[,j] <- thetas[,j-1]
      } else {
        # M-H check
        logpdf_old <- n*log(thetas[3,j-1]) - thetas[3,j-1]*vs[n,j] - n*log(thetas[2,j-1]-thetas[1,j-1]) + log(thetas[3,j-1])
        logpdf_prop <- n*log(theta.prop[3]) - theta.prop[3]*vs[n,j] - n*log(theta.prop[2]-theta.prop[1]) + eta.prop[3]
        alpha <- logpdf_prop - logpdf_old
        rmh <- log(runif(1))
        if (!is.finite(logpdf_old)) {
          print(logpdf_old)
          stop('logpdf_old not finite!') # should not happen
        }  
        if (!is.finite(logpdf_prop) || alpha < rmh) {
          # reject
          thetas[,j] <- thetas[,j-1]
        } else {
          # accept
          thetas[,j] <- theta.prop
        }
      }
    }
    
    if (check.constr && !is.all.constraint.ok(thetas[,j], y, x, vs[,j])) {
      print('2: constraints not ok!')
    }
    
    ## 3/3 Rate scale update
    ########################
    if (do.rate.upd) {
      if (j<=2) {cat('Note: Rate scale updates enabled!\n')}
      #c.rate <- 1.7 # scenario 1
      #c.rate <- 1.004 # scenario 2
      #c.rate <- opt$c.rate
      c.rate <- opt$c.rate*runif(1) # random c.rate
      
      z <- runif(1)
      z <- (z<=.5) - (z>.5) # 1 or -1 with equal prob.
      v.prop <- c.rate^z*vs[,j]
      eta3.prop <- log(thetas[3,j]) - z*log(c.rate)
      theta.prop <- c(thetas[1:2,j], exp(eta3.prop))
      if (is.rate.upd.constraint.ok(v.prop, y, x, theta.prop)) {
        logpdf_old <- n*log(thetas[3,j]) - thetas[3,j]*vs[n,j]  + log(thetas[3,j])
        logpdf_prop <- n*log(theta.prop[3]) - theta.prop[3]*v.prop + eta3.prop
        alpha <- logpdf_prop - logpdf_old + (z*n)*log(c.rate)
        rmh <- log(runif(1))
        if (!is.finite(logpdf_old)) {
          print(logpdf_old)
          stop('logpdf_old not finite in rate upd!') # should not happen
        }
        if (!is.finite(logpdf_prop) || alpha < rmh) {
          # reject, nothing to do
        } else {
          # accept
          thetas[3,j] <- theta.prop[3]
          vs[,j] <- v.prop
        }
      } else {
        # reject, nothing to do
      }
    }
    
    if (check.constr && !is.all.constraint.ok(thetas[,j], y, x, vs[,j])) {
      print('3: constraints not ok!')
    }
  }
  
  ## thin
  thetas <- thin.samples(thetas, nf, nb)
  vs <- thin.samples(vs, nf, nb)
  
  ## sample future data
  ypreds <- NULL; xpreds <- NULL; vpreds <- NULL
  if (n.pred > 0) {
    nff <- dim(thetas)[2]
    ypreds <- matrix(NA,n.pred,nff)
    xpreds <- matrix(NA,n.pred,nff)
    vpreds <- matrix(NA,n.pred,nff)
    for (j in 1:nff) {
      preds <- mg1.simul.future(vs[n,j], x[n], thetas[,j], n.pred)
      ypreds[,j] <- preds$ypreds
      xpreds[,j] <- preds$xpreds
      vpreds[,j] <- preds$vpreds
    }
  }
  return(list(thetas=thetas, vs=vs, ypreds=ypreds, xpreds=xpreds, vpreds=vpreds))
}

################################################################################

to.eta <- function(theta) c(theta[1],theta[2]-theta[1],log(theta[3]))

to.theta <- function(eta) c(eta[1],eta[1]+eta[2],exp(eta[3]))

runifbounds <- function(a,b,rij,i,j) {
  # Generates random numbers from Uniform([a,b]) and also checks the bounds a and b.
  if (0<=a && a<=b) {
    return(a + (b-a)*rij)
  } else {
    #print(i); print(j)
    if (a < 0) {
      print(a)
      stop('Something is wrong, incorrect bounds (a<0) for runif.')
    } else {
      print(a); print(b)
      stop('Something is wrong, incorrect bounds (a>b) for runif.')
    }
  }
}

is.theta.prior.constraint.ok <- function(theta) {
  # Returns true if theta satisfies the prior constraints and false otherwise.
  return(theta[1] >= 0 && theta[1] <= 10 && 
           theta[2] - theta[1] >= 0 && theta[2] - theta[1] <= 10 && 
           theta[3] >= 0 && theta[3] <= 1/3)
}

is.theta.constraint.ok <- function(theta, y, x, v) {
  # Returns true if all the conditions related to theta are satisfied. 
  n <- length(y)
  c0 <- is.theta.prior.constraint.ok(theta)
  #c1 <- (theta[2] - theta[1] >= 0)
  c2 <- (theta[1] <= y[1] - v[1])
  c3 <- (theta[2] >= y[1] - v[1])
  m <- y[2:n] - pmax(0, v[2:n]-x[1:(n-1)])
  c4 <- (theta[1] <= min(m))
  c5 <- (theta[2] >= max(m))
  return(all(c(c0, c2,c3,c4,c5)))
}

is.all.constraint.ok <- function(theta, y, x, v) {
  # Checks that the given parameters satisfy all constraints.
  cth <- is.theta.constraint.ok(theta, y, x, v)
  cv1 <- v[1]>=0 #all(v>=0)
  cv2 <- all(diff(v)>=0)
  cv3 <- all(x>=v)
  return(all(c(cth,cv1,cv2,cv3)))
}

is.rate.upd.constraint.ok <- function(v.pr, y, x, theta.pr) {
  # Returns true if the constraints related to the rate upd are satisfied.
  # Alternatively, 'is.all.constraint.ok' could be used but it includes some 
  # redundant checkings (e.g. rate upd cannot make v_1 negative).
  ind <- y>theta.pr[2]
  return(theta.pr[3] >= 0 && theta.pr[3] <= 1/3 && all(v.pr<x-theta.pr[1]) && 
           all(v.pr[ind]>x[ind]-theta.pr[2]))
}

################################################################################

init.theta.v <- function(y, x, theta0) {
  # Rough initialization. Not very useful...
  
  n <- length(y)
  v <- x - min(y)
  theta <- c(0,0,1/6)
  m <- y[2:n] - pmax(0, v[2:n]-x[1:(n-1)])
  theta[1] <- max(0,min(y[1]-v[1],min(m))-0.1)
  theta[2] <- max(y[1]-v[1],max(m))+0.1
  
  # # Initialization mentioned in the arXiv 2014 paper p.10:
  # # It does not seem to work, maybe I'm missing something here?
  # my <- min(y)
  # theta <- to.theta(c(my,5-my,-2.089)) #-2.089==mean(log(1/3*runif(n=10000)))
  # v <- x - my
  
  params <- list(theta = theta, v = v)
  return(params)
}


