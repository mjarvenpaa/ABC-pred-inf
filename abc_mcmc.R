simple.abcmcmc <- function(params0, sim.model, log.prior, discr, opt) {
  # A simple ABC-MCMC algorithm implementation that is also compatible with the 
  # ABC prediction. Internal thinning/burn-in removal is implemented to reduce 
  # memory usage and the size of the resulting samples. A Gaussian proposal is 
  # used whose cov. matrix can be adapted as in the Adaptive Metropolis 
  # algorithm by Haario et al. Bernoulli 2001. 
  # It is assumed sim.model function call returns either data output 'y' 
  # directly or a list containing both 'y' and latent variables 'lat'. It is 
  # assumed that 'y' and 'lat' have fixed size wrt. randomness and parameter and
  # that they are either vectors or matrices (not arrays).
  # TODO: Improved the cov.matrix adaption code
  
  ns <- opt$n.samples # how many iterations of ABC-MCMC
  nbi <- opt$burnin # how many samples neglected as "burn-in"
  nf <- opt$n.final # how many samples to be returned
  if (ns <= 0 || nbi < 0 || nf <= 0 || ns-nbi-nf < 0 || ns > 10^8 || 
      nf > 10^5 || (ns-nbi)%%nf != 0) {
    stop('Number of MCMC samples specified incorrectly.')
  }
  by <- (ns-nbi)/nf
  p <- length(params0)
  eps <- opt$eps # ABC threshold
  C.L <- t(chol(opt$C))
  
  sim.output0 <- sim.model(params0) # simulation at initial point
  output.is.list <- is.list(sim.output0) # whether sim.model returns a list or not
  if (output.is.list) {
    y.is.vec <- is.vector(sim.output0$y)
    lat.is.vec <- is.vector(sim.output0$lat)
    nn <- ifelsem(y.is.vec,length(sim.output0$y),dim(sim.output0$y)) # size of output data
    nl <- ifelsem(lat.is.vec,length(sim.output0$lat),dim(sim.output0$lat)) # size of output latent variables
    if (lat.is.vec) {
      lat.sims <- matrix(NA,nl,nf) # latent variables returned from the sim.model
    } else {
      lat.sims <- array(NA,c(nl,nf))
    }
  } else {
    y.is.vec <- is.vector(sim.output0)
    nn <- ifelsem(y.is.vec,length(sim.output0),dim(sim.output0)) # size of output data
    nl <- 0
    lat.sims <- NULL
  }
  if ((prod(nn)+prod(nl))*nf > 10^8) {
    stop('Too much memory needed.')
  }
  thetas <- matrix(NA,p,nf)
  if (y.is.vec) {
    y.sims <- matrix(NA,nn,nf)
  } else {
    y.sims <- array(NA,c(nn,nf))
  }
  theta.cur <- params0
  output.cur <- sim.output0
  
  # for simple covariance matrix adaption:
  cov.adapt <- !is.null(opt[["cov.adapt"]]) && opt$cov.adapt
  if (cov.adapt) {
    thetas.a <- matrix(NA,p,opt$cov.adapt.len)
    thetas.a[,1] <- params0
    sd <- 2.4^2/p
  }
  
  i <- 1
  nr.acc <- 0 # number of acceptances
  nr.acc.nobi <- 0 # number of acceptances after burn-in
  acc <- NA
  for (j in 2:ns) {
    #print(j)
    # adapt the covariance matrix of the proposal:
    if (cov.adapt && j > opt$cov.adapt.t0 && j <= opt$cov.adapt.len) {
      # TODO: Recursive computation or updating only at e.g. every 10th iter.
      C.L <- t(chol(sd*(cov(t(thetas.a[,1:(j-1),drop=F])) + opt$cov.adapt.e*diag(p))))
    }
    theta.prop <- theta.cur + C.L%*%rnorm(p)
    if (!is.finite(log.prior(theta.prop))) {
      # out of bounds -> reject
      acc <- FALSE
    } else {
      # inside bounds, test whether to accept or not based on the simulation
      output.prop <- sim.model(theta.prop)
      discr.prop <- discr(output.prop) # evaluates discrepancy (which can also be vector)
      if (length(discr.prop) != length(eps)) {
        stop('Incorrect ABC thresholding.')
      }
      if (all(discr.prop <= eps) && 
          log.prior(theta.prop)-log.prior(theta.cur) >= log(runif(1))) {
        acc <- TRUE # accept
      } else {
        acc <- FALSE # reject
      }
    }
    if (acc) {
      theta.cur <- theta.prop
      output.cur <- output.prop
      nr.acc <- nr.acc + 1
      if (j > nbi) {
        nr.acc.nobi <- nr.acc.nobi + 1
      }
    } else {
      # theta.cur and output.cur remain at their previous values
    }
    # whether to save the theta/sim.output to final sample set
    if (j - nbi > 0 && (j%%by) == 0) {
      thetas[,i] <- theta.cur
      if (output.is.list) {
        if (y.is.vec) {
          y.sims[,i] <- output.cur$y
        } else {
          y.sims[,,i] <- output.cur$y
        }
        if (lat.is.vec) {
          lat.sims[,i] <- output.cur$lat
        } else {
          lat.sims[,,i] <- output.cur$lat
        }
      } else {
        if (y.is.vec) {
          y.sims[,i] <- output.cur
        } else {
          y.sims[,,i] <- output.cur
        }
      }
      i <- i + 1
    }
    # save current theta for computing adaptive cov:
    if (cov.adapt && j <= opt$cov.adapt.len) {
      thetas.a[,j] <- theta.cur
    }
  }
  
  # handle output
  res <- list(thetas = thetas, y.sims = y.sims, lat.sims = lat.sims)
  res$acc.prob <- nr.acc/ns # acceptance rate
  res$acc.prob.nobi <- nr.acc.nobi/(ns-nbi) # acceptance rate after burn-in
  return(res)
}

################################################################################

thin.samples <- function(samples, n.final, n.burnin=0) {
  # Simple thinning function that removes burn-in and thins to the desired size.
  # samples must be #params x #samples matrix (or vector if #params==1).

  if (is.vector(samples)) {
    samples <- t(as.matrix(samples))
  }
  n.sa <- dim(samples)[2]
  if (n.burnin < 0 || n.final < 0 || n.final > n.sa - n.burnin) {
    stop('Incorrect burnin or final.size.')
  }
  inds <- unique(seq(n.burnin+1, n.sa, len=n.final))
  return(samples[,inds])
}

basic.abc.mcmc.check <- function(res, plt=F) {
  # Some basic convergence checkings using the coda package. Currently based on 
  # only a single chain. 
  # TODO: This should be sufficient for the simple illustrative examples 
  # considered here but for more realistic/complex inference tasks more careful
  # analysis of convergence might be needed. 
  
  library(coda)
  if (!is.null(res[["acc.prob"]])) {
    cat('Acceptance probability:\n')
    print(res$acc.prob)
  }
  if (!is.null(res[["acc.prob.nobi"]])) {
    cat('Acceptance probability (after burn-in):\n')
    print(res$acc.prob.nobi)
  }
  samples <- as.mcmc(t(res$thetas))
  cat('Effective sample sizes:\n')
  print(effectiveSize(samples))
  if (plt) {
    dev.new()
    autocorr.plot(samples, 20)
    #gelman.diag(samples, autoburnin=FALSE, multivariate=TRUE) # requires multiple chains
  }
}

################################################################################

simple.plot.mcmc.chain <- function(samples) {
  # Plots MCMC chains, samples must be #params x #samples matrix (or vector if #params==1).
  
  if (is.vector(samples)) {
    samples <- t(as.matrix(samples))
  }
  p <- dim(samples)[1]
  n.samples <- dim(samples)[2]
  if (p > 30 || n.samples <= 1) {
    stop('Incorrect sample set to simple chain plot.')
  }
  if (p <= 3) {
    par(mfrow=c(p,1))
  } else {
    par(mfrow=c(ceiling(p/3),3))
  }
  inds <- floor(seq(1,n.samples,len=1000))
  for (i in 1:p) {
    plot(inds,samples[i,inds],type='l',xlab = 'sample',ylab = '',xlim = c(1,n.samples))
  }
}

simple.pair.plot <- function(samples) {
  # A (very) simple pair plot for visualizing the MCMC samples.
  
  pairs(t(thin.samples(samples, 1000))) 
}


