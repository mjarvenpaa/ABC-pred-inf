abc.dist.invW <- function(sdim, sim.model, sstats, n=1000, cov.method='cov', 
                          sim.method='th', th=NA, prior.sim=NA, pr=T, sdim.not=0) {
  # Returns an inverse covariance matrix (or some other weighting matrix which 
  # may be diagonal) invW for Mahalanobis-type of distance function 
  # sqrt((s-s')*invW*(s-s')) to be used as ABC discrepancy.
  #
  # cov.method: Which type of covariance (or "covariance") matrix is computed:
  # - 'simple': identity matrix is selected producing usual Euclidean distance
  # - 'cov': full covariance matrix is estimated
  # - 'stdev': only standard deviations of indiv.summaries which lead to a diagonal matrix
  # - 'mad': as above but MAD used instead of stdev for robustness
  #
  # sim.method: How the simulations for computing the cov matrix are obtained:
  # - 'th': The provided parameter value (point estimate) 'th' is used
  # - 'prior': Sampled from prior predictive using the input function 'prior.sim'
  #
  # sdim.not: specifies how many of the summaries (counted from the end) are 
  # not used for computing the cov matrix, this allows special weighting 
  # treatment for some of the summaries
  
  if (sdim.not>=sdim) {
    stop('It seems summary dimensions are incorrect.')
  }
  if (pr) {
    cat('Summary stats dimension:\n')
    print(sdim)
    if (sdim.not>0) {
      cat('Number of special summaries:\n')
      print(sdim.not)
    }
  }
  if (cov.method=='simple') {
    invW <- diag(sdim-sdim.not)
  } else {
    # need to simulate summary statistics to e.g. estimate the cov matrix
    ss.sims <- matrix(NA,n,sdim-sdim.not)
    if (sim.method=='prior' || sim.method=='th') {
      for (i in 1:n) {
        if (sim.method=='prior') {
          th <- prior.sim() # new sample from the prior
        } 
        ss.sims[i,] <- sstats(sim.model(th))[1:(sdim-sdim.not)]
      }
    } else {
      stop('Incorrect simulation method.')
    }
    if (any(cov.method==c('mad','MAD'))) {
      invW <- diag(1/apply(ss.sims,2,mad)^2) # 1/MAD^2 values on diagonals
    } else if (cov.method=='cov') {
      W <- cov(ss.sims) # could also use robust cov matrix estimation methods
      invW <- solve(W)
    } else if (cov.method=='stdev') {
      invW <- diag(1/diag(cov(ss.sims))) # Diagonal entries i.e. variances of full cov matrix
    } else {
      stop('Incorrect method for computing the weight matrix.')
    }
  }
  if (pr) {
    cat('Distance weight matrix:\n')
    print(invW)
  }
  return(invW)
}

