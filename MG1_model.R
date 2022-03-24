mg1.simul <- function(theta, n) {
  # Simulates data from M/G/1 model.
  
  u <- theta[1] + runif(n)*(theta[2]-theta[1])
  w <- rexp(n, rate = theta[3])
  v <- cumsum(w)
  y <- rep(0,n)
  x <- rep(0,n)
  y[1] <- u[1] + w[1]
  x[1] <- y[1]
  res=.C("MG1_simul_loop",as.double(x), as.double(y), as.double(v), as.double(u), as.integer(n))
  res <- list(x=res[[1]], y=res[[2]], v=res[[3]])
  return(res)
}

mg1.simul.future <- function(v.old, x.old, theta, n.pred) {
  # Simulates *future* data from M/G/1 model given the past values. 
  
  u <- theta[1] + runif(n.pred)*(theta[2]-theta[1])
  w <- rexp(n.pred, rate = theta[3])
  v <- v.old[length(v.old)] + cumsum(w)
  y <- rep(0,n.pred)
  x <- rep(0,n.pred)
  y[1] <- u[1] + max(0, v[1] - x.old[length(x.old)])
  x[1] <- x.old[length(x.old)] + y[1]
  res=.C("MG1_simul_loop",as.double(x), as.double(y), as.double(v), as.double(u), as.integer(n.pred))
  pred <- list(xpreds=res[[1]], ypreds=res[[2]], vpreds=res[[3]])
  return(pred)
}

mg1.test.simul <- function(scenario=2) {
  # Test M/G/1 simulation.
  
  opt <- list()
  opt$root <- '/home/mj/work/ABC_PRED/simulations/code' # code location
  source(file.path(opt$root,'MG1_setup.R'))
  dyn.load(file.path(opt$root,"MG1_simul.so"))
  n <- 100;
  
  graphics.off()
  if (scenario==0) {
    # plot just iid Exponential with unit mean for illustration
    y <- rexp(n = n, rate = 1)
    plot(1:n, y, type = 'o', xlab = '', ylab = '', main = 'iid Exp')
  } else {
    theta <- mg1.get.true.theta(scenario)
    # simulate and collect output
    obs <- mg1.simul(theta, n)
    y <- obs$y # interdeparture times
    wts <- obs$x - obs$v # waiting times
    
    # plot the simulation output
    graphics.off()
    par(mfrow=c(2,1))
    plot(1:n, y, type = 'o', xlab = 'interdeparture times y', ylab = '')
    if (1) {
      lines(1:n, rep(median(y),n),type = 'l', col='red')
      lines(1:n, rep(quantile(y,probs = 0.75),n),type = 'l', col='orange')
      lines(1:n, rep(quantile(y,probs = 0.9),n),type = 'l', col='blue')
    }
    plot(1:n, wts, type = 'o', xlab = 'waiting times', ylab = '')
  }
}

