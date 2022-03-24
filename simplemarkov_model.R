simplemarkov.simul <- function(theta, n) {
  # Initial value X_0=0 is not returned
  
  cp <- theta[1]
  phi <- theta[2]
  sigma <- sqrt(theta[3])
  
  es <- rnorm(n,0,sigma)
  x <- rep(0,n)
  res=.C("simplemarkov_simul_loop",as.double(x), as.double(es), as.double(cp), as.double(phi), as.integer(n))
  x <- res[[1]]
  return(x)
}

simplemarkov.test_simul <- function() {
  # Test simulation
  
  opt <- list()
  opt$root <- '/home/mj/work/ABC_PRED/simulations/code' # code location
  dyn.load(file.path(opt$root,"simplemarkov_simul.so"))
  
  n <- 50;
  theta <- c(0.5, 0.9, 1^2) # c, phi, sigma^2
  y <- simplemarkov.simul(theta, n)
  
  # plot the timeseries data y
  graphics.off()
  plot(1:n, y, type = 'o')
}

