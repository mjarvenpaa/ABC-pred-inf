lv.simul <- function(th, xy0, tobs, obs.noise.stdev=0) {
  # 3/4-parameter version of the LV model based on the C-code.
  # X: prey, Y: predator
  
  np <- length(th)
  nn <- length(tobs)
  xy <- matrix(0,2,nn)
  x <- rep(0,nn); y <- rep(0,nn)
  res=.C("LV_simul_full", as.integer(x), as.integer(y), as.integer(xy0), 
         as.double(tobs), as.integer(nn), as.double(th), as.integer(np))
  xy[1,] <- res[[1]]; xy[2,] <- res[[2]]
  # add noise to the observations to simulate observation errors
  if (obs.noise.stdev > 0) {
    xyobs <- xy + matrix(rnorm(length(xy), 0, obs.noise.stdev), 2)
    return(list(xy=xy, xyobs=xyobs))
  }
  return(list(xy=NA, xyobs=xy))
}

lv.test.simul <- function() {
  # Test simulation of the 3/4-parameter LV model.
  
  opt <- list()
  opt$root <- '/home/mj/work/ABC_PRED/simulations/code' # code location
  dyn.load(file.path(opt$root,"LV_simul.so"))
  
  lv.mod <- 'LV3'
  #lv.mod <- 'LV4'
  
  t <- 30 # total simulation time
  dt <- 0.2 # timestep
  tobs <- seq(0,t,by=dt) # observation times
  xy0 <- c(100,50) # initial state
  if (lv.mod=='LV3') {
    th <- c(1, 0.005, 0.6)
  } else {
    th <- c(1, 0.01, 0.01, 0.5) # simulation parameter
  }

  # simulate
  obs <- lv.simul(th, xy0, tobs, obs.noise.stdev=10)
  
  # plot the simulation outputs
  graphics.off()
  #lv.simple.data.plot(obs$xy, tobs)
  lv.simple.data.plot(obs$xyobs, tobs)
}

lv.simple.data.plot <- function(xy,t) {
  # Simple plot of LV observations. 
  
  xcol <- 'blue' # prey
  ycol <- 'red' # predator
  plot(t,xy[1,],type='l',col=xcol,xlab='t',ylab='x,y',ylim=range(xy))
  lines(t,xy[2,],type='l',col=ycol)
}

