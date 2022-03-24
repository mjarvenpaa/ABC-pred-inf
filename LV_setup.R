lv.run.experiment <- function(seed.data=NA, seed.inf=NA, scenario=NA) {
  # Single run takes ~30min (~20min for the 'missing data' case 20x)
  
  # which inference experiment to run
  scenario <- 101
  #scenario <- 102
  
  #scenario <- 201
  ###scenario <- 202
  
  #scenario <- 301
  ###scenario <- 302
  
  seed.data <- 123456
  
  seed.inf <- 123456
  
  # get settings
  opt <- lv.setup(scenario, seed.data, seed.inf)
  
  source(file.path(opt$root,'LV_inference.R'))
  source(file.path(opt$root,'LV_model.R'))
  source(file.path(opt$root,'abc_mcmc.R'))
  source(file.path(opt$root,'abc_distance.R'))
  source(file.path(opt$root,'help_functions.R'))
  
  # run the inference and plot the results
  # inf.task: whether to run 1) ABC ['baseline' sumstats] 2) ABC [sumstats for pred.]
  # 3) simulation specific for 'missing data' case (neglected otherwise)
  lv.inference(opt, inf.task=c(1,1,1)-1, plot.task=1)
  
  invisible()
}

lv.setup <- function(scenario, seed.data, seed.inf) {
  # Get settings for an ABC inference task with LV model. 
  
  opt <- list(scenario=scenario)
  opt$root <- '/home/mj/work/ABC_PRED/simulations/code' # code location
  opt$save.loc <- file.path(opt$root,'..','results','LV', paste0('sce_',scenario)) # where output and plottings saved
  opt$fn.ext <- paste0('_sce',scenario,'_sd',seed.data,'_si',seed.inf)
  opt$fn.samples <- rep(NA,2)
  for (a in 1:2) {
    opt$fn.samples[a] <- file.path(opt$save.loc, paste0('abc_data_',a,opt$fn.ext,'.RData'))
  }
  sce <- list(seed.data=seed.data, seed.inf=seed.inf); 
  opt.abc <- vector('list',2); plt <- list()
  
  ## COMMON SETTINGS
  sce$lv.mod <- 'LV3' # 3 or 4 parameter version
  sce$xy0 <- c(100,50) # initial state (assumed known)
  sce$theta.true <- c(1, 0.005, 0.6)
  
  #sce$log.prior.lb <- -Inf
  #sce$log.prior.ub <- Inf
  sce$log.prior.lb <- rep(-6,3) # Unif([-6,2]) priors
  sce$log.prior.ub <- rep(2,3)
  
  ##############################################################################
  if (scenario == 101 || scenario == 102) {
    ## LV3, prediction, x&y both observed, no noise
    
    # MODEL AND INFERENCE SETTINGS
    sce$xy.ind <- as.logical(c(1,1)) # whether prey and/or predator observed
    sce$obs.noise.stdev <- 0 # stdev of Gaussian observation noise (assumed known)
    sce$t <- 24 # simulation time with observations
    sce$t.pred <- 21 # time for prediction/missing data scenario
    sce$t2 <- 0 # additional simulation time with observations, 0 if usual prediction task
    sce$dt <- 0.3 # timestep

    # ABC INFERENCE SETTINGS
    # original summary case:
    #opt.abc[[1]]$n.samples <- 0.5*10^5 # TEST
    opt.abc[[1]]$n.samples <- 10^6 # FINAL
    opt.abc[[1]]$C <- diag(c(0.05, 0.05, 0.1)^2)
    opt.abc[[1]]$burnin <- 10^4
    opt.abc[[1]]$n.final <- 10^4
    opt.abc[[1]]$theta0 <- sce$theta.true
    opt.abc[[1]]$cov.adapt <- T
    opt.abc[[1]]$cov.adapt.len <- opt.abc[[1]]$burnin
    opt.abc[[1]]$cov.adapt.t0 <- 5*10^2
    opt.abc[[1]]$cov.adapt.e <- 10^-6
    
    opt.abc[[1]]$summary <- 'orig'
    opt.abc[[1]]$d.n <- 1000 # for scaling of summaries
    opt.abc[[1]]$d.cov.method <- 'mad' # 'cov', 'stdev', 'mad'
    if (scenario == 101) {
      opt.abc[[1]]$d.sim.method <- 'th' # 'th', 'prior'
    } else { # 102 case:
      opt.abc[[1]]$d.sim.method <- 'prior' # did not work well!!
    }
    
    # pred summary case:
    opt.abc[[2]] <- opt.abc[[1]]
    opt.abc[[2]]$C <- diag(c(0.05, 0.05, 0.1)^2)
    opt.abc[[2]]$summary <- 'pred'
    
    # manually adjust suitable ABC threshold
    # acc prob ~5%
    if (scenario == 101) {
      if (seed.data == 123456) {
        opt.abc[[1]]$eps <- 2.0
        opt.abc[[2]]$eps <- c(2.5, 100)
      } else if (seed.data == 123456+1) {
        opt.abc[[1]]$eps <- 3.5
        opt.abc[[2]]$eps <- c(4.5, 100)
      } else if (seed.data == 123456+2) {
        opt.abc[[1]]$eps <- 2
        opt.abc[[2]]$eps <- c(3.5, 100)
      } else if (seed.data == 123456+3) { # acc prob only ~1% in this special case
        opt.abc[[1]]$eps <- 6
        opt.abc[[2]]$eps <- c(15, 200)
      }
      
    #*******************************
    } else { # 102 case:
      if (seed.data == 123456) { # acc prob only ~2%, this approach did not work well
        # This is likely because the behaviour of the model is very different with
        # different parameters of the prior support.
        opt.abc[[1]]$eps <- 10000
        opt.abc[[2]]$eps <- c(10000, 100)
      }
    }
    
    ############################################################################
  } else if (scenario == 201 || scenario == 202) {
    ## LV3, prediction, **missing data** case, x&y observed, no noise
    
    # MODEL AND INFERENCE SETTINGS
    sce$xy.ind <- as.logical(c(1,1)) # whether prey and/or predator observed
    sce$obs.noise.stdev <- 0 # stdev of Gaussian observation noise (assumed known)
    sce$t <- 15 # simulation time with observations
    sce$t.pred <- 21 # time for prediction/missing data scenario
    sce$t2 <- 15 # additional simulation time with observations, 0 if usual prediction task
    sce$dt <- 0.3 # timestep
    
    # settings for computing the approximate true result (specific to this case)
    sce$mis.obs <- list(n.samples = 2*10^5, n.final = 2*10^3)
    sce$mis.obs$fn.samples <- file.path(opt$save.loc, paste0('abc_mis_approx_true',opt$fn.ext,'.RData'))

    # ABC INFERENCE SETTINGS
    #opt.abc[[2]]$n.samples <- 0.5*10^5 # TEST
    opt.abc[[2]]$n.samples <- 10^6 # FINAL
    opt.abc[[2]]$C <- diag(c(0.05, 0.05, 0.1)^2)
    opt.abc[[2]]$burnin <- 10^4
    opt.abc[[2]]$n.final <- 10^4
    opt.abc[[2]]$theta0 <- sce$theta.true
    opt.abc[[2]]$cov.adapt <- T
    opt.abc[[2]]$cov.adapt.len <- opt.abc[[2]]$burnin
    opt.abc[[2]]$cov.adapt.t0 <- 5*10^2
    opt.abc[[2]]$cov.adapt.e <- 10^-6
    
    opt.abc[[2]]$summary <- 'misdata'
    opt.abc[[2]]$d.n <- 1000 # for scaling of summaries
    opt.abc[[2]]$d.cov.method <- 'mad' # 'cov', 'stdev', 'mad'
    if (scenario == 201) {
      opt.abc[[2]]$d.sim.method <- 'th' # 'th', 'prior'
    } else { # 102 case:
      opt.abc[[2]]$d.sim.method <- 'prior'
    }
    
    # manually adjust suitable ABC threshold
    if (scenario == 201) {
      if (seed.data == 123456) { # acc prob ~2%
        opt.abc[[2]]$eps <- c(2.5, 50)
      } else if (seed.data == 123456+1) { # acc prob ~1%
        opt.abc[[2]]$n.samples <- 10^6
        opt.abc[[2]]$eps <- c(2.5, 50)
      }
      
      #*******************************
    } else { # 202 case:
      # ...
    }
    
    ############################################################################
  } else if (scenario == 301 || scenario == 302) {
    ## LV3, prediction, **only x** observed, no noise
    
    # MODEL AND INFERENCE SETTINGS
    sce$xy.ind <- as.logical(c(1,0)) # whether prey and/or predator observed
    sce$obs.noise.stdev <- 0 # stdev of Gaussian observation noise (assumed known)
    sce$t <- 24 # simulation time with observations
    sce$t.pred <- 21 # time for prediction/missing data scenario
    sce$t2 <- 0 # additional simulation time with observations, 0 if usual prediction task
    sce$dt <- 0.3 # timestep
    
    # ABC INFERENCE SETTINGS
    # original summary case:
    #opt.abc[[1]]$n.samples <- 0.5*10^5 # TEST
    opt.abc[[1]]$n.samples <- 10^6 # FINAL
    opt.abc[[1]]$C <- diag(c(0.1, 0.1, 0.1)^2)
    opt.abc[[1]]$burnin <- 10^4
    opt.abc[[1]]$n.final <- 10^4
    opt.abc[[1]]$theta0 <- sce$theta.true
    opt.abc[[1]]$cov.adapt <- T
    opt.abc[[1]]$cov.adapt.len <- opt.abc[[1]]$burnin
    opt.abc[[1]]$cov.adapt.t0 <- 5*10^2
    opt.abc[[1]]$cov.adapt.e <- 10^-6
    
    opt.abc[[1]]$summary <- 'orig'
    opt.abc[[1]]$d.n <- 1000 # for scaling of summaries
    opt.abc[[1]]$d.cov.method <- 'mad' # 'cov', 'stdev', 'mad'
    if (scenario == 301) {
      opt.abc[[1]]$d.sim.method <- 'th' # 'th', 'prior'
    } else { # 302 case:
      opt.abc[[1]]$d.sim.method <- 'prior'
    }
    
    # pred summary case:
    opt.abc[[2]] <- opt.abc[[1]]
    opt.abc[[2]]$C <- diag(c(0.1, 0.1, 0.1)^2)
    opt.abc[[2]]$summary <- 'pred'
    
    # manually adjust suitable ABC threshold
    # acc prob ~5%
    if (scenario == 301) {
      if (seed.data == 123456) {
        opt.abc[[1]]$eps <- 1.1
        opt.abc[[2]]$eps <- c(1.4, 25)
      } else if (seed.data == 123456+1) {
        #...
      }
      
      #*******************************
    } else { # 302 case:
      #...
    }
  }
  
  ## COMMON PLOTTING SETTINGS
  plt$q <- c(0.10) # which CI level(s) to the plot, e.g. 0.10==90% central CI
  plt$pr <- T # whether to print some info on screen in inf task
  plt$noisypred <- F # TBD
  
  return(c(opt, list(sce=sce, opt.abc=opt.abc, plt=plt)))
}


# LV 4
#sce$theta.true <- c(1, 0.01, 0.01, 0.5)
#opt.abc$C <- diag(c(0.2, 0.01, 0.01, 0.1)^2)
#sce$log.prior.lb <- rep(-6,4) # Unif([-6,2]) priors
#sce$log.prior.ub <- rep(2,4)

