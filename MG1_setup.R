mg1.run.experiment <- function(seed.data=NA, seed.inf=NA, scenario=NA) {
  # Single run takes ~8min
  
  # which inference experiment to run
  scenario <- 101
  #scenario <- 201 
  
  seed.data <- 123456
  
  seed.inf <- 123456
  
  # get settings
  opt <- mg1.setup(scenario, seed.data, seed.inf)
  
  source(file.path(opt$root,'MG1_inference.R'))
  source(file.path(opt$root,'MG1_model.R'))
  source(file.path(opt$root,'MG1_mcmc.R'))
  source(file.path(opt$root,'abc_mcmc.R'))
  source(file.path(opt$root,'abc_distance.R'))
  source(file.path(opt$root,'help_functions.R'))
  
  # run the inference and plot the results
  # inf.task: whether to run 1) MCMC 2) ABC ['baseline' sumstats] 3) ABC [sumstats for pred.]
  mg1.inference(opt, inf.task=c(1,1,1)-1, plot.task=1)
  
  invisible()
}

mg1.setup <- function(scenario, seed.data, seed.inf) {
  # Get settings for an MCMC/ABC inference task with M/G/1 model.
  
  opt <- list(scenario=scenario)
  opt$root <- '/home/mj/work/ABC_PRED/simulations/code' # code location
  opt$save.loc <- file.path(opt$root,'..','results','MG1', paste0('sce_',scenario)) # where output and plottings saved
  opt$fn.ext <- paste0('_sce',scenario,'_sd',seed.data,'_si',seed.inf)
  opt$fn.samples.mcmc <- file.path(opt$save.loc, paste0('mcmc_data_',opt$fn.ext,'.RData'))
  opt$fn.samples <- rep(NA,2)
  for (a in 1:2) {
    opt$fn.samples[a] <- file.path(opt$save.loc, paste0('abc_data_',a,opt$fn.ext,'.RData'))
  }
  sce <- list(seed.data=seed.data, seed.inf=seed.inf)
  opt.mcmc <- list()
  opt.abc <- vector('list',2); plt <- list()
  
  if (scenario == 101 || scenario == 102) { 
    # param.scenario 1 (increasing queue - approx. linear grow)
    
    ## MODEL AND INFERENCE SETTINGS
    sce$param.scenario <- 1
    sce$n <- 100
    sce$n.pred <- 50
    sce$theta.true <- mg1.get.true.theta(sce$param.scenario)
    
    ## MCMC BASELINE SETTINGS
    opt.mcmc$n.samples <- 5*10^5
    opt.mcmc$burnin <- 10^4
    opt.mcmc$n.final <- 10^4
    opt.mcmc$C <- diag(c(0.05, 0.05, 0.2)^2) # proposal cov for MCMC
    
    ## ABC INFERENCE SETTINGS
    # original summary case:
    #opt.abc[[1]]$n.samples <- 0.2*10^6 # TEST
    opt.abc[[1]]$n.samples <- 2*10^6 # FINAL
    opt.abc[[1]]$burnin <- 10^4
    opt.abc[[1]]$n.final <- 10^4
    opt.abc[[1]]$theta0 <- sce$theta.true
    opt.abc[[1]]$cov.adapt <- TRUE
    opt.abc[[1]]$cov.adapt.len <- opt.abc[[1]]$burnin
    opt.abc[[1]]$cov.adapt.t0 <- 10^3
    opt.abc[[1]]$cov.adapt.e <- 10^-6
    opt.abc[[1]]$C <- diag(c(0.05, 0.3, 0.075)^2)
    
    opt.abc[[1]]$summary <- 'orig' # 'orig', 'pred', 'pred2'
    opt.abc[[1]]$d.n <- 1000 # for scaling of summaries
    opt.abc[[1]]$d.cov.method <- 'cov' # 'cov', 'stdev', 'mad'
    if (scenario == 101) {
      opt.abc[[1]]$d.sim.method <- 'th' # 'th', 'prior'
    } else { # 102 case:
      opt.abc[[1]]$d.sim.method <- 'prior'
    }
    
    # pred summary case:
    opt.abc[[2]] <- opt.abc[[1]]
    opt.abc[[2]]$summary <- 'predf' # 'orig', 'pred', 'pred2'
    
    # manually adjust suitable ABC threshold
    # acc prob ~2%
    if (scenario == 101) {
      if (seed.data == 123456) {
        opt.abc[[1]]$eps <- 1.6
        opt.abc[[2]]$eps <- c(7, 5)
      } else if (seed.data == 123456+1) {
        opt.abc[[1]]$eps <- 1.4
        opt.abc[[2]]$eps <- c(2.4, 5)
      } else if (seed.data == 123456+2) {
        opt.abc[[1]]$eps <- 1.3
        opt.abc[[2]]$eps <- c(3, 7)
      } else if (seed.data == 123456+3) {
        opt.abc[[1]]$eps <- 1.9
        opt.abc[[2]]$eps <- c(6, 5)
      } else if (seed.data == 123456+4) {
        opt.abc[[1]]$eps <- 1.8
        opt.abc[[2]]$eps <- c(5, 5)
      }
    }
    
    ############################################################################
  } else if (scenario == 201 || scenario == 202) { 
    # param.scenario 2 (queue variable)

    ## MODEL AND INFERENCE SETTINGS
    sce$param.scenario <- 2
    sce$n <- 100
    sce$n.pred <- 50
    sce$theta.true <- mg1.get.true.theta(sce$param.scenario)
    
    ## MCMC BASELINE SETTINGS
    opt.mcmc$n.samples <- 5*10^5
    opt.mcmc$burnin <- 10^4
    opt.mcmc$n.final <- 10^4
    opt.mcmc$C <- diag(c(0.02, 0.05, 0.05)^2) # proposal cov for MCMC
    
    ## ABC INFERENCE SETTINGS
    # original summary case:
    #opt.abc[[1]]$n.samples <- 0.2*10^6 # TEST
    opt.abc[[1]]$n.samples <- 2*10^6 # FINAL
    opt.abc[[1]]$burnin <- 10^4
    opt.abc[[1]]$n.final <- 10^4
    opt.abc[[1]]$theta0 <- sce$theta.true
    opt.abc[[1]]$cov.adapt <- TRUE
    opt.abc[[1]]$cov.adapt.len <- opt.abc[[1]]$burnin
    opt.abc[[1]]$cov.adapt.t0 <- 10^3
    opt.abc[[1]]$cov.adapt.e <- 10^-6
    opt.abc[[1]]$C <- diag(c(0.02, 0.25, 0.06)^2)
    
    opt.abc[[1]]$summary <- 'orig' # 'orig', 'pred', 'pred2'
    opt.abc[[1]]$d.n <- 1000 # for scaling of summaries
    opt.abc[[1]]$d.cov.method <- 'cov' # 'cov', 'stdev', 'mad'
    if (scenario == 201) {
      opt.abc[[1]]$d.sim.method <- 'th' # 'th', 'prior'
    } else { # 202 case:
      opt.abc[[1]]$d.sim.method <- 'prior'
    }
    
    # pred summary case:
    opt.abc[[2]] <- opt.abc[[1]]
    opt.abc[[2]]$summary <- 'predf' # 'orig', 'pred', 'pred2'
    
    # manually adjust suitable ABC threshold
    # acc prob ~2%
    if (scenario == 201) {
      if (seed.data == 123456) {
        opt.abc[[1]]$eps <- 1.6
        opt.abc[[2]]$eps <- c(2.5, 5)
      } else if (seed.data == 123456+1) {
        opt.abc[[1]]$eps <- 2.2
        opt.abc[[2]]$eps <- c(3.2, 8) # we actually set this larger than 5 
        # since otherwise meeting the 2% acc.prob. seems problematic
      } else if (seed.data == 123456+2) {
        # MCMC seems to have some difficulties with multimodality of 
        # the posterior of \theta_2 in this case
        opt.abc[[1]]$eps <- 1.2
        opt.abc[[2]]$eps <- c(2.5, 8)
      } else if (seed.data == 123456+3) {
        opt.abc[[1]]$eps <- 2
        opt.abc[[2]]$eps <- c(4, 6)
      } else if (seed.data == 123456+4) {
        opt.abc[[1]]$eps <- 1.5
        opt.abc[[2]]$eps <- c(6, 5)
      }
    }
    
    ############################################################################
  }
  
  ## COMMON PLOTTING SETTINGS
  plt$q <- c(0.10,0.25) # which CI level(s) to the plot, e.g. 0.10==90% central CI
  plt$pp.fig <- c(1,10,sce$n.pred) # where to predict in the figure
  plt$pr <- T # whether to print some info on screen in inf task
  
  ## COMPLETE/CHECK CERTAIN SETTINGS
  if (any(sce$param.scenario==c(1,4))) { 
    opt.mcmc$c.rate <- 1.7 # enables scale rate updates
  }
  plt$pp.fig <- unique(pmin(sce$n.pred,plt$pp.fig))
  
  return(c(opt, list(sce=sce, opt.mcmc=opt.mcmc, opt.abc=opt.abc, plt=plt)))
}

################################################################################

mg1.get.true.theta <- function(scenario) {
  # Returns some parameter value used for testing the ABC inference etc.
  
  ths <- list(c(8, 16, 0.15), # case 1; queue full and waiting times grow approx. linearly
              c(4, 7, 0.15),  # case 2; non-increasing and somewhat varying queue
              c(1, 2, 0.01),  # case 3; queue mostly empty and interarrival times mostly determine y's
              c(5, 12, 0.15)) # additional case; like case 1 but slightly more varying queue length
  if (!any(scenario==1:length(ths))) {
    stop('Incorrect test case parameter requested.')
  }
  return(ths[[scenario]])
}

