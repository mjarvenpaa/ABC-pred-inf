simplemarkov1D.inference <- function() {
  # 1D version, only c considered as a parameter to be estimated.
  # This code is used to generate the numerical results for the first Example 
  # of the paper, in **Figure 1** and **Figure C.1** in appendix. 
  # TODO: Could use other approach for weighting the summaries in 'ABC2' case. 
  #
  # This function compares the following methods:
  # - exact predictive density based on the true parameter c
  # - exact (Gaussian) posterior for c and for prediction
  # - ABC1: 1D summary (bary) [sufficient only for inferring c]
  # - ABC2: 2D summary (bary,y_n) [sufficient for all inference tasks here]
  # - ABC3: 1D summary (bary+phi*y_n) [sufficient only for the first step prediction]
  # - ABCF: (any of the above summaries)
  
  library('matrixStats')
  
  seed.data <- 123456
  seed.inf <- 123456
  
  inf.task <- 0
  plot.task <- 1
  
  opt <- list()
  opt$scenario <- 1 # ***initial phi=0.5 case in the main paper***
  ###opt$scenario <- 2 # ***second example; phi=0.99 case in the appendix***
  opt$root <- '/home/mj/work/ABC_PRED/simulations/code' # code location
  opt$save.loc <- file.path(opt$root,'..','results','simplemarkov1D') # where output and plottings saved
  opt$fn.samples <- file.path(opt$save.loc, paste0('abc_data','_sce',opt$scenario,'.RData'))
  
  source(file.path(opt$root,'simplemarkov_model.R'))
  source(file.path(opt$root,'simplemarkov_inference.R'))
  source(file.path(opt$root,'abc_mcmc.R'))
  source(file.path(opt$root,'abc_distance.R'))
  source(file.path(opt$root,'help_functions.R'))
  dyn.load(file.path(opt$root,"simplemarkov_simul.so"))
  
  n <- 100 # obs.data size
  n.pred <- 100 # prediction size
  theta.true <- 1 # true value of the parameter c to be estimated
  if (opt$scenario==1) {
    phi <- 0.5 # fixed phi
  } else {
    phi <- 0.99 # fixed phi
  }
  s2 <- 1^2 # fixed sigma^2
  q <- 0.10 # which CI level to the plot, e.g. 0.10==90% central CI (not actually used here)
  pp.fig <- c(1, 10, 100) # where to predict in the figure
  
  if (opt$scenario==1) {
    epss <- c(0.037, 0.9, 0.23) # NOTE: different threshold for each ABC method!
  } else {
    epss <- c(0.035, 1.7, 2.25) # NOTE: different threshold for each ABC method!
  }
  opt.abc <- list()
  opt.abc$n.samples <- 10^6
  opt.abc$burnin <- 10^4
  opt.abc$n.final <- 10^4
  opt.abc$theta0 <- theta.true
  opt.abc$cov.adapt <- T
  opt.abc$cov.adapt.len <- opt.abc$burnin
  opt.abc$cov.adapt.t0 <- 5*10^2
  opt.abc$cov.adapt.e <- 10^-6
  opt.abc$C <- 0.1^2
  opt.abc$d.n <- 500 # how many simulations for scaling of summaries (when summary stats dimension>1)
  opt.abc$d.cov.method <- 'cov' # 'cov', 'stdev', 'mad'
  opt.abc$d.sim.method <- 'th' # 'th', 'prior'
  
  abc.rep <- 3 # how many ABC summary stats cases
  pp <- 1:n.pred
  pp.fig <- unique(pmin(n.pred,pp.fig)) 
  
  # generate 'observed' data from the model:
  set.seed(seed.data)
  nn <- n + n.pred
  y <- simplemarkov.simul(c(theta.true,phi,s2), nn)
  y.obs <- y[1:n]
  #y.pred.obs <- y[(n+1):nn]
  
  if (inf.task) {
    # infer parameters using ABC-MCMC and with different summary stats
    ##################################################################
    # set model
    sim.model <- function(theta) simplemarkov.simul(c(theta,phi,s2), nn)
    
    # set prior
    prior.lb <- -10
    prior.ub <- 10
    log.prior <- function(theta) logpmunif(theta, prior.lb, prior.ub)
    prior.sim <- function() rmunif(1, prior.lb, prior.ub)
    
    res.abc <- vector('list',abc.rep)
    for (i in 1:abc.rep) {
      # set summary statistics and discrepancy:
      if (i==1) {
        sumstats <- function(y) {
          return(((1-phi)*sum(y[1:(n-1)])+y[n])/n) # (bary)
        }
      } else if (i==2) {
        sumstats <- function(y) {
          return(c(((1-phi)*sum(y[1:(n-1)])+y[n])/n, y[n])) # (bary,y_n)
        }
      } else {
        sumstats <- function(y) {
          return(((1-phi)*sum(y[1:(n-1)])+y[n])/n + phi*y[n]) # (bary+phi*y_n)
        }
      }
      sumstats.obs <- sumstats(y.obs)
      sdim <- length(sumstats.obs)
      if (sdim==1) {
        invW <- 1 # weighting procedure unnecessary when only one summary
      } else {
        invW <- abc.dist.invW(sdim, sim.model, sumstats, opt.abc$d.n, opt.abc$d.cov.method, 
                              opt.abc$d.sim.method, th=theta.true, prior.sim, pr=T)
        if (opt$scenario==2) {
          invW[sdim,sdim] <- 1.5*invW[sdim,sdim] # some more weight to y_n-summary!
        }
      }
      
      discr <- function(y.sim) {
        d <- sumstats(y.sim) - sumstats.obs
        return(sqrt(d%*%invW%*%d))
      }
      
      # adjustments for the current ABC method
      opt.abc$eps <- epss[i]
      
      # Run ABC-MCMC!
      cat(paste0('Running ABC, ',i,'/',abc.rep,'...\n'))
      set.seed(seed.inf)
      res.abc[[i]] <- simple.abcmcmc(opt.abc$theta0, sim.model, log.prior, discr, opt.abc)
      cat('ABC done!\n')
    }
    save(res.abc, file = opt$fn.samples)
  }
  
  if (plot.task) {
    load(file = opt$fn.samples)
    set.seed(seed.data)
    
    graphics.off()
    for (i in 1:abc.rep) {
      cat(paste0('ABC, ',i,':\n'))
      basic.abc.mcmc.check(res.abc[[i]])
      cat('\n')
    }
    
    # exact (Gaussian) posterior for parameter c
    res.cpost <- simplemarkov1D.cpost.stats(phi, s2, y.obs, n)
    
    # predictive distribution based on true parameter c
    y.pred.true <- simplemarkov.pred.post.stats.true(c(theta.true,phi,s2), pp, y.obs, n, q)
    
    # exact posterior predictive distribution
    y.pred.post <- simplemarkov1D.pred.post.stats(phi, s2, pp, y.obs, n, q)
    
    # Compute ABC mean prediction and predictive intervals from the output
    y.pred.abc <- vector('list',abc.rep)
    for (i in 1:abc.rep) {
      y.pred.abc[[i]] <- pred.post.stats(res.abc[[i]]$y.sims, q)
    }
    
    # ABC-F approach
    ###y.pred.abcf <- simplemarkov1D.pred.post.stats.mcmc(phi, s2, samples, pp, y.obs, n, q) # Not used currently
    
    ## PLOT RESULTS
    ###############
    for (i in 1:abc.rep) {
      if (i>1) {dev.new()}
      simple.plot.mcmc.chain(res.abc[[i]]$thetas) # ABC
    }
    
    # plot posterior of parameters and predictions
    simplemarkov1D.plot.params.and.pred(res.cpost, y.pred.true, y.pred.post, res.abc, phi, s2, y.obs, n, pp.fig, opt)
    
    # TODO: could also plot mean, quantiles etc. as a function of time for additional illustration
    #...
  }
  invisible()
}

################################################################################

simplemarkov1D.cpost.stats <- function(phi, s2, y.obs, n) {
  # Computes mean and sd of the exact (Gaussian) posterior for parameter c.
  
  res <- list(mean=((1-phi)*sum(y.obs[1:(n-1)])+y.obs[n])/n, sd=sqrt(s2/n))
  return(res)
}

simplemarkov1D.pred.post.stats <- function(phi, s2, pp, y.obs, n, q) {
  # Computes mean, sd and q-quantiles of the exact posterior predictive distribution.
  
  bary <- ((1-phi)*sum(y.obs[1:(n-1)])+y.obs[n])/n
  pred <- list()
  pred$mean <- (1-phi^pp)/(1-phi)*bary + phi^pp*y.obs[n]
  pred$sd <- sqrt(pmax(0,s2*(1-phi^(2*pp))/(1-phi^2) + (1-phi^pp)^2/(1-phi)^2*s2/n))
  pred$u <- pred$true.mean + qnorm(1-q/2)*pred$true.sd
  pred$l <- pred$true.mean + qnorm(q/2)*pred$true.sd
  return(pred)
}

simplemarkov1D.pred.post.stats.mcmc <- function(phi, s2, samples, pp, y.obs, n, q) {
  # Computes mean and q-quantiles of the predictive distribution based on MCMC 
  # (or ABC) samples (currently no weights). 
  
  cs <- samples
  n.pred <- length(pp)
  pred <- list(mean = rep(0,n.pred), u = rep(0,n.pred), l = rep(0,n.pred))
  for (p in pp) {
    mes <- cs*(1-phi^p)/(1-phi) + phi^p*y.obs[n]
    sds <- sqrt(pmax(0,s2*(1-phi^(2*p))/(1-phi^2)))
    # mean of Gaussian mixture:
    pred$mean[p] <- mean(mes)
    # quantiles of Gaussian mixture:
    pred$u[p] <- gaussian.mixture.quantile(1-q/2,1/length(cs),mes,sds)
    pred$l[p] <- gaussian.mixture.quantile(q/2,1/length(cs),mes,sds)
  }
  return(pred)
}

simplemarkov1D.pred.post.dens.mcmc <- function(xp, phi, s2, samples, pp, y.obs, n) {
  # Computes full predictive distribution at points 'xp' based on MCMC 
  # (or ABC) samples (currently no weights)
  # Returns density values in #xp x #pp matrix.
  
  cs <- samples
  gmp <- matrix(0,dim(xp)[1],length(pp))
  for (i in 1:length(pp)) {
    p <- pp[i]
    mes <- cs*(1-phi^p)/(1-phi) + phi^p*y.obs[n]
    sds <- sqrt(pmax(0,s2*(1-phi^(2*p))/(1-phi^2)))
    for (j in (1:dim(xp)[1])) {
      gmp[j,i] <- mean(dnorm(xp[j,i], mean = mes, sd = sds))
    }
  }
  return(gmp)
}

################################################################################
# Plotting etc. functions:

simplemarkov1D.plot.params.and.pred <- function(res.cpost, y.pred.true, y.pred.post, res.abc, phi, s2, y.obs, n, pp.fig, opt) {
  # Plots posterior densities obtained using ABC and the baselines.
  
  library(latex2exp)
  ins <- c('True param.','Exact post.',TeX(c('ABC, $s^{(1)}$','ABC, $s^{(2)}$','ABC, $s^{(3)}$')))
  cols <- c('black','red','blue','orange','magenta','orange')
  #lts <- c('solid','solid','solid','dashed','dotdash','solid')
  lts <- c('solid','solid','solid','solid','solid','solid')
  #ltys <- c(1,1,1,2,4,1)
  ltys <- rep(1,6)
  ins2 <- c('True param.','Exact post.',TeX(c('ABC-P, $s^{(1)}$','ABC-P, $s^{(2)}$','ABC-P, $s^{(3)}$', 'ABC-F, $s^{(3)}$')))
  pn <- 'c'
  #y.name <- '\\tilde{y}'
  y.name <- 'y'
  lw <- 1.05
  
  d <- 1 + length(pp.fig)
  abc.rep <- length(res.abc)
  kde <- function(samples) density(samples, adjust = 1)
  fn <- file.path(opt$save.loc, paste0('params_pred_plot','_sce',opt$scenario,'.pdf'))
  pdf(file=fn, width = min(8,1.6*d), height = 2)
  par(mfrow=c(1,d))
  par(mai=c(0.45,0.05,0.05,0.05), mgp=c(1.8,0.5,0))
  pabc <- vector('list',abc.rep)
  for (i in 1:d) {
    if (i==1) {
      # parameter c
      xla <- TeX(paste0('$',pn,'$'))
      for (j in 1:abc.rep) {
        pabc[[j]] <- kde(res.abc[[j]]$thetas)
      }
      adx <- c(0,0.1)
      if (opt$scenario==1) {adx <- c(0.25,0.5)} # finetune plot range so that legend fits better
      rax <- c(min(pabc[[1]]$x)-adx[1], max(pabc[[1]]$x)+adx[2])
      xgrid <- seq(rax[1],rax[2],len=1000)
      d.cpost <- dnorm(xgrid, mean=res.cpost$mean, sd=res.cpost$sd)
      ray <- c(0, 1.05*max(d.cpost))
      plot(rep(y.pred.true$theta.true[1],2),lwd=lw,c(ray[1],1.1*ray[2]),col=cols[1],lty=lts[1],type='l',
           ylab = '',xlab = xla, main = '',yaxt='n', xlim=rax, ylim=ray) # true c as vertical line
      #lines(rep(y.pred.true$theta.true[1],2),lwd=lw,c(ray[1],1.1*ray[2]),col=cols[1],lty=lts[1]) # true c as vertical line
      lines(xgrid, d.cpost, col=cols[2],lty=lts[2],lwd=lw)
      for (j in 1:abc.rep) {
        lines(pabc[[j]],col=cols[2+j],lty=lts[2+j],lwd=lw)
      }
      legend(x='topright', inset = c(0.02,0.02), legend=ins[1:5], col=cols[1:5], lty=ltys[1:5], lwd=rep(lw,5), bg = "white", cex=0.65)
      
    } else {
      # prediction
      xla <- TeX(paste0('$',y.name,'_{',n+pp.fig[i-1],'}$'))
      for (j in 1:abc.rep) {
        pabc[[j]] <- kde(res.abc[[j]]$y.sim[n+pp.fig[i-1],])
      }
      adx <- c(0,0)
      if (i==2) {adx <- c(0.7,0.1)} # finetune plot range so that legend fits better
      rax <- c(min(pabc[[1]]$x)-adx[1], max(pabc[[1]]$x)+adx[2])
      xgrid <- seq(rax[1],rax[2],len=1000)
      d.true <- dnorm(xgrid, mean=y.pred.true$mean[pp.fig[i-1]], sd=y.pred.true$sd[pp.fig[i-1]])
      d.post <- dnorm(xgrid, mean=y.pred.post$mean[pp.fig[i-1]], sd=y.pred.post$sd[pp.fig[i-1]])
      ray <- c(0, 1.05*max(d.true, d.post))
      plot(xgrid, d.true,col=cols[1],lty=lts[1],type='l',lwd=lw,ylab = '',xlab = xla, main = '',yaxt='n', xlim=rax, ylim=ray)
      
      # ABC-F: 
      # NOTE: We use the *third* summary, in other cases similar results as the true posterior were obtained
      # and these are not plotted for clarity!
      ###gmp.abcf <- simplemarkov1D.pred.post.dens.mcmc(as.matrix(xgrid), phi, s2, res.abc[[3]]$thetas, pp.fig[i-1], y.obs, n)
      ###lines(xgrid, gmp.abcf, col=cols[6],lty=lts[6])
      
      #lines(xgrid, d.true, col=cols[1],lty=lts[1])
      lines(xgrid, d.post, col=cols[2],lty=lts[2],lwd=lw)
      
      for (j in 1:abc.rep) {
        lines(pabc[[j]], col=cols[2+j],lty=lts[2+j],lwd=lw)
      }
      
      locl <- 'topright'
      if (opt$scenario==1) {locl <- 'topleft'} # finetune legend location
      if (i==2) {
        legend(x=locl, inset = c(0.02,0.02), legend=ins2[1:5], col=cols[1:5], lty=ltys, lwd=rep(lw,5), bg = "white", cex=0.65)
      }
    }
  }
  dev.off()
}

