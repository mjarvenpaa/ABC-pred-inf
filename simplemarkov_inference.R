simplemarkov.inference <- function() {
  # Additional simple Markov model example not used for the paper. 
  
  library('matrixStats')
  
  seed.data <- 123456
  seed.inf <- 123456
  
  opt <- list()
  opt$root <- '/home/mj/work/ABC_PRED/simulations/code' # code location
  opt$save.loc <- file.path(opt$root,'..','results','simplemarkov') # where output and plottings saved
  opt$fn.samples <- file.path(opt$save.loc, 'mcmc_abc_data.RData')
  
  source(file.path(opt$root,'simplemarkov_model.R'))
  source(file.path(opt$root,'abc_mcmc.R'))
  source(file.path(opt$root,'abc_distance.R'))
  source(file.path(opt$root,'help_functions.R'))
  dyn.load(file.path(opt$root,"simplemarkov_simul.so"))
  
  inf.task <- 0
  plot.task <- 1
  
  n <- 200 # obs.data size
  n.pred <- 100 # prediction size
  theta.true <- c(1, 0.99, 1^2) # c, phi, sigma^2
  q <- 0.10 # which CI level to the plot, e.g. 0.10==90% central CI
  pp.fig <- c(1, 10, 100) # where to predict in the figure
  
  opt.mcmc <- list()
  opt.mcmc$n.samples <- 0.5*10^5
  opt.mcmc$theta0 <- theta.true
  
  opt.abc <- list()
  opt.abc$eps <- 1.2
  opt.abc$n.samples <- 1*10^6
  opt.abc$burnin <- 10^4
  opt.abc$n.final <- 10^4
  opt.abc$theta0 <- theta.true
  opt.abc$cov.adapt <- T
  opt.abc$cov.adapt.len <- opt.abc$burnin
  opt.abc$cov.adapt.t0 <- 5*10^2
  opt.abc$cov.adapt.e <- 10^-6
  # use the cov matrix from standard MCMC run for ABC-MCMC
  #opt.abc$C <- diag(c(.1,.01,.1)^2)
  opt.abc$d.n <- 500 # for scaling of summaries
  opt.abc$d.cov.method <- 'cov' # 'cov', 'stdev', 'mad'
  opt.abc$d.sim.method <- 'th' # 'th', 'prior'
  
  pp <- 1:n.pred
  pp.fig <- unique(pmin(n.pred,pp.fig)) 
  
  # generate 'observed' data from the model:
  set.seed(seed.data)
  nn <- n + n.pred
  y <- simplemarkov.simul(theta.true, nn)
  y.obs <- y[1:n]
  y.pred.obs <- y[(n+1):nn]
  
  if (inf.task) {
    # 1/3 infer parameters using MCMC
    #################################
    cat('Running MCMC...\n')
    set.seed(seed.inf)
    res.mcmc <- simplemarkov.mcmc(opt.mcmc$theta0, y.obs, opt.mcmc$n.samples)
    save(res.mcmc, file = opt$fn.samples)
    cat('MCMC done!\n')
    
    # 2/3 infer parameters using ABC-MCMC
    #####################################
    # set model
    sim.model <- function(theta) simplemarkov.simul(theta, nn)
    
    # set prior
    prior.lb <- c(-5,-2,0)
    prior.ub <- c(5,2,5)
    log.prior <- function(theta) logpmunif(theta, prior.lb, prior.ub)
    prior.sim <- function() rmunif(1, prior.lb, prior.ub)
    
    # set summary statistics and discrepancy:
    sumstats <- function(y) {
      y <- y[1:n]
      #return(c(sum(y)/n, sum(y^2)/n, sum(y[2:n]*y[1:(n-1)])/n))
      #return(c(sum(y)/n, y[n]))
      #return(c(y[n]))
      return(c(sum(y)/n, sum(y^2)/n, sum(y[2:n]*y[1:(n-1)])/n, y[n]))
    }
    
    sumstats.obs <- sumstats(y.obs)
    sdim <- length(sumstats.obs)
    invW <- abc.dist.invW(sdim, sim.model, sumstats, opt.abc$d.n, opt.abc$d.cov.method, 
                          opt.abc$d.sim.method, th=theta.true, prior.sim, pr=T)
    #invW[sdim,sdim] <- 1.1*invW[sdim,sdim] # some more weight to y_n-summary!
    
    discr <- function(y.sim) {
      d <- sumstats(y.sim) - sumstats.obs
      return(sqrt(d%*%invW%*%d))
    }
    
    opt.abc$C <- cov(t(res.mcmc$thetas))
    
    # Run ABC-MCMC!
    cat('Running ABC...\n')
    res.abc <- simple.abcmcmc(opt.abc$theta0, sim.model, log.prior, discr, opt.abc)
    cat('ABC done!\n')
    save(res.mcmc, res.abc, file = opt$fn.samples)
  }
  
  if (plot.task) {
    load(file = opt$fn.samples)
    set.seed(seed.data)
    
    graphics.off()
    basic.abc.mcmc.check(res.mcmc)
    basic.abc.mcmc.check(res.abc)
    
    # posterior predictive (MCMC)
    y.pred.mcmc <- simplemarkov.pred.post.stats.mcmc(res.mcmc$thetas, pp, y.obs, n, q)
    
    # posterior predictive based on true parameters
    y.pred.true <- simplemarkov.pred.post.stats.true(theta.true, pp, y.obs, n, q)
    
    # posterior predictive (ABC)
    y.pred.abc <- pred.post.stats(res.abc$y.sims, q)
    
    # 3/3 ABC with exact predictive lik.
    ####################################
    # We use the same ABC posterior as above -> No separate parameter fitting or plotting
    y.pred.abcf <- simplemarkov.pred.post.stats.mcmc(res.abc$thetas, pp, y.obs, n, q)
    
    ## PLOT RESULTS
    ###############
    # plot MCMC chains
    simple.plot.mcmc.chain(res.mcmc$thetas) # MCMC
    dev.new()
    simple.plot.mcmc.chain(res.abc$thetas) # ABC
    
    # plot posterior of parameters
    simplemarkov.plot.params(res.mcmc, res.abc, opt)
    
    # plot data + predictive densities
    simplemarkov.plot.data.and.true.pred(y, y.obs, y.pred.obs, y.pred.true, y.pred.mcmc, y.pred.abc, y.pred.abcf, opt)
    
    # another illustration of predictive densities
    simplemarkov.plot.example.preds(res.mcmc, res.abc, res.abcf, y.pred.true, y.obs, n, pp.fig, opt)
  }
  invisible()
}

################################################################################
simplemarkov.mcmc <- function(theta0, y, ns) {
  # A simple Gibbs sampler for the 3 parameters of the simple Markov model. 
   
  n <- length(y)
  cs <- c(theta0[1],rep(0,ns-1))
  phis <- c(theta0[2],rep(0,ns-1))
  sigma2s <- c(theta0[3],rep(0,ns-1))
  
  # precomputing:
  rns <- matrix(rnorm(2*ns,0,1),2,ns)
  sumy <- sum(y)
  sumynon <- sum(y[1:(n-1)])
  sumy2 <- sum(y^2)
  sumy2non <- sum(y[1:(n-1)]^2)
  sumcy <- sum(y[2:n]*y[1:(n-1)])
  
  # Run Gibbs sampler
  for (j in (2:ns)) {
    # c:
    cs[j] <- ((1-phis[j-1])*sumynon + y[n])/n + sqrt(sigma2s[j-1]/n)*rns[1,j]
    
    # phi:
    phis[j] <- (sumcy - cs[j]*sumynon)/sumy2non + sqrt(1/sumy2non)*rns[2,j]
    
    # sigma2:
    bj <- (sumy2 + phis[j]^2*sumy2non - 2*phis[j]*sumcy - 2*cs[j]*sumy + 
           2*cs[j]*phis[j]*sumynon + n*cs[j]^2)/2
    #if (bj < 0) {
    #  stop('Negative rate parameter.')
    #}
    sigma2s[j] <- 1/rgamma(n = 1, shape = n/2-1, rate = bj)
  }
  
  return(list(thetas=rbind(cs,phis,sigma2s)))
}

simplemarkov.pred.post.stats.mcmc <- function(samples, pp, y.obs, n, q) {
  # Computes mean and q-quantiles of the predictive distribution based on MCMC 
  # (or ABC) samples (currently no weights). 
  
  cs <- samples[1,]
  phis <- samples[2,]
  sigma2s <- samples[3,]
  n.pred <- length(pp)
  pred <- list(mean = rep(0,n.pred), u = rep(0,n.pred), l = rep(0,n.pred))
  for (p in pp) {
    mes <- cs*(1-phis^p)/(1-phis) + phis^p*y.obs[n]
    sds <- sqrt(pmax(0,sigma2s*(1-phis^(2*p))/(1-phis^2)))
    # mean of Gaussian mixture:
    pred$mean[p] <- mean(mes)
    # quantiles of Gaussian mixture:
    pred$u[p] <- gaussian.mixture.quantile(1-q/2,1/length(cs),mes,sds)
    pred$l[p] <- gaussian.mixture.quantile(q/2,1/length(cs),mes,sds)
  }
  return(pred)
}

simplemarkov.pred.post.stats.true <- function(theta.true, pp, y.obs, n, q) {
  # Computes mean, sd and q-quantiles of the predictive distribution based on the true parameter.
  
  pred <- list(theta.true=theta.true)
  pred$mean <- theta.true[1]*(1-theta.true[2]^pp)/(1-theta.true[2]) + theta.true[2]^pp*y.obs[n]
  pred$sd <- sqrt(pmax(0,theta.true[3]*(1-theta.true[2]^(2*pp))/(1-theta.true[2]^2)))
  pred$u <- pred$mean + qnorm(1-q/2)*pred$sd
  pred$l <- pred$mean + qnorm(q/2)*pred$sd
  return(pred)
}

simplemarkov.pred.post.dens.mcmc <- function(xp, samples, pp, y.obs, n) {
  # Computes full predictive distribution at points 'xp' based on MCMC 
  # (or ABC) samples (currently no weights)
  # Returns density values in #xp x #pp matrix.
  
  cs <- samples[1,]
  phis <- samples[2,]
  sigma2s <- samples[3,]
  gmp <- matrix(0,dim(xp)[1],length(pp))
  for (i in 1:length(pp)) {
    p <- pp[i]
    mes <- cs*(1-phis^p)/(1-phis) + phis^p*y.obs[n]
    sds <- sqrt(pmax(0,sigma2s*(1-phis^(2*p))/(1-phis^2)))
    for (j in (1:dim(xp)[1])) {
      gmp[j,i] <- mean(dnorm(xp[j,i], mean = mes, sd = sds))
    }
  }
  return(gmp)
}

################################################################################

simplemarkov.plot.params <- function(res.mcmc, res.abc, opt) {
  # Plots posterior densities obtained using MCMC and ABC.
  
  ins <- c('MCMC','ABC-P/F')
  cols <- c('black','blue')
  pns <- c('c','\\phi','\\sigma^2')
  
  library(latex2exp)
  d <- dim(res.mcmc$thetas)[1]
  kde <- function(samples) density(samples, adjust = 1.5)
  fn <- file.path(opt$save.loc, 'params_plot.pdf')
  pdf(file=fn, width = 2*d, height = 2)
  par(mfrow=c(1,d))
  par(mai=c(0.45,0.2,0.05,0.05), mgp=c(2,0.5,0))
  for (i in 1:d) {
    pmcmc <- kde(res.mcmc$thetas[i,])
    pabc <- kde(res.abc$thetas[i,])
    xla <- TeX(paste0('$',pns[i],'$'))
    rax <- range(pmcmc$x, pabc$x)
    ray <- range(pmcmc$y, pabc$y)
    plot(pabc,col=cols[2],ylab = '',xlab = xla, main = '',yaxt='n', xlim=rax, ylim=ray)
    lines(pmcmc,col=cols[1])
    if (i==1) {
      legend(x='topright', inset = 0.02, legend=ins, col=cols, lty=rep(1,2), bg = "white", cex=0.8)
    }
  }
  dev.off()
}

simplemarkov.plot.data.and.true.pred <- function(y, y.obs, y.pred.obs, y.pred.true, y.pred.mcmc, y.pred.abc, y.pred.abcf, opt) {
  # Plots observed data and predictions.
  # TODO: Better adjustment of plotting limits. 
  
  lw <- 1.5
  ins <- c('True param.','MCMC','ABC-P','ABC-F')
  cols <- c('orange','black','blue','magenta')
  col.data <- 'black'
  col.fdata <- 'orange'
  
  n <- length(y.obs)
  nn <- length(y)
  t.pred <- (n+1):nn
  t.all <- 1:nn
  
  for (i in 1:2) {
    fn <- file.path(opt$save.loc, paste0('ypred_plot_',i,'.pdf'))
    pdf(file=fn, width = 6, height = 6)
    par(mai=c(0.4,0.6,0.05,0.05), mgp=c(2,0.5,0))
    if (i==1) {
      # plot all observed data
      plot(1:n, y.obs, type='p', lwd=0.5, col = col.data, xlab = 't', pch=1, cex=0.6, ylab = 'y', 
           xlim = c(0,nn), ylim = c(min(y),1.1*max(y)))
    } else {
      # plot only part of observed data -> better visualization
      n.hist <- unique(pmax(1,(n-100):n))
      plot(n.hist, y.obs[n.hist], type='p', lwd=0.5, col = col.data, xlab = 't', pch=1, cex=0.6, 
           ylab = 'y', xlim = c(n.hist[1],nn), ylim = c(min(y[n.hist]),1.1*max(y[n.hist])))
    }
    lines(t.pred, y.pred.obs, type='p', lwd=0.5, col=col.fdata, pch=1, cex=0.6) # 'future' data
    
    # prediction based on true parameter
    simplemarkov.plot.single.pred(t.pred, y.pred.true, cols[1], lw)
    # prediction (MCMC)
    simplemarkov.plot.single.pred(t.pred, y.pred.mcmc, cols[2], lw)
    # prediction (ABC)
    simplemarkov.plot.single.pred(t.all, y.pred.abc, cols[3], lw)
    # prediction (ABC + exact pred.lik.)
    simplemarkov.plot.single.pred(t.pred, y.pred.abcf, cols[4], lw)
    
    legend(x='topleft', inset = 0.02, legend=ins, col=cols, lty=rep(1,4), lwd=lw, bg = "white", cex=0.8)
    dev.off()
  }
}

simplemarkov.plot.single.pred <- function(t, stats, col, lw, type='l') {
  # Plotting help function:
  lines(t, stats$mean, type = type, col = col, lty = 'solid', lwd=lw) # mean prediction
  lines(t, stats$u, type = type, col = col, lty = 'dashed', lwd=lw) # 95% CI
  lines(t, stats$l, type = type, col = col, lty = 'dashed', lwd=lw) # 95% CI
}

simplemarkov.plot.example.preds <- function(res.mcmc, res.abc, res.abcf, y.pred.true, y.obs, n, pp.fig, opt) {
  # Plots full posterior predictive densities in some example cases.
  
  ins <- c('True param.','MCMC','ABC-P','ABC-F')
  cols <- c('orange','black','blue','magenta')
  
  library(latex2exp)
  # MCMC: compute and plot the Gaussian mixture, ABC: plot the kde from simulated samples
  kde <- function(samples) density(samples, adjust = 1.5)
  # determine suitable plotting region
  nx <- 400 # how many grid points
  npp <- length(pp.fig)
  xp <- matrix(0,nx,npp)
  for (i in 1:npp) {
    xp[,i] <- seq(0.99*min(res.abc$y.sim[n+pp.fig[i],]),1.01*max(res.abc$y.sim[n+pp.fig[i],]),len=nx)
  }
  gmp.mcmc <- simplemarkov.pred.post.dens.mcmc(xp, res.mcmc$thetas, pp.fig, y.obs, n)
  gmp.abcf <- simplemarkov.pred.post.dens.mcmc(xp, res.abc$thetas, pp.fig, y.obs, n)
  
  fn <- file.path(opt$save.loc, 'ypred_plot_v2.pdf')
  pdf(file=fn, width = 2*npp, height = 2)
  par(mfrow=c(1,npp))
  par(mai=c(0.45,0.2,0.05,0.05), mgp=c(2,0.5,0))
  for (i in 1:npp) {
    # 'true' pred:
    xla <- TeX(paste0('$\\tilde{y}_{',n+pp.fig[i],'}$'))
    gp.true <- dnorm(xp[,i], mean = y.pred.true$mean[pp.fig[i]], sd = y.pred.true$sd[pp.fig[i]])
    ray <- c(0, 1.01*max(gp.true,gmp.mcmc[,i],gmp.abcf[,i]))
    plot(xp[,i],gp.true,col=cols[1], type='l' ,ylab = '',xlab = xla, main = '', yaxt='n', ylim=ray)
    
    # MCMC pred:
    lines(xp[,i],gmp.mcmc[,i],col=cols[2], type='l')
    
    # ABC pred:
    pabc <- kde(res.abc$y.sim[n+pp.fig[i],])
    lines(pabc,col=cols[3], type='l')
    
    # ABC pred with exact pred.lik:
    lines(xp[,i],gmp.abcf[,i],col=cols[4], type='l')
    
    if (i==1) {
      legend(x='topright', inset = 0.02, legend=ins, col=cols, lty=rep(1,4), bg = "white", cex=0.7)
    }
  }
  dev.off()
}

