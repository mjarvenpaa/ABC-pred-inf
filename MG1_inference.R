mg1.inference <- function(opt, inf.task, plot.task) {
  # Runs the inference algorithms and plots the results for the M/G/1 example 
  # using the options given it the corresponding setup-file.
  
  # Extract settings
  opt.mcmc <- opt$opt.mcmc
  opt.abc <- opt$opt.abc
  opt.plt <- opt$plt
  for (i in 1:length(opt$sce)) {assign(names(opt$sce)[i],opt$sce[[i]])}
  
  library('matrixStats')
  dyn.load(file.path(opt$root,"MG1_simul.so"))
  
  nn <- n + n.pred
  # generate 'observed' data from the M/G/1 model:
  set.seed(seed.data)  
  obs <- mg1.simul(theta.true, nn)
  y.obs <- obs$y[1:n]; x.obs <- obs$x[1:n]; v.obs <- obs$v[1:n]
  th1.ub <- min(y.obs)
    
  if (inf.task[1]) {
    # 1/3 infer parameters using MCMC
    #################################
    # first determine initialization for MCMC
    set.seed(seed.inf)
    if (0) {
      # crude initialization
      params0 <- init.theta.v(y.obs, x.obs, theta.true)
    } else {
      # initialize using true params and corresponding simulated v values
      params0 <- list(theta = theta.true, v = v.obs)
    }
    opt.mcmc$params0 <- params0
    # check that the initialization - no matter how determined - is feasible.
    if (!is.all.constraint.ok(params0$theta, y.obs, x.obs, params0$v)) {
      stop('Initialization not OK.')
    }
    
    if (opt.plt$pr) {cat('Running MCMC...\n')}
    res.mcmc <- mg1.mcmc(params0, y.obs, opt.mcmc, n.pred)
    if (opt.plt$pr) {cat('MCMC done!\n')}
    if (!dir.exists(opt$save.loc)) {dir.create(opt$save.loc)} # create folder
    save(res.mcmc, obs, opt, file = opt$fn.samples.mcmc)
  }
  
  # Run ABC analyses
  for (a in 1:length(opt.abc)) {
    if (!inf.task[1+a]) {
      next
    }
    set.seed(seed.data)
    # 2/3 infer parameters using ABC-MCMC (produces also ABC-P)
    ###########################################################
    if (is.null(opt.abc[[a]]$C) || is.na(opt.abc[[a]]$C)) {
      # option to use the cov matrix from standard MCMC run for ABC-MCMC:
      if (!inf.task[1]) {
        load(file = opt$fn.samples.mcmc) # load MCMC samples if they were not computed
      }
      opt.abc[[a]]$C <- cov(t(res.mcmc$thetas))
    }
    
    # set model
    sim.model <- function(theta) {
      obs <- mg1.simul(theta, nn)
      return(list(y.sims = obs$y, lat = obs$v)) # NOTE: v as latent variables
    }
    
    # set prior
    #th1.ub <- Inf # additional constraint not otherwise explicitly involved in ABC analysis
    log.prior <- function(theta) ifelse(is.theta.prior.constraint.ok(theta)&&theta[1]<th1.ub,0,-Inf)
    #prior.sim <- function() rmunif(1, log.prior.lb, log.prior.ub) # TODO:
    
    # set summary statistics and discrepancy:
    if (opt.abc[[a]]$summary=='orig') {
      sumstats <- function(y) {
        if (is.list(y)) {y <- y$y.sims}
        y <- y[1:n]
        return(quantile(y, probs=c(0, 0.25, 0.5, 0.75, 1), names=F)) # "original summary"
        #return(quantile(y, probs=c(0, 0.2, 0.35, 0.5, 0.65, 0.8, 1), names=F)) # more quantiles
      }
    } else if (opt.abc[[a]]$summary=='pred') {
      # First attempt
      sumstats <- function(y) {
        if (is.list(y)) {y <- y$y.sims}
        y <- y[1:n]
        x <- cumsum(y)
        return(c(quantile(y, probs=c(0, 0.25, 0.5, 0.75), names=F), x[(n-1):n])) # "fairly good" summary for predicting v
        #return(c(quantile(y, probs=c(0, 0.25, 0.5, 0.75, 1), names=F), y[n]))
        #return(c(quantile(y, probs=c(0, 0.25, 0.5, 0.75, 1), names=F), y[(n-1):n]))
      }
    } else if (opt.abc[[a]]$summary=='pred2') {
      # A bit ad-hoc but seems more suitable for predictive ABC than 'pred':
      sumstats <- function(y) {
        if (is.list(y)) {y <- y$y.sims}
        y <- y[1:n]
        return(c(quantile(y, probs=c(0, 0.25, 0.5, 0.75, 1), names=F), range(y[(n-10):n])))
      }
    } else if (opt.abc[[a]]$summary=='predf') {
      # Quite rudimentary summary stats but suitable for predictive ABC:
      #qs.obs <- quantile(y.obs, probs = c(0.5, 0.75, 0.9), names=F) # 0.5 case actually not that useful...
      qs.obs <- quantile(y.obs, probs = c(0.7, 0.8, 0.9), names=F)
      sumstats <- function(y) {
        if (is.list(y)) {y <- y$y.sims}
        y <- y[1:n]
        sq1 <- max(1,which(y>=qs.obs[1]))
        sq2 <- max(1,which(y>=qs.obs[2]))
        sq3 <- max(1,which(y>=qs.obs[3]))
        return(c(quantile(y, probs=c(0, 0.25, 0.5, 0.75, 1), names=F), sq1, sq2, sq3))
      }
    } else {
      stop('Incorrect summary statistics.')
    }
    
    # determine how many special summaries related to prediction
    sdim.not <- 0
    if (opt.abc[[a]]$summary=='pred') {
      sdim.not <- 2 # NOTE: Need to be adjusted if summary stats are changed!
    } else if (opt.abc[[a]]$summary=='pred2') {
      sdim.not <- 2 # NOTE: Need to be adjusted if summary stats are changed!
    } else if (opt.abc[[a]]$summary=='predf') {
      sdim.not <- 3 # NOTE: Need to be adjusted if summary stats are changed!
    }
    
    sumstats.obs <- sumstats(y.obs) # observed summary
    print(sumstats.obs)
    sdim <- length(sumstats.obs)
    invW <- abc.dist.invW(sdim, sim.model, sumstats, opt.abc[[a]]$d.n, opt.abc[[a]]$d.cov.method, 
                          opt.abc[[a]]$d.sim.method, th=theta.true, prior.sim, opt.plt$pr, sdim.not)
    
    discr <- function(obs) {
      d <- sumstats(obs$y.sims) - sumstats.obs
      if (sdim.not==0) {
        return(sqrt(d%*%invW%*%d))
      }
      # special weighting for summaries related to prediction:
      d1 <- d[1:(sdim-sdim.not)]
      dpred <- d[(sdim-sdim.not+1):sdim]
      return(c(sqrt(d1%*%invW%*%d1), max(abs(dpred)))) # 2D output
    }
    
    # Run ABC-MCMC!
    if (opt.plt$pr) {cat('Running ABC...\n')}
    set.seed(seed.inf)
    res.abc <- simple.abcmcmc(opt.abc[[a]]$theta0, sim.model, log.prior, discr, opt.abc[[a]])
    if (opt.plt$pr) {cat('ABC done!\n')}
    save(res.abc, obs, opt, file = opt$fn.samples[a])
  }
  
  
  # Compute the rest of the ABC results (fast) and then plot/analyze them
  if (plot.task) {
    graphics.off()
    set.seed(seed.data)
    
    # load samples and check them
    load(file = opt$fn.samples.mcmc)
    basic.abc.mcmc.check(res.mcmc)
    res.abca <- rep(list(NULL),length(opt.abc))
    for (a in 1:length(opt.abc)) {
      load(file = opt$fn.samples[a])
      res.abca[[a]] <- res.abc
      basic.abc.mcmc.check(res.abca[[a]])
    }
    
    # MCMC-based predictive posterior for y
    res.mcmc$pred.y.stats <- pred.post.stats(res.mcmc$ypreds, opt.plt$q, T, 1:n.pred)
    
    # MCMC-based predictive posterior for waiting times
    res.mcmc$pred.wts <- res.mcmc$xpreds - res.mcmc$vpreds
    res.mcmc$pred.wts.stats <- pred.post.stats(res.mcmc$pred.wts, opt.plt$q, T, 1:n.pred)
    
    # predictive posterior based on true theta and true v
    res.true <- mg1.pred(y.obs, x.obs, v.obs, theta.true, n.pred, opt.plt$q, true.pred = T)
    
    for (a in 1:length(opt.abc)) {
      # ABC-based predictive posterior for y (ABC-P)
      res.abca[[a]]$pred.y.stats <- pred.post.stats(res.abca[[a]]$y.sims, opt.plt$q, T, 1:nn)
      
      # ABC-based predictive posterior for waiting times (ABC-P)
      #res.abc$pred.wts <- colCumsums(res.abca[[a]]$y.sims)[(n+1):nn,] - res.abca[[a]]$lat[(n+1):nn,] # 'matrixStats'
      res.abca[[a]]$pred.wts <- colCumsums(res.abca[[a]]$y.sims) - res.abca[[a]]$lat # 'matrixStats'
      res.abca[[a]]$pred.wts.stats <- pred.post.stats(res.abca[[a]]$pred.wts, opt.plt$q, T, 1:nn)
    }
    
    # 3/3 ABC with exact predictive lik. simulation given latent variable v[n] (ABC-L)
    ##################################################################################
    # We use the same ABC posterior as above -> No separate parameter fitting.
    # NOTE: We consider the second set of summary stats only.
    res.abclat <- mg1.pred(y.obs, x.obs, res.abca[[2]]$lat.sims, res.abca[[2]]$thetas, n.pred, opt.plt$q)
    
    # plot MCMC chains
    simple.plot.mcmc.chain(rbind(res.mcmc$thetas,res.mcmc$vs[c(1:3,(n-2):n),]))
    for (a in 1:length(opt.abc)) {
      dev.new()
      simple.plot.mcmc.chain(res.abca[[a]]$thetas)
    }
    
    # plot parameter posterior
    mg1.plot.params(n, res.mcmc, res.abca, opt, plot.eta=T, th1.ub=th1.ub)
    mg1.plot.params(n, res.mcmc, res.abca, opt, th1.ub=th1.ub)
    cat('min(y):\n')
    print(min(y.obs))
    
    # plot predictions
    mg1.plot.pred.dens(n, res.mcmc, res.abca, res.abclat, res.true, opt) 
    mg1.plot.pred.dens(n, res.mcmc, res.abca, res.abclat, res.true, opt, plot.wts=T)
    
    mg1.plot.pred(obs, n, n.pred, res.mcmc, res.abca, res.abclat, res.true, opt)
    mg1.plot.pred(obs, n, n.pred, res.mcmc, res.abca, res.abclat, res.true, opt, plot.wts=T)
    
    mg1.plot.just.data(obs, n, opt)
  }
  invisible()
}

################################################################################

mg1.pred <- function(y.obs, x.obs, vs, thetas, n.pred, q, true.pred=F) {
  # Computes predictive densities using either ABC/MCMC samples or true values 
  # of v and theta by simulating from the exact pred.lik given y and the latent
  # variable v[n].
  
  if (true.pred) {
    n.sim <- 10000
  } else {
    n.sim <- dim(thetas)[2] # use all provided v and theta samples
  }
  y.preds <- matrix(NA,n.pred,n.sim); x.preds <- matrix(NA,n.pred,n.sim)
  v.preds <- matrix(NA,n.pred,n.sim)
  n <- length(y.obs)
  for (j in 1:n.sim) {
    if (true.pred) {
      # vs and thetas assumed to contain only single value, the true parameter
      preds <- mg1.simul.future(vs[n], x.obs[n], thetas, n.pred)
    } else {
      preds <- mg1.simul.future(vs[n,j], x.obs[n], thetas[,j], n.pred)
    }
    y.preds[,j] <- preds$ypreds; x.preds[,j] <- preds$xpreds
    v.preds[,j] <- preds$vpreds
  }
  pred <- list(ypreds = y.preds)
  pred$pred.y.stats <- pred.post.stats(y.preds, q, T, 1:n.pred)
  pred$pred.wts <- x.preds - v.preds
  pred$pred.wts.stats <- pred.post.stats(pred$pred.wts, q, T, 1:n.pred)
  return(pred)
}

################################################################################
# Plotting etc. functions:
# NOTE: All plotting functions below currently assume two set of summary stats 
# were considered i.e. length of res.abca is 2 and both elements contain results. 

mg1.plot.params <- function(n, res.mcmc, res.abca, opt, plot.eta=F, th1.ub=Inf) {
  # Plots posterior densities obtained using MCMC and ABC.
  
  library(latex2exp)
  ins <- c('MCMC',TeX(c('ABC-P, $s^{(0)}$','ABC-P, $s^{(1)}$')))
  cols <- c('magenta','red','blue')
  v.name <- 'v' # interarrival time
  
  # 1/2: theta (or eta)
  d <- dim(res.mcmc$thetas)[1]
  fn <- file.path(opt$save.loc, paste0('params_plot',ifelse(plot.eta,'_eta',''),opt$fn.ext,'.pdf'))
  pdf(file=fn, width = 2*d, height = 2)
  par(mfrow=c(1,d))
  par(mai=c(0.35,0.15,0.01,0.07), mgp=c(1.8,0.5,0))
  if (plot.eta) {
    f <- function(th,i) if(i==1){th[1,]} else if(i==2){th[2,]-th[1,]} else{log(th[3,])}
  } else {
    f <- function(th,i) th[i,]
  }
  for (i in 1:d) {
    if (i==1) {
      kde <- function(samples) density.bound(samples, ub=th1.ub, adjust=1.3) # take account ub for theta_1
    } else {
      kde <- function(samples) density(samples, adjust=1.3)
    }
    pmcmc <- kde(f(res.mcmc$thetas,i))
    pabc1 <- kde(f(res.abca[[1]]$thetas,i))
    pabc2 <- kde(f(res.abca[[2]]$thetas,i))
    xla <- TeX(paste0('$',ifelse(plot.eta,'\\eta','\\theta'),'_',i,'$'))
    rax <- range(pmcmc$x, pabc1$x)
    ray <- range(pmcmc$y, pabc1$y)
    plot(pabc1,col=cols[2],ylab = '',xlab = xla, main = '',yaxt='n',xlim=rax, ylim=ray)
    lines(pabc2,col=cols[3])
    lines(pmcmc,col=cols[1])
    if (i==1) {
      legend(x='topleft', inset = 0.02, legend=ins, col=cols, lty=rep(1,3), bg = "white", cex=0.8)
    }
  }
  dev.off()
  if (plot.eta) {
    return() # latents v not plotted in this case
  }
  
  # 2/2: v[n-2],v[n-1],v[n]
  d <- 3
  fn <- file.path(opt$save.loc, paste0('params_plot_v',opt$fn.ext,'.pdf'))
  pdf(file=fn, width = 2*d, height = 2)
  par(mfrow=c(1,d))
  par(mai=c(0.35,0.15,0.01,0.07), mgp=c(1.8,0.5,0))
  kde <- function(samples) density(samples, adjust=1.2)
  for (i in 1:d) {
    ind <- n-d+i
    pmcmc <- kde(res.mcmc$vs[ind,])
    pabc1 <- kde(res.abca[[1]]$lat.sims[ind,])
    pabc2 <- kde(res.abca[[2]]$lat.sims[ind,])
    xla <- TeX(paste0('$',v.name,'_{',ind,'}$'))
    #rax <- range(pmcmc$x, pabc1$x)
    #ray <- range(pmcmc$y, pabc1$y)
    rax <- range(pabc1$x, pabc2$x)
    ray <- range(pabc1$y, pabc2$y)
    plot(pabc1,col=cols[2],ylab = '',xlab = xla, main = '',yaxt='n',xlim=rax, ylim=ray)
    lines(pabc2,col=cols[3])
    lines(pmcmc,col=cols[1])
    if (i==1) {
      legend(x='topleft', inset = 0.02, legend=ins, col=cols, lty=rep(1,2), bg = "white", cex=0.8)
    }
  }
  dev.off()
}

mg1.plot.pred.dens <- function(n, res.mcmc, res.abca, res.abclat, res.true, opt, plot.wts=F) {
  # Plots prediction density of interarrival times y or the corresponding waiting times, 
  # at some individual future times. 
  
  library(latex2exp)
  ins <- c('True param.','MCMC',TeX(c('ABC-P, $s^{(0)}$','ABC-P, $s^{(1)}$','ABC-L, $s^{(1)}$')))
  cols <- c('black','magenta','red','blue','orange')
  #y.name <- '\\tilde{y}' # new interdeparture time
  y.name <- 'y' # new interdeparture time
  wt.name <- '\\tilde{\\omega}' # new waiting time
  #wt.name <- '\\omega' # new waiting time. (This was originally used so that tilde was missing.)
  
  pp.fig <- opt$plt$pp.fig
  npp <- length(pp.fig)
  fn <- file.path(opt$save.loc, paste0(ifelse(plot.wts,'pred_plot_1_wts','pred_plot_1'),opt$fn.ext,'.pdf'))
  pdf(file=fn, width = npp*2, height = 2)
  par(mfrow=c(1,npp))
  par(mai=c(0.35,0.15,0.01,0.07), mgp=c(1.8,0.5,0))
  if (plot.wts) {
    kde <- function(samples) density.bound(samples, lb = 0, adjust = 1)
    get.samples <- function(res,i,o=0) res$pred.wts[o+pp.fig[i],]
  } else {
    kde <- function(samples) density.bound(samples, lb = 0, adjust = 1)
    get.samples <- function(res,i,o=0) ifelsem(o>0,res$y.sims[o+pp.fig[i],],res$ypreds[pp.fig[i],])
  }
  for (i in 1:npp) {
    pmcmc <- kde(get.samples(res.mcmc,i))
    pabc1 <- kde(get.samples(res.abca[[1]],i,n))
    pabc2 <- kde(get.samples(res.abca[[2]],i,n))
    pabclat <- kde(get.samples(res.abclat,i))
    ptrue <- kde(get.samples(res.true,i))
    maxx <- sort(c(max(pmcmc$x), max(pabc1$x), max(pabc2$x), max(pabclat$x)),T)[2] # useful when one case has long tail
    if (plot.wts && opt$scenario==201 && npp==3) {
      maxx <- min(maxx, c(100,150,200)[i]) # additional ad-hoc adjustment to plotting limits
    }
    rax <- c(min(pmcmc$x, pabc1$x, pabc2$x, pabclat$x), maxx)
    ray <- range(pmcmc$y, pabc1$y, pabc2$y, pabclat$y)
    xla <- TeX(paste0('$',ifelse(plot.wts,wt.name,y.name),'_{',n + pp.fig[i],'}$'))
    plot(pabc1,col=cols[3],ylab = '',xlab = xla,main = '',yaxt='n', xlim=rax, ylim=ray)
    lines(pabc2,col=cols[4])
    lines(ptrue,col=cols[1])
    lines(pmcmc,col=cols[2])
    lines(pabclat,col=cols[5])
    if (i==1) {
      leloc <- 'topright'
      # ad hoc fix for legend placement:
      if (plot.wts && opt$scenario==101) {leloc <- 'topleft'}
      legend(x=leloc, inset = 0.02, legend=ins, col=cols, lty=rep(1,5), bg = "white", cex=0.8)
    }
  }
  dev.off()
}

mg1.plot.pred <- function(obs, n, n.pred, res.mcmc, res.abca, res.abclat, res.true, opt, plot.wts=F) {
  # Plots predictions as a function of time. 
  
  library(latex2exp)
  ins <- c('True param.','MCMC',TeX(c('ABC-P, $s^{(0)}$','ABC-P, $s^{(1)}$','ABC-L, $s^{(1)}$')))
  cols <- c('black','magenta','red','blue','orange')
  col.data <- 'black'
  col.fdata <- 'orange'
  #y.name <- 'y_{i}' # interdeparture time
  y.name <- 'y'
  #wt.name <- '\\omega_{i}' # waiting time
  wt.name <- 'Waiting time'
  x.name <- 'Customer'
  lw <- 1.5
  
  nn <- n + n.pred
  t.obs <- 1:n
  t.pred <- (n+1):nn
  t.all <- 1:nn
  fn <- file.path(opt$save.loc, paste0(ifelse(plot.wts,'pred_plot_2_wts','pred_plot_2'),opt$fn.ext,'.pdf'))
  pdf(file=fn, width = 3.4, height = 3)
  par(mai=c(0.6,0.6,0.05,0.05), mgp=c(1.8,0.5,0))
  
  rax <- c(1,nn)
  if (plot.wts) {
    # waiting times:
    wts <- obs$x - obs$v
    ray <- c(1,1.05)*range(wts, res.mcmc$pred.wts.stats$l1, res.mcmc$pred.wts.stats$u1)
    plot(t.obs, wts[t.obs], type = 'o', cex=0.8, xlab = x.name, ylab=wt.name, col=col.data, xlim=rax, ylim=ray) # observed
    if (1) {
      # gray line where predictions start
      lines(rep(t.pred[1],2), ray+100*c(-1,1), type = 'l', col = 'gray')
    }
    ###lines(t.pred, wts[t.pred], type = 'o', cex=0.8, col = col.fdata) # 'future' data
    
    # prediction (true)
    mg1.plot.single.pred(t.pred, res.true$pred.wts.stats, col=cols[1], lw=lw)
    # prediction (MCMC)
    mg1.plot.single.pred(t.pred, res.mcmc$pred.wts.stats, col=cols[2], lw=lw)
    # prediction (ABC)
    mg1.plot.single.pred(t.all, res.abca[[1]]$pred.wts.stats, col=cols[3], lw=lw)
    mg1.plot.single.pred(t.all, res.abca[[2]]$pred.wts.stats, col=cols[4], lw=lw)
    mg1.plot.single.pred(t.pred, res.abclat$pred.wts.stats, col=cols[5], lw=lw)
  } else {
    # interarrival times y:
    ray <- range(obs$y, res.mcmc$pred.y.stats$l1, res.mcmc$pred.y.stats$u1)
    plot(t.obs, obs$y[t.obs], type = 'o', cex=0.8, xlab = x.name, ylab=y.name, col=col.data, xlim=rax, ylim=ray) # observed
    if (1) {
      # gray line where predictions start
      lines(rep(t.pred[1],2), ray+100*c(-1,1), type = 'l', col = 'gray')
    }
    lines(t.pred, obs$y[t.pred], type = 'o', cex=0.8, col = col.fdata) # 'future' data
    
    # prediction (true)
    mg1.plot.single.pred(t.pred, res.true$pred.y.stats, col=cols[1], lw=lw)
    # prediction (MCMC)
    mg1.plot.single.pred(t.pred, res.mcmc$pred.y.stats, col=cols[2], lw=lw)
    # prediction (ABC)
    mg1.plot.single.pred(t.all, res.abca[[1]]$pred.y.stats, col=cols[3], lw=lw)
    mg1.plot.single.pred(t.all, res.abca[[2]]$pred.y.stats, col=cols[4], lw=lw)
    mg1.plot.single.pred(t.pred, res.abclat$pred.y.stats, col=cols[5], lw=lw)
  }
  legend(x='topleft', inset = 0.02, legend=ins, col=cols, lty=rep(1,5), lwd=lw, bg = "white", cex=0.8)
  dev.off()
}

mg1.plot.single.pred <- function(t, stats, col, lw, type='l') {
  # Plotting help function:
  #lines(t, stats$mean, type = type, col = col, lty = 'dotted', lwd=lw)
  lines(t, stats$med, type = type, col = col, lty = 'solid', lwd=lw)
  lines(t, stats$u1, type = type, col = col, lty = 'dashed', lwd=lw)
  lines(t, stats$l1, type = type, col = col, lty = 'dashed', lwd=lw)
  if (!is.null(stats[["u2"]])) {
    lines(t, stats$u2, type = type, col = col, lty = 'dotdash', lwd=lw)
    lines(t, stats$l2, type = type, col = col, lty = 'dotdash', lwd=lw)
  }
}

mg1.plot.just.data <- function(obs, n, opt) {
  # Plots data nicely (and nothing else).
  
  #y.name <- 'y'
  y.name <- 'Interdeparture time'
  x.name <- 'Customer'
  cols <- 'black'
  ds <- 0.7
  
  fn <- file.path(opt$save.loc, paste0('data_plot',opt$fn.ext,'.pdf'))
  pdf(file=fn, width = 4.2, height = 3)
  par(mai=c(0.6,0.6,0.05,0.07), mgp=c(1.8,0.5,0))
  
  plot(1:n, obs$y[1:n], type = 'p', pch=16, cex=ds, xlab = x.name, ylab=y.name, col=cols, xlim=c(1,n)) # observed data
  lines(1:n, obs$y[1:n], type = 'l', col=cols) # adds lines between points
  dev.off()
}

