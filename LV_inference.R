lv.inference <- function(opt, inf.task, plot.task) {
  # Runs the inference algorithms and plots the results for the L-V example 
  # using the options given it the corresponding setup-file.
  # TODO: Currently future population observations only predicted, add option 
  # to predict non-noisy population sizes.
  
  # Extract settings
  opt.abc <- opt$opt.abc
  opt.plt <- opt$plt
  for (i in 1:length(opt$sce)) {assign(names(opt$sce)[i],opt$sce[[i]])}
  
  library('matrixStats')
  dyn.load(file.path(opt$root,"LV_simul.so"))
  
  # set-up times
  t.sum <- t+t.pred+t2
  t.all <- seq(0,t.sum,by=dt) # all observation times
  nn <- length(t.all)
  n <- t/dt+1
  ind.obs <- 1:n
  ind.obs2 <- ifelsem(t2>0,(nn-t2/dt):nn,NULL)
  ind.pred <- setdiff(1:nn,union(ind.obs,ind.obs2))
  t.end <- t.all[n] # last observation time
  
  # generate 'observed' data from the LV model:
  set.seed(seed.data)
  obs <- lv.simul(theta.true, xy0, t.all, obs.noise.stdev)
  
  for (a in 1:length(opt.abc)) {
    if (!inf.task[a] || is.null(opt.abc[[a]])) {
      next
    }
    
    set.seed(seed.data)
    # 1/3 infer parameters using ABC-MCMC
    #####################################
    # set model
    # NOTE: parameters are in **log-domain**
    # NOTE2: simulated latent variables (future data, possible missing observations)
    # are included in obs$xyobs (and not in a separate 'lat'-list variable) and are
    # distinguished from the actual observed data only in sumstats-function!
    sim.model <- function(log.theta) {
      # Returns 2 x nn matrix of simulated populations. 
      # Full noiseless simulated populations as additional latents.
      obs <- lv.simul(exp(log.theta), xy0, t.all, obs.noise.stdev)
      return(list(y.sims=obs$xyobs, lat=obs$xy)) # NOTE: xy==NA in non-noisy case
    }
    
    # set the uniform prior in **log-domain**
    log.prior <- function(log.theta) logpmunif(log.theta, log.prior.lb, log.prior.ub)
    prior.sim <- function() rmunif(1, log.prior.lb, log.prior.ub)
    
    # select summary statistics and discrepancy:
    if (opt.abc[[a]]$summary=='misdata') {
      # "MISSING DATA" CASE
      sumstats <- function(xy) {
        if (is.list(xy)) {xy <- xy$y.sims}
        xy1 <- xy[xy.ind,ind.obs,drop=F] # DETERMINES OBSERVED DATA!
        xy2 <- xy[xy.ind,ind.obs2,drop=F]
        # Simplified summaries!!:
        # mean and stdev:
        s <- NULL
        s <- c(rowMeans(xy1), sqrt(rowVars(xy1)), rowMeans(xy2), sqrt(rowVars(xy2)))
        # summaries for *predicting missing data*:
        s <- c(s, xy1[,ncol(xy1)], xy2[,1])
        return(s)
      }
    } else if (any(opt.abc[[a]]$summary==c('orig','pred'))) {
      # INFERENCE/PREDICTION CASE
      sumstats <- function(xy) {
        if (is.list(xy)) {xy <- xy$y.sims}
        xy <- xy[xy.ind,1:n,drop=F] # DETERMINES OBSERVED DATA!
        # mean, stdev, lag1 and 2 autocorrelations:
        s <- NULL
        for (i in 1:dim(xy)[1]) {
          s <- c(s, lv.mean.sd.acf12(xy[i,]))
        }
        if (all(xy.ind)) {
          s <- c(s, cor(xy[1,],xy[2,])) # cross-correlation
        }
        # summaries for *prediction*:
        if (opt.abc[[a]]$summary=='pred') {
          s <- c(s, xy[,n])
        }
        return(s)
      }
    } else {
      stop('Incorrect summary statistics.')
    }
    
    # determine how many special summaries related to prediction
    sdim.not <- 0
    if (opt.abc[[a]]$summary=='misdata') {
      sdim.not <- 2*sum(xy.ind)
    } else if (opt.abc[[a]]$summary=='pred') {
      sdim.not <- sum(xy.ind)
    }
    
    sumstats.obs <- sumstats(obs$xyobs) # observed summary
    sdim <- length(sumstats.obs)
    invW <- abc.dist.invW(sdim, sim.model, sumstats, opt.abc[[a]]$d.n, opt.abc[[a]]$d.cov.method, 
                          opt.abc[[a]]$d.sim.method, th=log(theta.true), prior.sim, opt.plt$pr, sdim.not)
    
    discr <- function(obs) {
      d <- sumstats(obs) - sumstats.obs
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
    res.abc <- simple.abcmcmc(log(opt.abc[[a]]$theta0), sim.model, log.prior, discr, opt.abc[[a]])
    # transform back from log-domain:
    res.abc$thetas <- exp(res.abc$thetas)
    if (opt.plt$pr) {cat('ABC done!\n')}
    if (!dir.exists(opt$save.loc)) {dir.create(opt$save.loc)} # create folder
    save(res.abc, obs, opt, file = opt$fn.samples[a])
  }
  
  # Special simulations for computing the approx. true result in 'mis.data' case
  ##############################################################################
  if (inf.task[length(opt.abc)+1] && t2 > 0 && obs.noise.stdev == 0) {
    if (opt.plt$pr) {cat('Computing approx. true result...\n')}
    set.seed(seed.inf)
    res.atrue <- mis.data.approx.true(obs$xyobs[,n], theta.true, t.end, t.all[ind.pred], t.all[ind.obs2[1]], 
                                      obs, ind.obs2, obs.noise.stdev, mis.obs, opt.plt$q)
    save(res.atrue, obs, opt, file = mis.obs$fn.samples)
    if (opt.plt$pr) {cat('Computation done!\n')}
  }
  
  if (plot.task) {
    graphics.off()
    set.seed(seed.data)
    
    res.abca <- rep(list(NULL),length(opt.abc))
    for (a in 1:length(opt.abc)) {
      if (!is.null(opt.abc[[a]])) {
        load(file = opt$fn.samples[a])
        res.abca[[a]] <- res.abc
        
        basic.abc.mcmc.check(res.abca[[a]])
        
        # predictive posterior
        res.abca[[a]]$pred.stats.x <- pred.post.stats(res.abca[[a]]$y.sims[1,,], opt.plt$q, T, t.all)
        res.abca[[a]]$pred.stats.y <- pred.post.stats(res.abca[[a]]$y.sims[2,,], opt.plt$q, T, t.all)
      } 
    }
    
    # Predictive distribution based on true theta and true population sizes at time n.
    # This applies for the standard prediction case. An approximate solution to 
    # the 'missing data' case is also implemented. 
    res.true <- NA
    if (t2 == 0) {
      if (any(is.na(obs$xy))) {
        xy0p <- obs$xyobs[,n] # non-noisy case
      } else {
        xy0p <- obs$xy[,n] # noisy case
      }
      res.true <- lv.pred(xy0p, theta.true, t.end, t.all[ind.pred], obs.noise.stdev, opt.plt$q)
    } else if (t2 > 0 && obs.noise.stdev == 0) {
      load(file = mis.obs$fn.samples) # 'res.atrue' is pre-computed
      res.true <- res.atrue # allows to conveniently use existing plotting code
      cat('max distance:\n')
      print(res.true$d.eps)
      cat('number of accepted simulations:\n')
      print(res.true$n)
    }
    
    # 2/3 ABC with exact predictive simulation
    ##########################################
    # We use the ABC posterior with 'ordinary' summary -> No separate parameter fitting
    # This applies for the non-noisy standard prediction case when prey/predator both observed.
    res.abcf <- NA
    if (t2 == 0 && obs.noise.stdev == 0 && all(xy.ind)) {
      xy0p <- obs$xyobs[,n]
      res.abcf <- lv.pred(xy0p, res.abca[[1]]$thetas, t.end, t.all[ind.pred], obs.noise.stdev, opt.plt$q)
    }
    
    # 3/3 ABC with exact predictive simulation and simulated latent populations
    ###########################################################################
    # We use the ABC posterior with 'predictive' summary -> No separate parameter fitting
    # This is now implemented only for the case where only prey or predator is observed
    # and when the evaluations are noiseless.
    res.abclat <- NA
    if (t2 == 0 && obs.noise.stdev == 0 && !all(xy.ind)) {
      if (xy.ind[1]) {
        xy0p <- rbind(obs$xyobs[1,n], res.abca[[2]]$y.sims[2,n,])
      } else {
        xy0p <- rbind(res.abca[[2]]$y.sims[1,n,], obs$xyobs[2,n])
      }
      res.abclat <- lv.pred(xy0p, res.abca[[2]]$thetas, t.end, t.all[ind.pred], obs.noise.stdev, opt.plt$q)
    }
    
    ## PLOT RESULTS
    ###############
    # plot MCMC chains
    for (a in 1:length(opt.abc)) {
      if (!is.null(opt.abc[[a]])) {
        dev.new()
        simple.plot.mcmc.chain(res.abca[[a]]$thetas)
      }
    }
    
    # plot parameter posterior
    lv.plot.params(res.abca, theta.true, opt, log.th=F)
    lv.plot.params(res.abca, theta.true, opt, log.th=T)
    
    # plot prediction
    lv.plot.pred(obs, xy.ind, t.all, ind.obs, ind.obs2, ind.pred, res.abca, res.abcf, res.abclat, res.true, opt)
  }
  invisible()
}

################################################################################

lv.mean.sd.acf12 <- function(x) {
  # Computes the mean, standard deviation, lag1 and 2 autocorrelation estimates
  # for the summary statistics, given data vector 'x'.
  me <- mean(x)
  va <- var(x)
  y <- x-me
  n <- length(y)
  cs <- c(sum(y[2:n]*y[1:(n-1)]), sum(y[3:n]*y[1:(n-2)]))/(n*va)
  return(c(me, sqrt(va), cs))
  #return(c(acf(x,lag.max=2,plot=F)$acf[2:3])) # seems slower
}

mis.data.approx.true <- function(xy0p, theta.true, t.end, t.preds, t2.start, obs, 
                                 ind.obs2, obs.noise.stdev, mis.obs, q) {
  # Computes approximately the true prediction for the 'mis.data' case. It is 
  # assumed that the true parameter value is known and the prediction is matched
  # exactly at the beginning and approximately at the end prediction point. 
  
  n.pred <- length(t.preds)
  if (n.pred*mis.obs$n.samples > 10^8) {
    stop('Too much memory needed.') # TODO: all simulations currently saved to memory
  }
  x.preds <- matrix(NA,n.pred,mis.obs$n.samples)
  y.preds <- matrix(NA,n.pred,mis.obs$n.samples)
  ds <- rep(NA,mis.obs$n.samples)
  # first simulate all
  for (j in 1:mis.obs$n.samples) {
    preds <- lv.simul(theta.true, xy0p, c(t.end,t.preds,t2.start), obs.noise.stdev)
    ds[j] <- max(abs(preds$xyobs[,n.pred+2] - obs$xyobs[,ind.obs2[1]]))
    x.preds[,j] <- preds$xyobs[1,2:(n.pred+1)] # first and last timepoint neglected
    y.preds[,j] <- preds$xyobs[2,2:(n.pred+1)]
  }
  # select those simulations with the smallest distances (at the last timepoint)
  n.final <- min(mis.obs$n.final, mis.obs$n.samples)
  d.eps <- sort(ds)[n.final]
  inds <- which(ds<=d.eps) # can actually produce >n.final samples but this is ok
  x.preds <- x.preds[,inds]
  y.preds <- y.preds[,inds]
  res.atrue <- list(d.eps=d.eps, n=length(inds), x.preds=x.preds,y.preds=y.preds)
  # compute the quantiles for plotting already here:
  res.atrue$pred.stats.x <- pred.post.stats(x.preds, q, T, t.preds)
  res.atrue$pred.stats.y <- pred.post.stats(y.preds, q, T, t.preds)
  return(res.atrue)
}

lv.pred <- function(xyns, thetas, t.end, t.preds, obs.noise.stdev, q) {
  # Computes the predictive density of the population sizes at some future time 
  # points 't.preds' given
  # 1) the population size (observed value or some samples of it) 'xyns' at the
  # last observation time 't.end' and
  # 2) the theta parameter (true value or ABC samples of it) 'thetas'.
  
  if (length(t.end)!=1 || t.end>=t.preds[1]) {
    stop('Incorrect initial time.') # check just in case
  }
  if (is.vector(thetas) && is.vector(xyns)) {
    n.sim <- 10000
  } else {
    n.sim <- c(dim(thetas)[2],dim(xyns)[2])[1] # use all provided samples
  }
  n.pred <- length(t.preds)
  x.preds <- matrix(NA,n.pred,n.sim); y.preds <- matrix(NA,n.pred,n.sim)
  th <- thetas
  xy0p <- xyns
  for (j in 1:n.sim) {
    if (!is.vector(thetas)) {
      th <- thetas[,j]
    } 
    if (!is.vector(xyns)) {
      xy0p <- xyns[,j]
    }
    # Note: the last observation time 't.end' need to be included to the future
    # prediction times and finally neglected from the output
    preds <- lv.simul(th, xy0p, c(t.end,t.preds), obs.noise.stdev)
    x.preds[,j] <- preds$xyobs[1,2:(n.pred+1)]
    y.preds[,j] <- preds$xyobs[2,2:(n.pred+1)]
  }
  pred <- list()
  pred$pred.stats.x <- pred.post.stats(x.preds, q, T, t.preds)
  pred$pred.stats.y <- pred.post.stats(y.preds, q, T, t.preds)
  return(pred)
}

################################################################################

lv.plot.params <- function(res.abca, theta.true, opt, log.th=F) {
  # Plots posterior densities obtained using ABC. Plots either the parameters 
  # or log parameters.
  
  library(latex2exp)
  ins <- c('True param.',TeX(c('ABC-P, $s^{(0)}$','ABC-P, $s^{(1)}$')))
  cols <- c('black','red','blue')
  
  d <- length(theta.true)
  kde <- function(samples) density(samples, adjust = 1.2)
  fn <- file.path(opt$save.loc, paste0(ifelse(log.th,'params_plot_log','params_plot'),opt$fn.ext,'.pdf'))
  pdf(file=fn, width = 2*d, height = 2)
  par(mfrow=c(1,d))
  par(mai=c(0.45,0.2,0.05,0.05), mgp=c(2,0.5,0))
  if (log.th) {
    f <- function(x) log(x)
  } else {
    f <- function(x) x
  }
  m <- c(T,!sapply(res.abca,is.null)) # which methods were computed
  pabc <- vector('list',length(res.abca))
  for (i in 1:d) {
    rax <- NULL; ray <- NULL
    for (j in 1:length(res.abca)) {
      if (m[j+1]) {
        pabc[[j]] <- kde(f(res.abca[[j]]$thetas[i,]))
        rax <- range(rax,pabc[[j]]$x)
        ray <- range(ray,pabc[[j]]$y)
      }
    }
    xla <- TeX(paste0('$',ifelse(log.th,'\\log\\,',''),'\\theta_',i,'$'))
    # true value as horizontal line:
    plot(rep(f(theta.true[i]),2),ray*c(1,1.1),type='l',col=cols[1],ylab = '',xlab = xla, main = '', xlim=rax, ylim=ray, yaxt='n')
    for (j in 1:length(res.abca)) {
      if (m[j+1]) {
        lines(pabc[[j]],col=cols[j+1])
      }
    }
    if (i==1) {
      legend(x='topleft', inset = c(0.02,0.02), legend=ins[m], col=cols[m], lty=rep(1,sum(m)), bg = "white", cex=0.6)
    }
  }
  dev.off()
}

lv.plot.pred <- function(obs, xy.ind, t.all, ind.obs, ind.obs2, ind.pred, res.abca, res.abcf, res.abclat, res.true, opt) {
  # Plots prediction as a function of time. 
  
  fn <- file.path(opt$save.loc, paste0('pred_plot_1',opt$fn.ext,'.pdf'))
  pdf(file=fn, width = 6, height = 5)
  par(mfrow=c(2,1))
  # plot prey:
  par(mai=c(0.5,0.6,0.1,0.05), mgp=c(1.5,0.5,0))
  lv.plot.popul(1, obs, xy.ind, t.all, ind.obs,ind.obs2,ind.pred, res.abca, res.abcf, res.abclat, res.true)
  # plot predator:
  par(mai=c(0.5,0.6,0.05,0.05), mgp=c(1.5,0.5,0))
  lv.plot.popul(2, obs, xy.ind, t.all, ind.obs,ind.obs2,ind.pred, res.abca, res.abcf, res.abclat, res.true)
  dev.off()
}

lv.plot.popul <- function(id, obs, xy.ind, t.all, ind.obs, ind.obs2, ind.pred, res.abca, res.abcf, res.abclat, res.true) {
  # Plotting help function that plots either the prey or predator population: 
  
  library(latex2exp)
  ins <- c('True param.',TeX(c('ABC-P, $s^{(0)}$','ABC-P, $s^{(1)}$','ABC-L, $s^{(1)}$','ABC-F, $s^{(0)}$')))
  cols <- c('black','red','blue','orange','orange')
  col.data <- 'black'
  col.fdata <- 'orange'
  lw <- 1.5
  ds <- 0.6 # data point size
  s <- c('pred.stats.x','pred.stats.y')
  ms <- c('res.true','res.abca[[1]]','res.abca[[2]]','res.abclat','res.abcf')
  
  # 1/2 plot observations:
  ml <- methods.and.limits(id, obs, t.all, res.abca, res.abcf, res.abclat, res.true, s, ms)
  if (id == 1) {
    if (xy.ind[1]) {
      plot(t.all[ind.obs], obs$xyobs[id,ind.obs], type = 'p', pch=16, cex=ds, col = col.data, xlab = '', #xaxt='n'
           ylab = 'Prey population', xlim = ml$rax, ylim = ml$ray)
    } else {
      plot(NULL, NULL, xlab = '', ylab = 'Prey population', xlim = ml$rax, ylim = ml$ray)
    }
  } else {
    if (xy.ind[2]) {
      plot(t.all[ind.obs], obs$xyobs[id,ind.obs], type = 'p', pch=16, cex=ds, col = col.data, xlab = 't', 
           ylab = 'Predator population', xlim = ml$rax, ylim = ml$ray)
    } else {
      plot(NULL, NULL, xlab = 't', ylab = 'Predator population', xlim = ml$rax, ylim = ml$ray)
    }
  }
  if (!is.null(ind.obs2) && xy.ind[id]) {
    # second set of observed data
    lines(t.all[ind.obs2], obs$xyobs[id,ind.obs2], type = 'p', pch=16, cex=ds, col = col.data)
  }
  if (1) {
    # gray line(s) where predictions start(/end)
    lines(rep(t.all[ind.pred[1]],2), ml$ray+100*c(-1,1), type = 'l', col = 'gray')
    if (!is.null(ind.obs2)) {
      lines(rep(t.all[ind.pred[length(ind.pred)]],2), ml$ray+100*c(-1,1), type = 'l', col = 'gray')
    }
  }
  ###lines(t.all[ind.pred], obs$xyobs[id,ind.pred], type = 'p', pch=16, cex=0.6, col = col.fdata) # 'future' data
  #lines(t.all,rep(0,length(t.all)),'o')
  
  # 2/2 plot predictions:
  for (i in 1:length(ms)) {
    if (ml$m[i]) {
      ti <- t.all[ind.pred]
      if (any(i==2:3)) {ti <- t.all} # all points prediction in ABC-P case
      lv.plot.single.pred(ti, eval(parse(text=paste0(ms[i],'$',s[id]))), col=cols[i], lw=lw)
    }
  }
  if (id==1) {
    legend(x='topleft', inset = c(0.02,0.02), legend=ins[ml$m], col=cols[ml$m], lty=rep(1,sum(ml$m)), 
           ncol = min(2,sum(ml$m)), bg = "white", cex=0.65)
  }
}

methods.and.limits <- function(id, obs, t.all, res.abca, res.abcf, res.abclat, res.true, s, st) {
  # Help function that checks which methods were computed and determines suitable plotting limits.
  
  m <- rep(F,length(st)) # which methods were computed
  maxy <- max(obs$xyobs[id,]) # max value for y-axis
  for (i in 1:length(st)) {
    if (!is.na(eval(parse(text=st[i]))) && !is.null(eval(parse(text=st[i])))) {
      m[i] <- T
      maxy <- max(maxy, eval(parse(text=paste0(st[i],'$',s[id],'$u1'))))
    }
  }
  maxy <- min(maxy,1.5*max(obs$xyobs[id,])) # needed for a rare exp-growing case
  pd <- 1
  return(list(m=m, rax=c(pd,max(t.all)-pd), ray=c(0,maxy)))
}

lv.plot.single.pred <- function(t, stats, col, lw, type='l') {
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

