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
  
  # Run ABC analyses
  for (a in 1:length(opt.abc)) {
    if (!inf.task[a] || is.null(opt.abc[[a]])) {
      next
    }
    
    set.seed(seed.data)
    # 1/4 infer parameters using ABC-MCMC (produces also ABC-P)
    ###########################################################
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
    n.discr <- 6 # for 'pred2' summary only
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
    } else if (any(opt.abc[[a]]$summary==c('orig','pred','pred2'))) {
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
        } else if (opt.abc[[a]]$summary=='pred2') {
          s <- c(s, xy[,(n-n.discr+1):n])
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
    } else if (opt.abc[[a]]$summary=='pred2') {
      sdim.not <- n.discr*sum(xy.ind)
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
  
  # 2/4 Compute the approx. true result in 'mis.data' and one obs. popul. case
  # This can be expensive so we precompute it here.
  ############################################################################
  if(1 && inf.task[length(opt.abc)+1] && t2 == 0 && obs.noise.stdev == 0 && !all(xy.ind)) {
    # One observed population case
    if (opt.plt$pr) {cat('Computing approx. true result...\n')}
    set.seed(seed.inf)
    res.atrue <- part.data.approx.true(xy0, theta.true, t.all, n, obs, xy.ind, obs.noise.stdev, part.opt, opt.plt$q)
    
    save(res.atrue, obs, opt, file = part.opt$fn.samples)
    if (opt.plt$pr) {cat('Computation done!\n')}
    
  } else if (inf.task[length(opt.abc)+1] && t2 > 0 && obs.noise.stdev == 0) {
    # 'Missing data' case
    if (opt.plt$pr) {cat('Computing approx. true result...\n')}
    set.seed(seed.inf)
    res.atrue <- mis.data.approx.true(obs$xyobs[,n], theta.true, t.end, t.all[ind.pred], t.all[ind.obs2[1]], 
                                      obs, ind.obs2, obs.noise.stdev, mis.opt, opt.plt$q)
    save(res.atrue, obs, opt, file = mis.opt$fn.samples)
    if (opt.plt$pr) {cat('Computation done!\n')}
  }
  
  
  # Compute the rest of the ABC results (fast) and then plot/analyze them
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
    
    # Get the (approx.) 'true' predictive distribution.
    res.true <- NA
    if (1 && t2 == 0 && obs.noise.stdev == 0 && !all(xy.ind)) {
      # One observed population case
      load(file = part.opt$fn.samples) # 'res.atrue' is pre-computed
      res.true <- res.atrue # allows to conveniently use existing plotting code
      cat('max distance:\n')
      print(res.true$d.eps)
      cat('number of accepted simulations:\n')
      print(res.true$n)
      
    } else if (t2 > 0 && obs.noise.stdev == 0) {
      # 'Missing data' case
      load(file = mis.opt$fn.samples) # 'res.atrue' is pre-computed
      res.true <- res.atrue # allows to conveniently use existing plotting code
      cat('max distance:\n')
      print(res.true$d.eps)
      cat('number of accepted simulations:\n')
      print(res.true$n)
      
    } else if (t2 == 0) {
      if (obs.noise.stdev == 0) {
        xy0p <- obs$xyobs[,n] # non-noisy case
      } else {
        xy0p <- obs$xy[,n] # noisy case, implemented also here although xy0p unknown
      }
      res.true <- lv.pred(xy0p, theta.true, t.end, t.all[ind.pred], obs.noise.stdev, opt.plt$q)
    }
    
    # 3/4 ABC with exact predictive simulation (ABC-F)
    ##################################################
    # We use the ABC posterior with 'ordinary' summary -> No separate parameter fitting.
    # This applies for the non-noisy standard prediction case when prey/predator both observed.
    res.abcf <- NA
    if (t2 == 0 && obs.noise.stdev == 0 && all(xy.ind)) {
      xy0p <- obs$xyobs[,n]
      res.abcf <- lv.pred(xy0p, res.abca[[1]]$thetas, t.end, t.all[ind.pred], obs.noise.stdev, opt.plt$q)
    }
    
    # 4/4 ABC with exact predictive simulation and simulated latent populations (ABC-L)
    ###################################################################################
    # We use the ABC posterior with 'predictive' summary -> No separate parameter fitting.
    # This is implemented only for the non-noisy case with one observed population.
    res.abclat <- NA
    if (t2 == 0 && obs.noise.stdev == 0 && !all(xy.ind)) {
      if (xy.ind[1]) {
        xy0p <- rbind(obs$xyobs[1,n], res.abca[[2]]$y.sims[2,n,])
      } else {
        xy0p <- rbind(res.abca[[2]]$y.sims[1,n,], obs$xyobs[2,n])
      }
      res.abclat <- lv.pred(xy0p, res.abca[[2]]$thetas, t.end, t.all[ind.pred], obs.noise.stdev, opt.plt$q)
    }
    
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
    
    # plot prediction error (added for v2)
    lv.plot.pred.err(obs, t.all, ind.obs, ind.obs2, ind.pred, res.abca, res.abcf, res.abclat, res.true, opt)
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

part.data.approx.true <- function(xy0, theta.true, t.all, n, obs, xy.ind, obs.noise.stdev, part.opt, q) {
  # Computes approximately the true prediction when only one of the populations
  # is observed. It is assumed that the true parameter is known and some last 
  # observed populations are matched approximately. 
  
  t.preds <- t.all[(n+1):length(t.all)]
  n.pred <- length(t.preds)
  n.sa <- part.opt$n.samples
  # determine how many of the latest observations are matched to simulations
  #n.discr <- 1 # used initially
  n.discr <- 6
  n.discr <- min(n.discr,n-1)
  ind.discr <- (n-n.discr+1):n
  
  if (n.pred*n.sa > 10^8) {
    stop('Too much memory needed.') # TODO: all simulations currently saved to memory
  }
  x.preds <- matrix(NA,n.pred,n.sa)
  y.preds <- matrix(NA,n.pred,n.sa)
  ds <- rep(NA,n.sa)
  # first simulate all
  for (j in 1:n.sa) {
    # we simulate only at those timepoints that are needed
    preds <- lv.simul(theta.true, xy0, c(t.all[1],t.all[ind.discr],t.preds), obs.noise.stdev)
    ds[j] <- max(abs(preds$xyobs[xy.ind,2:(n.discr+1)] - obs$xyobs[xy.ind,ind.discr]))
    x.preds[,j] <- preds$xyobs[1,(n.discr+2):(n.pred+n.discr+1)] # observed part neglected
    y.preds[,j] <- preds$xyobs[2,(n.discr+2):(n.pred+n.discr+1)]
  }
  # select those simulations with the smallest distances
  n.final <- min(part.opt$n.final, n.sa)
  return(handle.simul.dist(x.preds, y.preds, t.preds, ds, n.final, q))
}

mis.data.approx.true <- function(xy0p, theta.true, t.end, t.preds, t2.start, obs, 
                                 ind.obs2, obs.noise.stdev, mis.opt, q) {
  # Computes approximately the true prediction for the 'mis.data' case. The true
  # parameter is assumed known and the prediction is matched exactly at the 
  # beginning and approximately at the end prediction point. 
  
  n.pred <- length(t.preds)
  n.sa <- mis.opt$n.samples
  if (n.pred*n.sa > 10^8) {
    stop('Too much memory needed.') # TODO: all simulations currently saved to memory
  }
  x.preds <- matrix(NA,n.pred,n.sa)
  y.preds <- matrix(NA,n.pred,n.sa)
  ds <- rep(NA,n.sa)
  # first simulate all
  for (j in 1:n.sa) {
    preds <- lv.simul(theta.true, xy0p, c(t.end,t.preds,t2.start), obs.noise.stdev)
    ds[j] <- max(abs(preds$xyobs[,n.pred+2] - obs$xyobs[,ind.obs2[1]]))
    x.preds[,j] <- preds$xyobs[1,2:(n.pred+1)] # first and last timepoint neglected
    y.preds[,j] <- preds$xyobs[2,2:(n.pred+1)]
  }
  # select those simulations with the smallest distances (at the last timepoint)
  n.final <- min(mis.opt$n.final, n.sa)
  return(handle.simul.dist(x.preds, y.preds, t.preds, ds, n.final, q))
}

handle.simul.dist <- function(x.preds, y.preds, t.preds, ds, n.final, q) {
  # Help function that handles the selection of smallest distances.
  
  d.eps <- sort(ds)[n.final]
  inds <- which(ds<=d.eps) # can actually produce >n.final samples but this is ok
  x.preds <- x.preds[,inds]
  y.preds <- y.preds[,inds]
  res.atrue <- list(d.eps=d.eps, n=length(inds), x.preds=x.preds, y.preds=y.preds)
  
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
  pred <- list(x.preds=x.preds, y.preds=y.preds)
  pred$pred.stats.x <- pred.post.stats(x.preds, q, T, t.preds)
  pred$pred.stats.y <- pred.post.stats(y.preds, q, T, t.preds)
  return(pred)
}

################################################################################
# Plotting etc. functions:

lv.plot.params <- function(res.abca, theta.true, opt, log.th=F) {
  # Plots posterior densities obtained using ABC. Plots either the parameters 
  # or log parameters.
  
  library(latex2exp)
  ins <- c('True param.',TeX(c('ABC-P, $s^{(0)}$','ABC-P, $s^{(1)}$')))
  cols <- c('black','red','blue')
  
  d <- length(theta.true)
  kde <- function(samples) density(samples, adjust = 1.3)
  fn <- file.path(opt$save.loc, paste0(ifelse(log.th,'params_plot_log','params_plot'),opt$fn.ext,'.pdf'))
  pdf(file=fn, width = 2*d, height = 2)
  par(mfrow=c(1,d))
  par(mai=c(0.35,0.15,0.01,0.07), mgp=c(1.8,0.5,0))
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
      leloc <- 'topleft'
      # ad hoc fix for legend placement:
      if (opt$scenario==201) {leloc <- 'topright'}
      legend(x=leloc, inset = c(0.02,0.02), legend=ins[m], col=cols[m], lty=rep(1,sum(m)), bg = "white", cex=0.6)
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

lv.plot.pred.err <- function(obs, t.all, ind.obs, ind.obs2, ind.pred, res.abca, res.abcf, res.abclat, res.true, opt) {
  # Plots prediction error as a function of time. Prints also corresponding mean errors. 
  # Quickly made for v2 of the paper. 
  
  only.print <- F # whether to only print the results and not plot the figure
  crit <- 'ae' # which error criterion
  critn <- 'Abs. error'
  #library(latex2exp)
  #ins <- TeX(c('ABC-P, $s^{(0)}$','ABC-P, $s^{(1)}$','ABC-L, $s^{(1)}$','ABC-F, $s^{(0)}$'))
  cols <- c('red','blue','orange','orange')
  lw <- 1.5
  popn <- c('Prey','Predator')
  
  # first compute errors
  err <- array(NA,c(4,2,length(ind.pred))) # size: #ABC methods x #populations x #timepoints
  for (id in 1:2) {
    err[1,id,] <- pred.dens.err(res.abca[[1]], res.true, ind.pred, id, crit)
    err[2,id,] <- pred.dens.err(res.abca[[2]], res.true, ind.pred, id, crit)
    err[3,id,] <- pred.dens.err(res.abclat, res.true, ind.pred, id, crit)
    err[4,id,] <- pred.dens.err(res.abcf, res.true, ind.pred, id, crit)
  }
    
  # print computed errors
  mean.err <- round(rowMeans(err),2) # averaged over populs and timepoints
  nms <- c('ABC-P (s0)', 'ABC-P (s1)', 'ABC-L (s1)', 'ABC-F (s0)')
  names(mean.err) <- nms
  cat('\n')
  print(crit)
  print(mean.err)
  
  # print also stdevs e.g. at the first prediction point
  if (1) {
    #id.stdev <- 2 # which population
    fut.pt.ind <- 1
    #fut.pt.ind <- 60
    for (id.stdev in 1:2) {
      stdevs <- rep(NA,5) # stdevs e.g. at the first pred point
      stdevs[1] <- pred.dens.stdev(res.true, ind.pred, id.stdev)[fut.pt.ind]
      stdevs[2] <- pred.dens.stdev(res.abca[[1]], ind.pred, id.stdev)[fut.pt.ind]
      stdevs[3] <- pred.dens.stdev(res.abca[[2]], ind.pred, id.stdev)[fut.pt.ind]
      stdevs[4] <- pred.dens.stdev(res.abclat, ind.pred, id.stdev)[fut.pt.ind]
      stdevs[5] <- pred.dens.stdev(res.abcf, ind.pred, id.stdev)[fut.pt.ind]
      stdevs <- round(stdevs,2)
      names(stdevs) <- c('true pred',nms)
      cat('\n')
      print(paste0('stdevs of predictive densities, popul. ',id.stdev,':'))
      print(stdevs)
    }
  }
  
  # plot prey/predator as different columns
  if (only.print) {
    return()
  }
  fn <- file.path(opt$save.loc, paste0('pred_plot_err',opt$fn.ext,'.pdf'))
  pdf(file=fn, width = 6, height = 2.5)
  par(mfrow=c(1,2))
  for (id in 1:2) {
    par(mai=c(0.5,0.6,0.05,0.05), mgp=c(1.5,0.5,0))
    xl <- c(t.all[ind.pred[1]], t.all[ind.pred[length(ind.pred)]])
    yl <- range(err[,,], na.rm = T)
    yla <- paste0(critn,' (',popn[id],' popul.)')
    plot(NULL, NULL, xlab = 't', ylab = yla, xlim = xl, ylim = yl)
    #title(main = pop[id])
    for (i in 1:4) {
      if (all(!is.na(err[i,id,]))) {
        lines(t.all[ind.pred], err[i,id,], type = 'l', col = cols[i], lty = 'solid', lwd=lw)
      }
    }
  }
  dev.off()
}

pred.dens.err <- function(res, res.true, ind.pred, id, crit='ae') {
  # Computes the error between 'true' predictive density and corresponding ABC density. 
  # NOTE: Computing of absolute error 'ae' currently only implemented. 
  
  if (is.na(res) || is.null(res) || is.na(res.true) || is.null(res.true)) {
    return(NA)
  }
  # ad-hoc fix to handle ABC-P case where also predictions for observations included
  n.pred <- length(ind.pred)
  n.pred.res <- length(res$pred.stats.x$med)
  inds <- 1:n.pred
  if (n.pred < n.pred.res) {
    inds <- ind.pred
  }
  if (id == 1) {
    return(abs(res$pred.stats.x$med[inds]-res.true$pred.stats.x$med))
  }
  return(abs(res$pred.stats.y$med[inds]-res.true$pred.stats.y$med))
}

pred.dens.stdev <- function(res, ind.pred, id) {
  # Returns the computed stdevs of the predictive densities.
  if (is.na(res) || is.null(res)) {
    return(NA)
  }
  n.pred <- length(ind.pred)
  inds <- 1:n.pred
  if (n.pred < length(res$pred.stats.x$stdev)) {
    inds <- ind.pred
  }
  if (id == 1) {
    return(res$pred.stats.x$stdev[inds])
  }
  return(res$pred.stats.y$stdev[inds])
}


