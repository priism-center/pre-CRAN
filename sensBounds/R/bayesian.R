# Bayesian-only functions

#---------HELPER FUNCTIONS-----------

#R function to mimic one prior in STAN.
my_rbeta <-function(n,a,b,scale) {
  scale*rbeta(n,a,b)*sample(c(1,-1),n,replace=T)
}

#naming conventions/helper functions:

makeFileName <- function(priorChoice=0L,countryNumber=1,prefix="bayesV4R0prior",suffix="Rdata",countryNames=c("italy","scotland", "sweden","finland","usa","england")) {
  if (countryNumber!=0) {
    prefixAdj <- paste(countryNames[countryNumber],prefix,sep=".",collapse="")
  } else prefixAdj <- prefix
  paste(paste(prefixAdj,as.character(priorChoice),sep=""),suffix,sep=".",collapse="")
}

#core code to sample from different product priors:
makeProdPrior <- function(n,kernel=normal,seed=432101234,...) {
  set.seed(seed)
  prior_zeta <-kernel(n,...)*kernel(n,...)
  prior_delta <-kernel(n,...)*kernel(n,...)
  list(prior_zeta=prior_zeta,prior_delta=prior_delta)
}

priorSample <- function(n,priorChoice=1L,seed=432101234,...) {

  if (priorChoice==0L) {
    samps <- makeProdPrior(n=n,kernel=my_rbeta,seed=seed,a=dat$a,b=dat$b,scale=dat$scale)
  } else if (priorChoice==1L) {
    samps <- makeProdPrior(n=n,kernel=rnorm,seed=seed, mean=0,sd=dat$sd)
  } else if (priorChoice==2L) {
    samps <- makeProdPrior(n=n,kernel=runif,seed=seed,min=-dat$lim,max=+dat$lim)
  }
  samps
}


#---------PLOTTING-----------

#helper for plots
scaledColors <- function(x,colorScheme=cm.colors,nPerSide,forceMax=NULL,rev=T) {

  #simplified version in which max value determines symmetric range
  # or is overridden by forceMax.  the NAs that dropped colors in unused range was removed.
  # See older versions to recover.

  nCols <- 2*nPerSide+1
  pal <- colorScheme(nCols)
  if (rev) pal <- rev(pal) #colder colors are positive?
  fact <- max(max(abs(x)),min(abs(x)))
  if (!is.null(forceMax)) fact <- forceMax #override to make a common scale!
  newx <- x/fact #still has a true zero  -- and one endpt has abs(newx)==1 when force==NULL
  #top/bottom code if abs(newx)>1
  fix.idx <- abs(newx)>1
  newx[fix.idx] <- sign(newx)[fix.idx] #replace w/ +/- 1
  idx <- round(1+nPerSide*(newx+1))
  ticks <- seq(-fact,fact,length=nCols)
  return(list(cols=pal[idx],ticks=ticks,palette=pal,index=1:nCols))
}

plotResults <- function(singleRun,fit.sims,dat,priorChoice=1L, nGridPoints=401, tau.max=1,
                        priorSampSize=20000, gpSize=Inf,seed=432101234, confEllipse=100,
                        delta.plot.range=1, zeta.plot.range=1,offset=0,forceMax=NULL, refLine=T,refPoints=T) {

  if (!is.null(singleRun)) plotParams <- prePlotParams(singleRun, nGridPoints=nGridPoints, tau.max=tau.max, gpSize=gpSize)

  samps <- priorSample(n=priorSampSize,priorChoice=priorChoice,seed=seed,zeta_sd=dat$zeta_sd, delta_sd=dat$delta_sd,zeta_lim=zeta_lim, delta_lim=delta_lim, a=dat$a,b=data$b,scale=dat$scale)

  prior_zeta <- samps$prior_zeta
  prior_delta <- samps$prior_delta

  plot(prior_delta,prior_zeta, pch=1,cex=.15,col="darkgrey",xlim=c(-delta.plot.range+.2+offset,delta.plot.range-.2-offset),ylim=c(-zeta.plot.range+.2+offset,zeta.plot.range-.2-offset), xlab=expression(delta^yz), ylab=expression(zeta^yz), cex.lab=1.75, cex.axis=1.75)

  colSet <- scaledColors(fit.sims$b_y_z,cm.colors,forceMax=forceMax,nPerSide=10)
  cols <- colSet$cols
  cmNeg <- cm.colors(21)[1]
  cmPos <- cm.colors(21)[21]

  #plot posterior
  points(fit.sims$delta_yz,fit.sims$zeta_yz,pch=1,cex=.5, col=cols)
  #draw ref. lines
  abline(v=0,col=1);abline(h=0,col=1); abline(a=0,b=1,col=1); abline(a=0,b=-1,col=1)

  if (refLine) {
    points(plotParams[[1]][,2],plotParams[[1]][,1],pch=19,cex=.4,col="salmon")
  }

  if (refPoints) {  #these are the ref pts we want...
    xNames <- dimnames(dat$x)[[2]]
    wNames <- dimnames(dat$w)[[2]]
    deltaIdx <-  grep(".mn",wNames,fixed=TRUE) #this should yield preds introduced by us.  If there is an X lacking between group variation, we don't handle it well at all.
    zetaIdx <- 1:length(xNames)
    zetaVals <- apply(fit.sims$b_y_x*fit.sims$b_z_x,2,mean)  #within
    deltaVals <- apply(fit.sims$b_y_w*fit.sims$b_z_w,2,mean)  #between
    #the next line could have "null" cases as needs more error trapping.
    obsDeltaZetas <- rbind(cbind(deltaVals[deltaIdx],zetaVals[zetaIdx]),cbind(deltaVals[-deltaIdx],0))
    points(obsDeltaZetas,pch=5, col="green",cex=1.5,lwd=2.5)
  }
  if (confEllipse<100 & confEllipse>0) { #draw a conf. ellipse.
    samples <- cbind(fit.sims$delta_yz,fit.sims$zeta_yz)
    # Finding the confEllipse% highest density / minimum volume ellipse
    fit <- cov.mve(samples, quantile.used = nrow(samples) * confEllipse / 100)
    points_in_ellipse <- samples[fit$best, ]
    ellipse_boundary <- predict(ellipsoidhull(points_in_ellipse))
    # Plotting it
    lines(ellipse_boundary, col="lightgreen", lwd=4)
    #legend("topleft", paste0(confEllipse,"%"), col = "lightgreen", lty = 1, lwd = 2,bty='n')
  }

  colorbar.plot(delta.plot.range-.18-offset,zeta.plot.range-offset,colSet$index,col=colSet$palette,adj.y=1,strip.width = 0.06,strip.length=6*(.06),horizontal=F)
  len <- length(colSet$index)
  skipEvery <- 2
  ticks <- round(rev(colSet$ticks),1)
  frac <- 7.95*.06 #approx size of strip
  for(i in 1:len) {
    if ((i-1)%%skipEvery==0) mtext(ticks[i],side=4,line=0.2,at=delta.plot.range-offset-0.14-frac*(i-1)/len,padj=1,las=1,outer=F,cex=1)
  }
  mtext(expression(tau),side=3,line=+0.3,at=dat$scale^2-.18-offset,cex=1.75)
}





