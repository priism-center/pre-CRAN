require(lattice)
require(multcomp)
require(nlme)
require(arm)
require(foreign)
require(plyr)
library(reshape)
library(MASS)
library(cluster)


#IMPLEMENTATION:
#  1. for all group-varying preds: GMC everything (perhaps at the level of the design matrix for ease)
#     so matrix is [Y,Zctrd,Zmean,Xctrd,Xmean,W], where W are group const
#assumes ONE X, potentially NO Ws...
#  2. set up 4 types of models:
#     Y~Zctrd+Zmean
#     Y~Zctrd+Zmean+Xctrd+Xmean+W
#     Z~1
#     Z~Xctrd+Xmean+W
#
#     this yields tau0.ols, tau1.ols
#  4. on both Y models, above, run MLM (with +(1|id) or random=~1|id)
#     this yields tau0.w, tau1.w, tau0.b, tau1.b, because we have group mean ctrg.
#  5. on both Z models, above, run MLM (with +(1|id) or random=~1|id)
#     this yields sigma.sq.bet[1:2] (really should be RENAMED, as it is within, between, for the model with no covars); the with covars model yields zd.rho, which is the reduction in var as you move to more complex models.  and note correspondance to c.b0, c.b1, c.w0, c.w1... - easier to identify
#   6. run Y~1+(1|id) to get var.alpha.y.ucm.
#   7. ols.bias.diff is orig - new (with covars)
#

#We want the following functions:
#runModels:
#    fmla: the full formula for the Y model regression (not MLM), Y~Z+X...
#    all as formulas:
#    outcome=~Y
#    treatment = ~Z
#    level1.pred = ~X
#    level2.pred = ~W
#    group = ~id
#    data: the dataframe

#helper fnc:
rm.na <- function(x) x[!is.na(x)]

#some of the next 3 functions return lists with one object.  This is due to prior need for returning multiple objs and the desire not to recode everything that depends on them.  Returned objs could be simplified if dependencies handled.
lmExtract <- function(lmObject,treatName) {
    #extract tau (coef) corresponding to treatment Z & s.e.(tau)
    # based only on 'treatName' entry
    tau <- coef(lmObject)[treatName]
    rse <- summary(lmObject)$sigma
    vcv0 <- summary(lmObject)$cov
    se <- (sqrt(diag(vcv0))*rse)[treatName]
    list(tau=tau,se.tau=se)
}

lmeExtract <- function(lmeObject,treatName,treatGpName) {
    #works on lmer obj
    #extract fixef of treatment (two forms - 2nd one is group mean)
    #fixef:
    tau <- fixef(lmeObject)[c(treatName,treatGpName)]
    #se
    list(tau=tau)
}


lmeExtract.varcomp <- function(lmeObject) {
    #extract variances of random intercept model & ln(vcv) of those:
    #currently works on lmer only objects
    #variances:
    vcomps <- c(as.numeric(VarCorr(lmeObject)[[1]]),sigma(lmeObject)^2)
    list(vcomps=vcomps)
}

#R function to mimic one prior in STAN.
my_rbeta <-function(n,a,b,scale) {
    scale*rbeta(n,a,b)*sample(c(1,-1),n,replace=T)
}

#helper function to get s.e. of param ests.
seExtract <- function(lmerObj) coef(summary(lmerObj))[,"Std. Error"]

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

scaleAll <- function(x) {
    
  # try to scale factors whenever possible

  if (is.numeric(x)) return(scale(x))
  if (is.factor(x)) return(as.numeric(levels(x)[x]))
  return(x) #else it's probably character data
}


printResults <- function(mdlFit,mdlName,digits=3,debug=FALSE) {
  # helper function to give the output from the models
  
  #extract variance comps & some bias diffs
  pObj <- extractParams(mdlFit)
  #used in ObsStudies paper:
  print(paste("Intermediary Model Fits for: ",mdlName))
  print("Multilevel model fit:")
  print(summary(mdlFit$mlm1.y))
  print("OLS regression model fit:")
  print(summary(mdlFit$ols1.y))
  print(paste("Table 1 for: ",mdlName))
  rsltTab1 <- rbind(cbind(pObj$tau.w,pObj$tau.b,pObj$tau.ols),pObj$bias.diffs)
  dimnames(rsltTab1) <- list(c("tau0","tau1","diff"),c("Within","Between","OLS"))
  print(round(rsltTab1,digits=digits))
  print("ICCs for models:")
  print(round(cbind(pObj$sigs,pObj$sigs[,2]/apply(pObj$sigs,1,sum)),digits=digits))
  if (debug) print(paste("gy,vy,t.gz ",paste(round(sqrt(c(pObj$bndProdList$gy.vw.gy,pObj$sds.y.ucm$sd.alpha.y.ucm^2,2^2*pObj$bndProdList$gz.vw.gz)),digits=digits),collapse=", "),sep=" "))
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


plotResults <- function(singleRun,fit.sims,dat,priorChoice=1L, nGridPoints=401, tau.max=1, 
                        priorSampSize=20000, gpSize=Inf,seed=432101234, confEllipse=100,
                        delta.plot.range=1, zeta.plot.range=1,offset=0,forceMax=NULL, refLine=T,refPoints=T) {
    
    plotParams <- prePlotParams(singleRun, nGridPoints=nGridPoints, tau.max=tau.max, gpSize=gpSize)
    par(mar=oldpar$mar+c(2,3,2,4),bg="grey98")
    
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
    
    if (refPts) {  #these are the ref pts we want...
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
    
    colorbar.plot(dat$scale^2-.18-offset,dat$scale^2-offset,colSet$index,col=colSet$palette,adj.y=1,strip.width = 0.06,strip.length=6*(.06),horizontal=F)
    len <- length(colSet$index)
    skipEvery <- 2
    ticks <- round(rev(colSet$ticks),1)
    frac <- 7.95*.06 #approx size of strip
    for(i in 1:len) {
        if ((i-1)%%skipEvery==0) mtext(ticks[i],side=4,line=0.2,at=dat$scale^2-offset-0.14-frac*(i-1)/len,padj=1,las=1,outer=F,cex=1)
    }
    mtext(expression(tau),side=3,line=+0.3,at=dat$scale^2-.18-offset,cex=1.75)
}

combine.formula <- function(fmla1,fmla2) {
  #works on RHS formulas, e.g., ~x1+x2...
  terms1 <- attr(terms.formula(fmla1),"term.labels")
  terms2 <- attr(terms.formula(fmla2),"term.labels")
  as.formula(paste0("~",paste(c(terms1,terms2),collapse="+")))
}

geoMeanGroupSize <- function(groupID) {prod(table(groupID))^(1/length(unique(groupID)))}

prepDat <- function(outcome, treatment, level1.pred, level2.pred, group, data) {

  renumber <- function(x) {
    #assume repeats, renumber to 1:length(unique(...))
    match(x, sort(unique(x)))
  }  
  groupName <- as.character(rhs(group))
  outcomeName <- as.character(rhs(outcome))
  treatName <- as.character(rhs(treatment))
  
  #among other things, forces id to be 1:M, not whatever it is in raw data.
  
  N <- dim(data)[1]
  x <- as.array(model.matrix(level1.pred,data)[,-1]) # drop const
  y <- data[,outcomeName]
  z.ctrd <- data[,treatName]
  z.mn <- data[,paste0(treatName,".mn")]
  w <- as.array(model.matrix(level2.pred,data)[,-1]) #drop const.
  p <- dim(x)[2]
  q <- dim(w)[2]
  id <- renumber(data[,groupName])
  M <- length(unique(id))
  gpSize <- geoMeanGroupSize(id)
  #fix z; current data version is centered, and needs to not be, going forward
  list(N=N,M=M,x=x,y=y,z=z.ctrd+z.mn,w=w,p=p,q=q,id=id,gpSize=gpSize)
}

runModels <- function(outcome=~Y,treatment = ~Z, level1.pred = ~X1+X2, level2.pred = ~W, group = ~id, data) {
    
    # all formula components are given using formula notation
    # group identifier
    # data: the dataframe
    # already group mean centered ...
    
    IDname <- as.character(group)[2]

    #formula set up
    Yname <- as.character(rhs(outcome))
    Zname <- as.character(rhs(treatment))
    IDname <- as.character(rhs(group))
    fmlaAll <- combine.formula(level1.pred,level2.pred)
    #for lm:
    Yfmla1.lm <- as.formula(paste0(Yname,fmlaAll,"+I(",Zname,"+",Zname,".mn)")) #Z should be uncentered
    Yfmla0.lm <- as.formula(paste0(Yname,"~+I(",Zname,"+",Zname,".mn)"))  #ditto
    #for lmer:
    Yfmla1.lmer <-as.formula(paste0(Yname,fmlaAll,"+",Zname,"+",Zname,".mn+(1|",IDname,")"))
    Yfmla0.lmer <-as.formula(paste0(Yname,"~",Zname,"+",Zname,".mn+(1|",IDname,")"))
    YfmlaUCM.lmer <-as.formula(paste0(Yname,"~1+(1|",IDname,")"))
    Zfmla1.lmer <-as.formula(paste0("I(",Zname,"+",Zname,".mn)",fmlaAll,"+(1|",IDname,")")) #model 1
    Zfmla0.lmer <-as.formula(paste0("I(",Zname,"+",Zname,".mn)","~1+(1|",IDname,")")) #model 0
    #fit models needed for parameters used in the paper:
    mlm0.y <- lmer(Yfmla0.lmer,data=data)
    mlm1.y <- lmer(Yfmla1.lmer,data=data)
    mlm0.z <- lmer(Zfmla0.lmer, data=data)
    mlm1.z <- lmer(Zfmla1.lmer, data=data)
        #need the next run for bounds, but not for system of eqns.
    ols0.y <- lm(Yfmla0.lm,data=data)
    ols1.y <- lm(Yfmla1.lm,data=data)
    mlm.ucm.y <- lmer(YfmlaUCM.lmer,data=data)
    
    #even though it's not req'd to run stan after runModels, this is useful for stan call and a few variance calcs, so we do it here and pass the results forward.
    pDat <- prepDat(outcome, treatment, level1.pred, level2.pred, group, data)

    #construct V(Z), V(X) and V(W) - easier to do here instead of in postproc phase:
    varW <- as.matrix(var(pDat$w))
    varX <- as.matrix(var(pDat$x))
    varZ <- var(pDat$z) #may need this value to rescale treatment effect bias
    list(ols0.y=ols0.y,ols1.y=ols1.y, mlm.ucm.y=mlm.ucm.y,mlm0.y=mlm0.y,mlm1.y=mlm1.y, mlm0.z=mlm0.z,mlm1.z=mlm1.z, 
         outcome=outcome,treatment=treatment, level2.pred=level2.pred,varX=varX,varW=varW,varZ=varZ, data=data, preppedData=pDat)
}

makeSpecProds <- function(level2.pred,fe_mlm1.z,fe_mlm1.y,varW,varX) {
   #helper function to extractParams to generate terms needed eventually for bounds
   #get gamm.z'V(w)gamm.z from model 1 (for Z)
   #add sd*se to gamma or beta [ref]
   
   wNames <- attr(terms(level2.pred),"term.labels")
   nmsZ <- names(fe_mlm1.z)
   if (is.null(nmsZ)) nmsZ <- dimnames(fe_mlm1.z)[[2]] #handles 2 cases
   idxW <- match(wNames,nmsZ)
   xNames <- names(fe_mlm1.z)[-idxW][-1] #drop 'intercept'
   idxX <- match(xNames,nmsZ)
   if (class(fe_mlm1.z)=="data.frame") {
       gamm.z <- as.matrix(fe_mlm1.z[,idxW])
       beta.z <- as.matrix(fe_mlm1.z[,idxX])
   } else{ #assume single row
       gamm.z <- t(as.matrix(fe_mlm1.z[idxW,drop=F]))
       beta.z <- t(as.matrix(fe_mlm1.z[idxX,drop=F]))
   }
   
   #next bit is based on biased parms, but worth doing for comparison purposes (and bounds):
   
   nmsY <- names(fe_mlm1.y)
   if (is.null(nmsY)) nmsY <- dimnames(fe_mlm1.y)[[2]] #handles 2 cases
   idxW <- match(wNames,nmsY)
   idxX <- match(xNames,nmsY) #need to exclude treatment from these.
   if (class(fe_mlm1.y)=="data.frame") {
       gamm.y <- as.matrix(fe_mlm1.y[,idxW])
       beta.y <- as.matrix(fe_mlm1.y[,idxX])
   } else{ #single row
       gamm.y <- t(as.matrix(fe_mlm1.y[idxW,drop=F]))
       beta.y <- t(as.matrix(fe_mlm1.y[idxX,drop=F]))
   }
   #fix cols of var mat to align with gamma & beta
   idxW <- match(wNames,dimnames(varW)[[1]])
   idxX <- match(xNames,dimnames(varX)[[1]])
   varW <- as.matrix(varW[idxW,idxW])
   varX <- as.matrix(varX[idxX,idxX])
  
   bndProdList <- list(gz.vw.gz = gamm.z%*%varW%*%t(gamm.z), bz.vx.bz = beta.z%*%varX%*%t(beta.z), gy.vw.gz = gamm.y%*%varW%*%t(gamm.z), gy.vw.gy = gamm.y%*%varW%*%t(gamm.y), by.vx.bz = beta.y%*%varX%*%t(beta.z), by.vx.by=beta.y%*%varX%*%t(beta.y))
  
   list(varX=varX,varW=varW,bndProdList=bndProdList)
   
}


extractParams <- function(runModelRslt) {
    
    #extract certain ests of var comps & bias diffs & taus (extraneous)
    #for use in sensitivity plots (and diagnostics)
    
     data<-runModelRslt$data #some vars are added to this data.frame, so best to use it instead of original data
     level2.pred <- runModelRslt$level2.pred
     treatment <- runModelRslt$treatment
     outcome <- runModelRslt$outcome
     #get gamm.z'V(w)gamm.z from model 1 (for Z)
     xprods <- makeSpecProds(level2.pred,fixef(runModelRslt$mlm1.z),fixef(runModelRslt$mlm1.y),runModelRslt$varW,runModelRslt$varX)
     varW <- xprods$varW #these are reordered to reflect proper variable ordering.
     varX <- xprods$varX
     
     outcome.varname <- attr(terms(outcome),"term.labels")
     varY <- var(data[,outcome.varname])
     #note: treatment.gp is already extracted from the formula form... could be changed to be consistent w/ other passing protocol.
     treatment.varname <- attr(terms(treatment),"term.labels")
     treatment.gpname <-paste0(treatment.varname,".mn")
     varZ <- var(data[,treatment.varname])
 
     #WITHIN:
     #there need to be two extract pgms.  ONE does fixefs and se, other does VC and lvcv...
     #BETWEEN:
     #most effects are extracted.  just need between est of tau (more elegant would be to reuse lmeExtract, but would be redundant and add comp time):
     
     #add s.e. as well
     extr0 <- lmeExtract(runModelRslt$mlm0.y,treatment.varname,treatment.gpname)
     tau0.w<-extr0$tau[1]
     tau0.b<-extr0$tau[2]
     #
     extr1 <- lmeExtract(runModelRslt$mlm1.y,treatment.varname,treatment.gpname)
     tau1.w<-extr1$tau[1]
     tau1.b<-extr1$tau[2]

     within.bias.diff <- tau0.w-tau1.w  # based on lme, fixef
     between.bias.diff <- tau0.b-tau1.b  #based on lme
     
     #OLS
     #more complicated to find correct coef to extract.  need to find the `I(treat+treat.mn)` entry.
     treatConstrName<- attr(terms(runModelRslt$ols0.y),"term.labels") #there is only one predictor in this model & it's the treatment.
    
     extr <- lmExtract(runModelRslt$ols0.y,treatConstrName)
     tau0.ols <- extr$tau
     extr <- lmExtract(runModelRslt$ols1.y,treatConstrName)
     tau1.ols <- extr$tau
     #
     ols.bias.diff <- tau0.ols - tau1.ols
     
     #####
     #varcorr:
     extr.v0 <- lmeExtract.varcomp(runModelRslt$mlm0.z)
     sigma.sq.within.0 <- extr.v0$vcomps[2]
     sigma.sq.between.0 <- extr.v0$vcomps[1]
     #
     extr.v1 <- lmeExtract.varcomp(runModelRslt$mlm1.z)
     sigma.sq.within.1 <- extr.v1$vcomps[2]
     sigma.sq.between.1 <- extr.v1$vcomps[1]
     #
     extr.y0 <- lmeExtract.varcomp(runModelRslt$mlm0.y)
     sigma.sq.between.y.0 <- extr.y0$vcomps[1]
     
     sigs <- matrix(c(sigma.sq.within.0,sigma.sq.within.1,sigma.sq.between.0,sigma.sq.between.1),2,2,byrow=F)

     #these are used for bounds, so sampling variability of these is less crucial:
     ucm.vc <- lmeExtract.varcomp(runModelRslt$mlm.ucm.y)
     sd.alpha.y.ucm <- sqrt(ucm.vc$vcomps[1])
     sd.eps.y.ucm <- sqrt(ucm.vc$vcomps[2])
    
    list(sds.y.ucm=list(sd.eps.y.ucm=sd.eps.y.ucm, sd.alpha.y.ucm=sd.alpha.y.ucm),
     sigma.sq.between.y.0=sigma.sq.between.y.0,sigs=sigs,
     bias.diffs=c(within.bias.diff,between.bias.diff,ols.bias.diff),
     tau.w=c(tau0.w,tau1.w),tau.b=c(tau0.b,tau1.b),
     tau.ols=c(tau0.ols,tau1.ols),
     bndProdList=xprods$bndProdList,varY=varY,varZ=varZ)
}

makeBnds <- function(paramObj,param='gamma') {
    
    #extract from paramObj - yields a function of tau that sets a bound on feasible eta;
    #param allows one to switch the 'open' param (see paper) and this yields a different set of bounds.

    # for convenience, extract these here:

    gz.vw.gz <- paramObj$bndProdList$gz.vw.gz
    gy.vw.gy <- paramObj$bndProdList$gy.vw.gy
    c.b.y <- paramObj$sds.y.ucm$sd.alpha.y.ucm
    bz.vx.bz <- paramObj$bndProdList$bz.vx.bz
    by.vx.by <- paramObj$bndProdList$by.vx.by
    c.w.y <- paramObj$sds.y.ucm$sd.eps.y.ucm
    
    if (param=='gamma') {
      f <- function(tau.max,bound=T) {
          if (bound) {
              return(as.numeric(c(sqrt(gz.vw.gz)*max(abs(tau.max)*sqrt(gz.vw.gz),sqrt(gy.vw.gy),c.b.y))))
          } else {
              return(as.numeric(sqrt(gy.vw.gy*gz.vw.gz)))
          }
      }
    } else {
        f <- function(tau.max,bound=T) {
            if (bound) {
                return(as.numeric(c(sqrt(bz.vx.bz)*max(abs(tau.max)*sqrt(bz.vx.bz),sqrt(by.vx.by),c.w.y))))
            } else {
                return(as.numeric(sqrt(by.vx.by*bz.vx.bz)))
            }
        }
    }
    return(f)
}


makeLinEqMat <- function(c_w0,c_w1,c_b0,c_b1,idx=4,gpSize=Inf) {
    
    #set up the matrix of the 4 linear equations.  See paper.
    m <- matrix(0,4,4)
    m[1,1] <- 1/c_w0-1/c_w1
    m[1,3] <- 1/c_w0
    
    if (is.finite(gpSize)) {
        m[2,1] <- (1/gpSize)/(c_w0/gpSize+c_b0) - (1/gpSize)/(c_w1/gpSize+c_b1)
        m[2,2] <- 1/(c_w0/gpSize+c_b0) - 1/(c_w1/gpSize+c_b1)
        m[2,3] <- (1/gpSize)/(c_w0/gpSize+c_b0)
        m[2,4] <- 1/(c_w0/gpSize+c_b0)
    } else {
        m[2,2] <- 1/c_b0-1/c_b1
        m[2,4] <- 1/c_b0
    }
    m[3,1] <- m[3,2] <- 1/(c_w0+c_b0)-1/(c_w1+c_b1)
    m[3,3] <- m[3,4] <- 1/(c_w0+c_b0)
    m[4,idx] <- 1
    m
}

recover <- function(paramObj,varyParm="gamma",bnd.f,nParms=201,tau.max=2,gpSize=Inf) {
    
    #solves a system of equations for 3 out of 4 of (zeta_1, delta_1,byxbz,gywgz) in terms of the 4th, based on values for sigs and differences in bias for 3 estimation methods
    #sigs: 2x2 matrix rows are (c.w0,c.b0)
    #                          (c.w1,c.b1)
    #deltaBiases: bias(w,b,ols) or c(tau0.w-tau1.w,tau0.b-tau1.b,tau0.ols-tau1.ols)
    
    #one param must be 'open' - zeta_1, delta_1,etc.  varyParm choose it by naming it from parms
    #IMPT: the BOUNDS are only derived for gamma (gywgz) as the omitted param.
    
    #extract from paramObj:
    sigs <- paramObj$sigs
    deltaBiases <-paramObj$bias.diffs
    
    parms <- c("zeta","delta","beta","gamma")

    idx <- match(varyParm,parms)
    m <- makeLinEqMat(sigs[1,1],sigs[2,1],sigs[1,2],sigs[2,2],idx,gpSize)
    
    #generate parmRange
    scaleFactor <- seq(-1,1,length=nParms)
    parmRange <- scaleFactor*as.numeric(bnd.f(tau.max))

    #iterate through the 4th (omitted) param
    len <- length(parmRange)
    r <- matrix(NA,len,4) # result
    delts <- c(deltaBiases,0) #placeholder for 3 deltas and one open param
    #error trapping
    if (class(try(condNum<-kappa(m,exact=T)))=="try-error") print("unable to compute condition number\n")
    mInv <- solve(m)
    for (i in 1:len) {
        delts[4] <- parmRange[i] #assume we know the omitted param.
        r[i,] <- mInv%*%delts #solve for 3 other params, given components and omitted param.
    }
    dimnames(r)[[2]] <- parms
    #zetaDeltaMat are the zetas and deltas consistent with the data; parmRange gives the open param for the given row of zds
    #ALT cond num:
    ##m <- makeLinEqMat(sigs[1,1],sigs[2,1],sigs[1,2],sigs[2,2],idx=3,gpSize)
    #if (class(try(condNum<-kappa(m,exact=T)))=="try-error") print("unable to compute condition number\n")
    list(zetaDeltaMat=r,plausible.taus=scaleFactor*tau.max,parmRange=parmRange,condNum=condNum)
}


rescaleSetVal <- function(extParms) {
    #set values for rescaling based on fitted models
    c(sqrt(extParms$varY),sqrt(extParms$varZ))
}


#helper functs (mostly for zdPlot; useful in other contexts

#functions that compute asymptotic absolute bias differences for different estimators
betMwin <- function(zeta,delta,cB=1,cW=1) {abs(delta/cB)-abs(zeta/cW)}
olsMwin <- function(zeta,delta,cB=1,cW=1) {abs((zeta+delta)/(cW+cB))-abs(zeta/cW)}
glsMwin <- function(zeta,delta,cB=1,cW=1,lambda=.5) {abs((zeta+lambda*delta)/(cW+lambda*cB))-abs(zeta/cW)}

winBias <- function(zeta,delta,cB=1,cW=1) {zeta/cW}
betBias <- function(zeta,delta,cB=1,cW=1,gpSize=Inf) {
    delta/(cW/gpSize+cB)  #relies on 1/Inf == 0    
}
olsBias <- function(zeta,delta,cB=1,cW=1) {(zeta+delta)/(cW+cB)}
glsBias <- function(zeta,delta,cB=1,cW=1,lambda=.5) {(zeta+lambda*delta)/(cW+lambda*cB)}

correctedTau.o <- function(tau.o,zeta,delta,cW,cB) {
  return(tau.o - olsBias(zeta,delta,cW,cB))
}


prePlotParams <- function(mdlFit,nGridPoints=201,tau.max=1,gpSize=Inf) {
   ################################
   # Bounds calcs 
   # params needed for plots
   # cond numb
   ################################
   # nGridPoints: number of points to build the confounder line.
   # tau.max: abs(tau) used in the bounds calculation (plausible value supplied by user)
   # gpSize: for balanced designs, the number of subjects/group.  For others, try geom. mean or Inf.
   ################################
  
   #FIND THE RIGHT open Param from the gamma- or beta-based value for eta.
   pObj <- extractParams(mdlFit)
   # corresponds to the bound when gamma is the open param:
   recovParmsGammaBased <- recover(pObj,varyParm="gamma",makeBnds(pObj,param="gamma"),nParm=nGridPoints,tau.max=tau.max, gpSize=gpSize) 
   # beta-based equiv.:
   recovParmsBetaBased <- recover(pObj,varyParm="beta",makeBnds(pObj,param="beta"),nParm=nGridPoints,tau.max=tau.max,gpSize=gpSize)
   #get the gamma corresp. to this beta by searching through the list of evaluated points and finding closest corresponding.
   #this is one of the target values to report
   gammaYZ.at.beta <- recovParmsBetaBased$zetaDeltaMat[which.min(abs(recovParmsBetaBased$zetaDeltaMat[,"beta"]- as.numeric(pObj$bndProdList$by.vx.bz))),"gamma"]
   #get other 'target' values:
   gammaYZ.at.zeta.0 <- recovParmsGammaBased$zetaDeltaMat[which.min(abs(recovParmsGammaBased$zetaDeltaMat[,"zeta"])),"gamma"]
   gammaYZ.at.delta.0 <- recovParmsGammaBased$zetaDeltaMat[which.min(abs(recovParmsGammaBased$zetaDeltaMat[,"delta"])),"gamma"]
   gammaYZ.CS <- pObj$bndProdList$gy.vw.gz  #bound derived from C-S Ineq.
   
   #reduce range if info from Beta warrants it
   nu.g.min <- min(recovParmsBetaBased$zetaDeltaMat[,"gamma"])
   nu.g.max <- max(recovParmsBetaBased$zetaDeltaMat[,"gamma"])
   #adjust range here
   b.in.range <- recovParmsGammaBased$zetaDeltaMat[,"gamma"] >= nu.g.min & recovParmsGammaBased$zetaDeltaMat[,"gamma"] <= nu.g.max
   #by this point, it might be better not to call it "gamma based"
   recovParmsGammaBased$zetaDeltaMat <- recovParmsGammaBased$zetaDeltaMat[b.in.range,]
   recovParmsGammaBased$parmRange <- recovParmsGammaBased$parmRange[b.in.range]
   recovParmsGammaBased$plausible.taus <- recovParmsGammaBased$plausible.taus[b.in.range]
   
   bndVals <- c(gammaYZ.at.beta,gammaYZ.CS,gammaYZ.at.zeta.0,gammaYZ.at.delta.0)
   names(bndVals) <- c("gYZbeta","gYZcs","gYZzeta","gYZdelta")
   
   list(zetaDeltaMat=recovParmsGammaBased$zetaDeltaMat,parmRange=recovParmsGammaBased$parmRange, 
        plausible.taus=recovParmsGammaBased$plausible.taus, condNum=recovParmsGammaBased$condNum,bndVals=bndVals )
}



##ACTION: change size of labels in this plot... (is CEX= enough?)
zdPlot <- function(zeta1,delta1,parmRange,rescaleParms=c(1,1),targetVals=c(0),targetPch=c(0),taus=NULL,offset=5,cW=sqrt(5),cB=sqrt(5),n.pts=9,N=201, autoAdjZeta=F, zetaRange=NULL,deltaRange=NULL,cex=1,zInflator=1,debug=F,...) {
    #the function that plots "danger zones"
    #use targetVals to show where 'upper bd' might be or where eta wd be if we had unbiased Y eqn.
    #recale sd parms (y.w,y.b,z.w,z.b) defaults to no rescale.  O/w feed s.d.s based on '0' models for z,y(ucm)
    #NOTE: since between ests aren't used in the plots, the gpSize is never used by this function and so is not a param.
    
    pfunct = function(x, y, z,cols0,target,targetPch,taus,cex,...) {
        panel.levelplot(x, y, z,...)
        panel.abline(h=0,col.line=1)
        panel.abline(v=0,col.line=1)
        panel.points(x=delta1,y=zeta1,pch=16,cex=1*cex,col.symbol=cols0)
        panel.points(x=target[,1],y=target[,2],pch=targetPch,lwd=3,cex=1.4*cex,col=1)
        if (!is.null(taus)) {
            len <- length(taus$taus[,1])
            for (i in 1:len) {
                if (!is.na(taus$switch[i])) {
                    if (taus$switch[i]) {
                        panel.text(x=target[i,1]+taus$offset,y=target[i,2],label=bquote(tilde(tau)[o] == .(round(taus$taus[,1],2)[i])),adj=0,cex=.75*cex)
                    } else panel.text(x=target[i,1]-taus$offset,y=target[i,2],label=bquote(tilde(tau)[w] == .(round(taus$taus[,2],2)[i])),adj=1,cex=.75*cex)
                }
            }
        }
    }
    
    zeta <- seq(min(zeta1),max(zeta1),length=N) #default
    if (!is.null(zetaRange)) { #override everything
        zeta <- seq(zetaRange[1],zetaRange[2],length=N)
    }
    
    delta <- seq(min(delta1),max(delta1),length=N) #default
    if (!is.null(deltaRange)) { #override
        delta <- seq(deltaRange[1],deltaRange[2],length=N)
    }
    bdiff <-outer(zeta,delta,FUN="olsMwin",cW=cW,cB=cB)*zInflator
    
    rng.bd <- max(abs(min(bdiff)),abs(max(bdiff)))
    at.pts <-round(seq(-rng.bd,+rng.bd,length=n.pts),2)
    at.pts2 <- seq(-rng.bd,+rng.bd,length=1+2*(n.pts-1))
    
    spec.lbl <- as.character(at.pts)
    spec.lbl[1] <- "OLS\nbetter"
    spec.lbl[length(spec.lbl)] <- "Within\nbetter"
    
    #use sigma for y,z, b & w, to scale zeta & delta product params
    zeta <- zeta/prod(rescaleParms[1])
    delta <- delta/prod(rescaleParms[2])
    #need to rescale zeta1 & delta1 too....
    zeta1 <- zeta1/prod(rescaleParms[1])
    delta1 <- delta1/prod(rescaleParms[2])
    #locate target point (for beta or gamma):
    nLocs <- length(targetVals)
    locs <- rep(NA,nLocs) #initialize
    gap <- parmRange[2]-parmRange[1]
    for (i in 1:nLocs) {
      locs[i] <- which.min(abs(parmRange-targetVals[i]))
      if (locs[i]==1 && targetVals[i] < parmRange[1]-gap) {
          locs[i] <- NA #boundary cond - drop
          cat("Below minimum range, eta=",targetVals[i],"point beyond",x[1],y[1])
      }
      if (locs[i]==N && targetVals[i] > parmRange[1]+gap) {
          locs[i] <- NA #boundary cond - drop
          cat("Above minimum range, eta=",targetVals[i],"point beyond",x[N],y[N])
      }
    }
    targetPts <- matrix(cbind(delta1[locs],zeta1[locs]),nLocs,2,byrow=F) #byrow=F for when you pass more than one point...
    #here is where we get vals for plausible taus for the plot - evaluate at targetVals
    if (!is.null(taus)) { #these are ests. of tau based on zeta,delta
        corr.tau.o <- taus$ols - olsBias(zeta[locs],delta[locs],cW=cW,cB=cB)
        corr.tau.w <- taus$win - winBias(zeta[locs],cW=cW)
        tau.switch <- sign(olsMwin(zeta[locs],delta[locs],cW=cW,cB=cB))<=0
        if (debug) {
          cat("Bias-corrected (model-based) tau evaluated at specified points (least abs bias indicated on plot):\n")
          dstr <- cbind(zeta[locs],delta[locs],corr.tau.o,corr.tau.w,tau.switch)
          dimnames(dstr) <- list(NULL,c("delta","zeta","OLS-based","Within-based","OLS-better?"))
          print(round(dstr,digits=2))
        }
        tau.to.plot <-list(taus=cbind(corr.tau.o,corr.tau.w),offset=offset*(delta[2]-delta[1]),switch=tau.switch)
    } else {
        tau.to.plot <- NULL #for passing to levelplot pfunct
    }
    df <- data.frame(x=rep(delta,N),y=rep(zeta,each=N),z=c(t(bdiff)))
    
    #alt version of cols0=grey(1-abs(parmRange)/max(abs(parmRange)))
    levelplot(z~x*y,data=df,at=at.pts2,colorkey=list(at=at.pts,labels=list(at=at.pts,label=spec.lbl,cex=.75*cex)), panel=pfunct,row.values=delta,column.values=zeta,xlab=list(label=expression(delta^{yz}),cex=.85*cex),ylab=list(label=expression(zeta^{yz}),cex=.85*cex),zeta1=zeta1,delta1=delta1,cols0=8,target=targetPts,targetPch=targetPch,taus=tau.to.plot,cex=cex,scales=list(x=list(cex=.85*cex),y=list(cex=.85*cex),cex=cex),...)
}



lwDelete <- function(fmla,data) {
    opt0 <- options()
    options(na.action="na.pass")
    mm <- model.matrix(fmla,data)
    bb<-apply(is.na(mm),1,sum)==0
    options(opt0)
    data[bb,]
}

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

idealizedParms <- function(varX,varW,tau,betaY,betaZ,gammaY,gammaZ,zetaY,zetaZ,deltaY,deltaZ,varAlphaY,varAlphaZ,varEpsilonY,varEpsilonZ) {
    
    #scalar x,w only, for now.
    #Betw. Group Size -> infty for now
    c_W0 <- betaZ^2*varX + zetaZ^2 + varEpsilonZ
    c_B0 <- gammaZ^2*varW + deltaZ^2 + varAlphaZ
    c_W1 <- c_W0 - betaZ^2*varX
    c_B1 <- c_B0 - gammaZ^2*varW
    
    covWinYZ <- betaZ*varX*betaY + zetaZ*zetaY
    covBetYZ <- gammaZ*varW*gammaY + deltaZ*deltaY
    
    sd.eps.y.ucm <- sqrt(tau^2*(c_W0) + betaY^2*varX + zetaY^2 + varEpsilonY + 2*tau*covWinYZ)
    sd.alpha.y.ucm <- sqrt(tau^2*(c_B0) + gammaY^2*varW + deltaY^2 + varAlphaZ + 2*tau*covBetYZ)
    
    zetaYZ <- zetaY*zetaZ
    deltaYZ <- deltaY*deltaZ
    
    ##replace w/ function calls?
    
    winBias0 <- (betaZ*varX*betaY+zetaYZ)/(betaZ^2*varX+c_W1)  #corrected 5Feb15; orig bZ(vX)bY - wrong.
    betBias0 <- (gammaZ*varW*gammaY+deltaYZ)/(gammaZ^2*varW+c_B1) #fixed gammaZ*varW*gammaY to proper denom as per text  -- group size -> Infty
    olsBias0 <- (betaZ*varX*betaY+gammaZ*varW*gammaY+zetaYZ+deltaYZ)/(betaZ^2*varX+gammaZ^2*varW+c_W1+c_B1) #same correction
    
    winBias1 <- zetaYZ/c_W1
    betBias1 <- deltaYZ/c_B1
    olsBias1 <- (zetaYZ+deltaYZ)/(c_W1+c_B1)
    sigs <- matrix(c(c_W0,c_W1,c_B0,c_B1),2,2,byrow=F)
    list(sigs=sigs,biasDiffs=c(winBias0-winBias1,betBias0-betBias1,olsBias0-olsBias1),winBias=c(winBias0,winBias1),betBias=c(betBias0,betBias1),olsBias=c(olsBias0,olsBias1), sd.y.ucm=c(sd.alpha.y.ucm,sd.eps.y.ucm))
}
