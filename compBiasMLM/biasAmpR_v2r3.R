require(lattice)
require(multcomp)
require(nlme)
require(arm)
require(foreign)
require(plyr)
library(reshape)


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
    #currently works on lme & lmer objects; could deprecate the lme part
    #variances:
    if (class(lmeObject)=="lmerMod") {
       vcomps <- c(as.numeric(VarCorr(lmeObject)[[1]]),sigma(lmeObject)^2)
    } else vcomps <- as.numeric(VarCorr(lmeObject)[,1]) #variances in the first column
    #vcv is in the log
    list(vcomps=vcomps)
}

#helper function to get s.e. of param ests.
seExtract <- function(lmerObj) coef(summary(lmerObj))[,"Std. Error"]


scaleAll <- function(x) {
    
  # try to scale factors whenever possible

  if (is.numeric(x)) return(scale(x))
  if (is.factor(x)) return(as.numeric(levels(x)[x]))
  return(x) #else it's probably character data
}


processData <- function(X,Z,W,Wfmla,YZdata,id.old,id.new,center,stdz) {
    
    #group mean centering, stdzing, etc.  Everything you meed to do to put the data in form corresp. to the DGP as described in the ObsStudies paper

    #helper function; inelegant, but gets the job done
    #Could be replaced by Map()
    apply.tapply <- function(m,MARGIN,INDEX,FUN,...){
        #call apply on tapply, but keep the arguments straight
        f<-function(m,INDEX) {
            tapply(m,INDEX,FUN=FUN)
        }
        apply(m,MARGIN,FUN='f',INDEX=INDEX)
    }
    
    if (stdz) { #standardize everything (including indicators, means, etc.)
        dnx.hold <- dimnames(X)
        X<-apply(X,2,scale)
        dimnames(X) <- dnx.hold
        dnz.hold <- dimnames(Z)
        Z<-apply(as.matrix(Z),2,scale)
        dimnames(Z) <- dnz.hold
        dnw.hold <- dimnames(W)
        W<-apply(W,2,scale)
        dimnames(W) <- dnw.hold
        idx <-match(id.old,dimnames(YZdata)[[2]])
        id.hold <- YZdata[,idx]
        YZdata<-apply(YZdata,2,scaleAll)
        YZdata[,idx] <-id.hold #don't gmc id
    }
    X.mns <- apply.tapply(X,2,id.new,'mean')
    idx.mns <- match(id.new,as.numeric(dimnames(X.mns)[[1]]))
    X.mn <- X.mns[idx.mns,,drop=F]
    newNmsX <- paste(dimnames(X)[[2]],".mn",sep="")
    dimnames(X.mn)[[2]] <- newNmsX
    Z.mns <- apply.tapply(Z,2,id.new,'mean')
    idx.mns <- match(id.new,as.numeric(dimnames(Z.mns)[[1]])) #repeat, just in case
    Z.mn <- Z.mns[idx.mns,,drop=F]
    newNmsZ <- paste(dimnames(Z)[[2]],".mn",sep="")
    dimnames(Z.mn)[[2]] <- newNmsZ
    Znew <- cbind(Z - Z.mn,Z.mn)
    Xfmla <- as.formula(paste0(c("~",paste(dimnames(X)[[2]],collapse="+")),collapse=""))
    if (center) {
        Xnew <- cbind(X - X.mn,X.mn)
        if (!is.null(Wfmla)) {
            Wfmla <- as.formula(paste0(c(Wfmla,"+",paste(newNmsX,collapse="+")),collapse=""))
        } else {
            Wfmla <- as.formula(paste0(c("~",paste(newNmsX,collapse="+")),collapse=""))
        }
    } else {
        Xnew <- X
    }
    return(list(X=X,Xnew=Xnew,newNmsX=newNmsX,Xfmla=Xfmla,Z=Z,Znew=Znew,newNmsZ=newNmsZ,W=W,Wfmla=Wfmla, YZdata=YZdata))
}

#helper funct. for checking results of try() on lmer calls, which can fail:
isErrClass <- function(x) {class(x)=="try-error"}

runModels <- function(outcome=~Y,treatment = ~Z, level1.pred = ~X1+X2, level2.pred = ~W, group = ~id, data, center=T, stdz=T) {
    
    # all formula components are given using formula notation
    # group identifier
    # data: the dataframe
    #center: group mean center (making this an option may lead to unstable results -- perhaps add a warning?)
    
    ##
    ## Initially written to run both lmer and lme; need to remove redundant code for lme and just use lmer
    ##
    
    nRecs <- dim(data)[1]
    IDname <- as.character(group)[2]
    id <- data[,IDname]
    
    #it is important (really just if using a bootstrap, but we will do it to be consistent) that all id's appear together.  So we sort by id here:
    #could be replaced by more elegant command from the tidyverse?
    ord <- order(id)
    data <- data[ord,] #everything needs to be reordered
    id <- data[,IDname] #refresh this value

    uniq.ids <- sort(unique(id)) #sorting is important at this stage - for tapply & GMC to work properly
    nUniqIds <- length(uniq.ids)
    
    idNew <- match(id,uniq.ids) # map back onto original dimensions
    
    
    #formula set up
    groupStr0 <- as.formula(paste("~1|",IDname,sep=""))
    groupStr <- as.formula("~1|idNew")  #this will override the group ID

    Yname <- as.character(outcome)[2]
    Yfmla <-as.formula(paste(Yname,"~.",sep=""))
    #for lmer:
    Yfmla.lmer <-as.formula(paste(Yname,"~.+(1|idNew)",sep=""))
    Yfmla1.lmer <-as.formula(paste(Yname,"~1+(1|idNew)",sep=""))

    Zname <- as.character(treatment)[2]
    #for lmer:
    Zfmla.lmer <-as.formula(paste(Zname,"~.+(1|idNew)",sep=""))
    Zfmla1.lmer <-as.formula(paste(Zname,"~1+(1|idNew)",sep=""))
    YZdata <- data[,c(Yname,Zname,IDname)] #this will be useful in simpler models.
    
    
    #design matrix/initial data processing
    X <- model.matrix(level1.pred,data)[,-1,drop=F] # first col is just 1s.
    W <- NULL
    if (!is.null(level2.pred)) W <- model.matrix(level2.pred,data)[,-1,drop=F]
    Z <- model.matrix(treatment,data)[,-1,drop=F] # first col is just 1s.
    varZ <- var(Z) #need this value to rescale treatment effect bias
    
    id.cts <- table(idNew)
    allIds <- rep(1:nUniqIds,id.cts)

    #store orig:
    X0<-X
    Z0<-Z
    W0<-W
    YZdata0<-as.matrix(YZdata)
    dimnames(X0) <- list(NULL,dimnames(X0)[[2]]) #to avoid the dup names problem
    dimnames(Z0) <- list(NULL,dimnames(Z0)[[2]]) #to avoid the dup names problem
    dimnames(W0) <- list(NULL,dimnames(W0)[[2]]) #to avoid the dup names problem
    dimnames(YZdata0) <- list(NULL,dimnames(YZdata0)[[2]]) #to avoid the dup names problem

    #on Y~Z and Y~Z+Xctrd+Xmean+W, run OLS.
    #for all group-varying preds: GMC everything (perhaps at the level of the design matrix for ease)
    #so matrix is [Y,Zctrd,Zmean,Xctrd,Xmean,W], where W are group const
    #assumes at least ONE X, potentially NO Ws...
 
 pd <- processData(X0,Z0,W0,Wfmla=level2.pred,YZdata0, IDname,idNew,center,stdz)  #results in a list data structure with most of the hard data manip accomplished
 #rename terms to make identification (by us) clearer as we fit models
    X<-pd$X
    Z<-pd$Z
    W<-pd$W
    Xnew <- pd$Xnew
    newNmsX <- pd$newNmsX
    Znew <- pd$Znew
    Xfmla <- pd$Xfmla
    Znew <- pd$Znew
    newNmsZ <- pd$newNmsZ
    Wfmla <- pd$Wfmla
    YZdata <- pd$YZdata
    
    #fit models needed for parameters used in the paper:

    mlm0.y <- lmer(formula(terms(Yfmla.lmer,data=Znew,allowDotAsName=T)), data=as.data.frame(cbind(cbind(YZdata[,c(Yname,IDname)],Znew),idNew)))
    mlm1.y <- lmer(formula(terms(Yfmla.lmer,data=cbind(Znew,Xnew,W),allowDotAsName=T)), data=as.data.frame(cbind(cbind(YZdata[,c(Yname,IDname)],Znew,cbind(Xnew,W)),idNew)))
    mlm0.z <- lmer(Zfmla1.lmer,data=as.data.frame(cbind(YZdata,idNew)))
    mlm1.z <- lmer(formula(terms(Zfmla.lmer, data=cbind(Xnew,W),allowDotAsName=T)), data=as.data.frame(cbind(cbind(YZdata,Xnew,W),idNew)))
        #need the next run for bounds, but not for system of eqns.
    ols0.y <- lm( formula(terms(Yfmla,data=Z,allowDotAsName=T)),data=as.data.frame(YZdata))
    ols1.y <- lm( formula(terms(Yfmla,data=cbind(Z,Xnew,W),allowDotAsName=T)), data=as.data.frame(cbind(YZdata,Xnew,W))) #cbind to W to avoid fail on bind to data.frame
    mlm.ucm.y <- lmer(Yfmla1.lmer,data=as.data.frame(YZdata))



    #data used in runs (for any post processing)
    newData <- as.data.frame(cbind(YZdata[,c(Yname,IDname)],Znew,cbind(Xnew,W),idNew))
    #construct V(X) and V(W) - easier to do here instead of in postproc phase:
    varW <- as.matrix(var(model.matrix(Wfmla,newData)[,-1,drop=F]))
    varX <- as.matrix(var(model.matrix(Xfmla,newData)[,-1,drop=F]))

    list(newData=newData,ols0.y=ols0.y,ols1.y=ols1.y, mlm.ucm.y=mlm.ucm.y,mlm0.y=mlm0.y,mlm1.y=mlm1.y, mlm0.z=mlm0.z,mlm1.z=mlm1.z,outcome=outcome,treatment =  treatment,treatment.gp=newNmsZ,level2.pred =Wfmla,varX=varX,varW=varW,varZ=varZ,Xfmla=Xfmla,Wfmla=Wfmla,Yfmla=Yfmla.lmer,Zfmla=Zfmla1.lmer,stdz=stdz)
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
    
     data<-runModelRslt$newData #some vars are added to this data.frame, so best to use it instead of original data
     level2.pred <- runModelRslt$level2.pred
     treatment <- runModelRslt$treatment
     outcome <- runModelRslt$outcome
     #get gamm.z'V(w)gamm.z from model 1 (for Z)
     xprods <- makeSpecProds(level2.pred,fixef(runModelRslt$mlm1.z),fixef(runModelRslt$mlm1.y),runModelRslt$varW,runModelRslt$varX)
     varW <- xprods$varW #these are reordered to reflect proper variable ordering.
     varX <- xprods$varX
     
     outcome.varname <- attr(terms(outcome),"term.labels")
     varY <- var(data[,outcome.varname])
     #note: treatment.gp is already extracted from the formula form... could be changed to be consistent w/ otgher passing protocol.
     treatment.gpname <- runModelRslt$treatment.gp
     treatment.varname <- attr(terms(treatment),"term.labels")
     varZ <- var(data[,treatment.varname])
 
     #WITHIN:
     #there need to be two extract pgms.  ONE does fixefs and se, other does VC and lvcv...
     #BETWEEN:
     #most effects are extracted.  just need between est of tau (more elegant would be to reuse lmeExtract, but would be redundant and add comp time):
     
     #add s.e. as well
     extr0 <- lmeExtract(runModelRslt$mlm0.y,treatment.varname,treatment.gpname)
     vcv.tau0 <- extr0$vcv.taus
     tau0.w<-extr0$tau[1]
     tau0.b<-extr0$tau[2]
     #
     extr1 <- lmeExtract(runModelRslt$mlm1.y,treatment.varname,treatment.gpname)
     vcv.tau1 <- extr1$vcv.taus
     tau1.w<-extr1$tau[1]
     tau1.b<-extr1$tau[2]

     within.bias.diff <- tau0.w-tau1.w  # based on lme, fixef
     between.bias.diff <- tau0.b-tau1.b  #based on lme
     
     #OLS
     extr <- lmExtract(runModelRslt$ols0.y,treatment.varname)
     tau0.ols <- extr$tau
     se.tau0.ols <- extr$se.tau
     extr <- lmExtract(runModelRslt$ols1.y,treatment.varname)
     tau1.ols <- extr$tau
     se.tau1.ols <- extr$se.tau
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
     list(sds.y.ucm=list(sd.eps.y.ucm=sd.eps.y.ucm,sd.alpha.y.ucm=sd.alpha.y.ucm),sigma.sq.between.y.0=sigma.sq.between.y.0,sigs=sigs, bias.diffs=c(within.bias.diff,between.bias.diff,ols.bias.diff),tau.w=c(tau0.w,tau1.w),tau.b=c(tau0.b,tau1.b),tau.ols=c(runModelRslt$ols0.y$coef[2],runModelRslt$ols1.y$coef[2]),bndProdList=xprods$bndProdList, varY=varY,varZ=varZ,vcv.tau0=vcv.tau0,vcv.tau1=vcv.tau1,seOLS=c(se.tau0.ols, se.tau1.ols))
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
              return(sqrt(gz.vw.gz)*max(abs(tau.max)*sqrt(gz.vw.gz),sqrt(gy.vw.gy),c.b.y))
          } else {
              return(sqrt(gy.vw.gy*gz.vw.gz))
          }
      }
    } else {
        f <- function(tau.max,bound=T) {
            if (bound) {
                return(sqrt(bz.vx.bz)*max(abs(tau.max)*sqrt(bz.vx.bz),sqrt(by.vx.by),c.w.y))
            } else {
                return(sqrt(by.vx.by*bz.vx.bz))
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
    parmRange <- scaleFactor*bnd.f(tau.max)

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
olsMwin <- function(zeta,delta,cW=cW,cB=cB) {abs((zeta+delta)/(cW+cB))-abs(zeta/cW)}
olsBiasOnly <- function(zeta,delta,cW=cW,cB=cB) {(zeta+delta)/(cW+cB)}
winBiasOnly <- function(zeta,cW=cW) {zeta/cW}

correctedTau.o <- function(tau.o,zeta,delta,cW,cB) {
  return(tau.o - olsBiasOnly(zeta,delta,cW,cB))
}


##
## there are confPts, confPtsCol, and perhaps other params that we've deprecated.
## clean code of that redundancy.
##

##ACTION: change size of labels in this plot... (is CEX= enough?)
zdPlot <- function(zeta1,delta1,parmRange,rescaleParms=c(1,1),confPts=NULL,confPtsCol=8,targetVals=c(0),targetPch=c(0),taus=NULL,offset=5,cW=sqrt(5),cB=sqrt(5),n.pts=9,N=201, autoAdjZeta=F, autoAdjDelta=F, autoAdjProbs=c(.1,.9),zetaRange=NULL,deltaRange=NULL,cex=1,zInflator=1,...) {
    #the function that plots "danger zones"
    #use targetVals to show where 'upper bd' might be or where eta wd be if we had unbiased Y eqn.
    #recale sd parms (y.w,y.b,z.w,z.b) defaults to no rescale.  O/w feed s.d.s based on '0' models for z,y(ucm)
    
    pfunct = function(x, y, z,cols0,target,targetPch,confPts,confPtsCol,taus,cex,...) {
        panel.levelplot(x, y, z,...)
        if (!is.null(confPts)) panel.points(x=confPts[,1],y=confPts[,2],pch='.',col=confPtsCol)
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
    
    if (autoAdjZeta && !is.null(confPts)) {
        qtl_5_95 <- quantile(confPts[,2],prob=autoAdjProbs,na.rm=T)
        zeta <- seq(qtl_5_95[1],qtl_5_95[2],length=N)
    } else {
        zeta <- seq(min(zeta1),max(zeta1),length=N) #default
    }
    if (!is.null(zetaRange)) { #override everything
        zeta <- seq(zetaRange[1],zetaRange[2],length=N)
    }
    
    ##
    ## uses confPts - anything with bounds can be removed.
    ##
    
    if (autoAdjDelta && !is.null(confPts)) {
        qtl_5_95 <- quantile(confPts[,1],prob=autoAdjProbs,na.rm=T)
        delta <- seq(qtl_5_95[1],qtl_5_95[2],length=N)
    } else {
        delta <- seq(min(delta1),max(delta1),length=N) #default
    }
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
        corr.tau.o <- taus$ols - olsBiasOnly(zeta[locs],delta[locs],cW=cW,cB=cB)
        corr.tau.w <- taus$win - winBiasOnly(zeta[locs],cW=cW)
        tau.switch <- sign(olsMwin(zeta[locs],delta[locs],cW=cW,cB=cB))<=0
        cat(corr.tau.o,"\n")
        cat(corr.tau.w,"\n")
        tau.to.plot <-list(taus=cbind(corr.tau.o,corr.tau.w),offset=offset*(delta[2]-delta[1]),switch=tau.switch)
    } else {
        tau.to.plot <- NULL #for passing to levelplot pfunct
    }
    df <- data.frame(x=rep(delta,N),y=rep(zeta,each=N),z=c(t(bdiff)))
    
    #alt version of cols0=grey(1-abs(parmRange)/max(abs(parmRange)))
    levelplot(z~x*y,data=df,at=at.pts2,colorkey=list(at=at.pts,labels=list(at=at.pts,label=spec.lbl,cex=.75*cex)), panel=pfunct,row.values=delta,column.values=zeta,xlab=list(label=expression(delta^{yz}),cex=.85*cex),ylab=list(label=expression(zeta^{yz}),cex=.85*cex),zeta1=zeta1,delta1=delta1,cols0=8,target=targetPts,targetPch=targetPch,confPts=confPts,confPtsCol=confPtsCol,taus=tau.to.plot,cex=cex,scales=list(x=list(cex=.85*cex),y=list(cex=.85*cex),cex=cex),...)
}



lwDelete <- function(fmla,data) {
    opt0 <- options()
    options(na.action="na.pass")
    mm <- model.matrix(fmla,data)
    bb<-apply(is.na(mm),1,sum)==0
    options(opt0)
    data[bb,]
}



scaledColors <- function(x,colorScheme=cm.colors,nPerSide,rev=T) {
    
    nCols <- 2*nPerSide+1
    pal <- colorScheme(nCols)
    if (rev) pal <- rev(pal) #colder colors are positive?
    fact <- max(max(abs(x)),min(abs(x)))
    newx <- x/fact #still has a true zero  -- and one endpt has abs(x)==1
    idx <- round(1+nPerSide*(1+newx))
    ticks <- seq(-fact,fact,length=nCols)
    ticks[ticks<min(x)] <- NA #take these out
    ticks[ticks>max(x)] <- NA #take out
    index <- sort(unique(idx)) #these are used
    return(list(cols=pal[idx],ticks=ticks,palette=pal[index],index=index))
}

