require(lattice)
require(multcomp)
require(nlme)
require(arm)
require(foreign)
require(plyr)
library(reshape)
#NOTES:
#    Attempt to clean up code on Oct. 5, 2015
#
# 1. Dropping bootstrap in this version.
# 2. Implement lmer in parallel to lme to allow resampling
# 3: wrapper for it all - single call, lots of params.


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

#R funtion to mimic one prior in STAN.
my_rbeta <-function(n,a,b,scale) {
    scale*rbeta(n,a,b)*sample(c(1,-1),n,replace=T)
}

lmExtract <- function(lmObject,treatName) {
    #extract 'Z' beta & s.e.(beta); only 'treatname' entry
    tau <- coef(lmObject)[treatName]
    rse <- summary(lmObject)$sigma
    vcv0 <- summary(lmObject)$cov
    se <- (sqrt(diag(vcv0))*rse)[treatName]
    list(tau=tau,se.tau=se)
}

lmeExtract <- function(lmeObject,treatName,treatGpName) {
    #works on lme or lmer obj
    #extract fixef of treatment & variances of random intercept model & ln(vcv) of those:
    #fixef:
    tau <- fixef(lmeObject)[c(treatName,treatGpName)]
    #se
    list(tau=tau)
}


lmeExtract.varcomp <- function(lmeObject) {
    #extract fixef of treatment & variances of random intercept model & ln(vcv) of those:
    #variances:
    if (class(lmeObject)=="lmerMod") {
       vcomps <- c(as.numeric(VarCorr(lmeObject)[[1]]),sigma(lmeObject)^2)
    } else vcomps <- as.numeric(VarCorr(lmeObject)[,1]) #variances in the first column
    #vcv is in the log
    list(vcomps=vcomps)
}

scaleAll <- function(x) { #will try to scale factors whenever possible

  if (is.numeric(x)) return(scale(x))
  if (is.factor(x)) return(as.numeric(levels(x)[x]))
  return(x) #else it's probably character data
}


processData <- function(X,Z,W,Wfmla,YZdata,id.old,id.new,center,stdz) {
    
    #group mean centering, stdzing, etc.  Everything you meed to do to put the data in form corresp. to the DGP

    #helper function; inelegant, but gets the job done
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

#helper funct:
isErrClass <- function(x) {class(x)=="try-error"}

#consider removing all internal bootstrap features, incl. the need to reorder, and any lme (not lmer) calls.  
runModels <- function(outcome=~Y,treatment = ~Z, level1.pred = ~X1+X2, level2.pred = ~W, group = ~id, data, center=T, stdz=T) {
    
    # all formula components are given using formula notation
    # group identifier
    # data: the dataframe
    #center: group mean center
    
    nRecs <- dim(data)[1]
    IDname <- as.character(group)[2]
    id <- data[,IDname]
    #it is important (really just for bootstrap, but we will do it to be consistent) that all id's appear together.  So we sort by id here:
    ord <- order(id)
    data <- data[ord,] #everything needs to be reordered
    id <- data[,IDname] #refresh this value

    uniq.ids <- sort(unique(id)) #sorting is important at this stage - for tapply & GMC to work properly
    nUniqIds <- length(uniq.ids)
    
    id.orig <- match(id,uniq.ids) # map back onto original dimensions
    groupStr0 <- as.formula(paste("~1|",IDname,sep=""))
    groupStr <- as.formula("~1|idNew")  #this will override the group ID

    Yname <- as.character(outcome)[2]
    Yfmla <-as.formula(paste(Yname,"~.",sep=""))
###    Yfmla1 <-as.formula(paste(Yname,"~1",sep=""))
    #for lmer:
    Yfmla.lmer <-as.formula(paste(Yname,"~.+(1|idNew)",sep=""))
    Yfmla1.lmer <-as.formula(paste(Yname,"~1+(1|idNew)",sep=""))

    Zname <- as.character(treatment)[2]
###    Zfmla <-as.formula(paste(Zname,"~.",sep=""))
###    Zfmla1 <-as.formula(paste(Zname,"~1",sep=""))
    #for lmer:
    Zfmla.lmer <-as.formula(paste(Zname,"~.+(1|idNew)",sep=""))
    Zfmla1.lmer <-as.formula(paste(Zname,"~1+(1|idNew)",sep=""))
    YZdata <- data[,c(Yname,Zname,IDname)] #this will be useful in simpler models.
    
    X <- model.matrix(level1.pred,data)[,-1,drop=F] # first col is just 1s.
    W <- NULL
    if (!is.null(level2.pred)) W <- model.matrix(level2.pred,data)[,-1,drop=F]
    Z <- model.matrix(treatment,data)[,-1,drop=F] # first col is just 1s.
    varZ <- var(Z) #need this value to rescale treatment effect bias
    
    id.cts <- table(id.orig)
    allIds <- rep(1:nUniqIds,id.cts)

    #rm.na function may be more efficient than na.omit.
    rm.na <- function(x) x[!is.na(x)]


    #store orig:
    X0<-X
    Z0<-Z
    W0<-W
    YZdata0<-as.matrix(YZdata)
    dimnames(X0) <- list(NULL,dimnames(X0)[[2]]) #to avoid the dup names problem
    dimnames(Z0) <- list(NULL,dimnames(Z0)[[2]]) #to avoid the dup names problem
    dimnames(W0) <- list(NULL,dimnames(W0)[[2]]) #to avoid the dup names problem
    dimnames(YZdata0) <- list(NULL,dimnames(YZdata0)[[2]]) #to avoid the dup names problem

    idNew <- id.orig  ##unnec. but keep for a while (o/w all idNew must become id.orig....)
    id.sel <- 1:nRecs
    #on Y~Z and Y~Z+Xctrd+Xmean+W, run OLS.
    #for all group-varying preds: GMC everything (perhaps at the level of the design matrix for ease)
    #so matrix is [Y,Zctrd,Zmean,Xctrd,Xmean,W], where W are group const
    #assumes ONE X, potentially NO Ws...
 
    pd <- processData(X0,Z0,W0,Wfmla=level2.pred,YZdata0, IDname,idNew,center,stdz)
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
        

    mlm0.y <- lmer(formula(terms(Yfmla.lmer,data=Znew,allowDotAsName=T)), data=as.data.frame(cbind(cbind(YZdata[,c(Yname,IDname)],Znew),idNew)))
    mlm1.y <- lmer(formula(terms(Yfmla.lmer,data=cbind(Znew,Xnew,W),allowDotAsName=T)), data=as.data.frame(cbind(cbind(YZdata[,c(Yname,IDname)],Znew,cbind(Xnew,W)),idNew)))
    mlm0.z <- lmer(Zfmla1.lmer,data=as.data.frame(cbind(YZdata,idNew)))
    mlm1.z <- lmer(formula(terms(Zfmla.lmer, data=cbind(Xnew,W),allowDotAsName=T)), data=as.data.frame(cbind(cbind(YZdata,Xnew,W),idNew)))
        #need the next run for bounds, but not for system of eqns.
    ols0.y <- lm( formula(terms(Yfmla,data=Z,allowDotAsName=T)),data=as.data.frame(YZdata))
    ols1.y <- lm( formula(terms(Yfmla,data=cbind(Z,Xnew,W),allowDotAsName=T)), data=as.data.frame(cbind(YZdata,Xnew,W))) #cbind to W to avoid fail on bind to data.frame
    mlm.ucm.y <- lmer(Yfmla1.lmer,data=as.data.frame(YZdata))



    #construct V(X) and V(W) - easier here:
    newData <- as.data.frame(cbind(YZdata[,c(Yname,IDname)],Znew,cbind(Xnew,W),idNew))
    varW <- as.matrix(var(model.matrix(Wfmla,newData)[,-1,drop=F]))
    varX <- as.matrix(var(model.matrix(Xfmla,newData)[,-1,drop=F]))

    list(newData=newData,ols0.y=ols0.y,ols1.y=ols1.y, mlm.ucm.y=mlm.ucm.y,mlm0.y=mlm0.y,mlm1.y=mlm1.y, mlm0.z=mlm0.z,mlm1.z=mlm1.z,outcome=outcome,treatment =  treatment,treatment.gp=newNmsZ,level2.pred =Wfmla,varX=varX,varW=varW,varZ=varZ,Xfmla=Xfmla,Wfmla=Wfmla,Yfmla=Yfmla.lmer,Zfmla=Zfmla1.lmer,stdz=stdz)
}



makeSpecProds <- function(level2.pred,fe_mlm1.z,fe_mlm1.y,varW,varX) {
   #get gamm.z'V(w)gamm.z from model 1 (for Z)
   #add sd*se to gamma or beta
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
   #this is based on biased parms, but worth doing for comparison purposes (and bounds):
   
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
   gz.vw.gz<-gamm.z%*%varW%*%t(gamm.z)
   bz.vx.bz<-beta.z%*%varX%*%t(beta.z)
   gy.vw.gz<-gamm.y%*%varW%*%t(gamm.z)
   gy.vw.gy<-gamm.y%*%varW%*%t(gamm.y)
   by.vx.bz<-beta.y%*%varX%*%t(beta.z)
   by.vx.by<-beta.y%*%varX%*%t(beta.y)
   ## if (dim(fe_mlm1.z)[1]>1) { #signals passing a sim result in - force a scalar per iteration
   ##   gz.vw.gz<-diag(gz.vw.gz)
   ##   bz.vx.bz<-diag(bz.vx.bz)
   ##   gy.vw.gz<-diag(gy.vw.gz)
   ##   gy.vw.gy<-diag(gy.vw.gy)
   ##   by.vx.bz<-diag(by.vx.bz)
   ##   by.vx.by<-diag(by.vx.by)
   ##}
   list(varX=varX,varW=varW,gz.vw.gz=gz.vw.gz, bz.vx.bz=bz.vx.bz, gy.vw.gz=gy.vw.gz, gy.vw.gy=gy.vw.gy, by.vx.bz=by.vx.bz, by.vx.by=by.vx.by)
}

seExtract <- function(lmerObj) coef(summary(lmerObj))[,"Std. Error"]

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
     gz.vw.gz <- xprods$gz.vw.gz
     bz.vx.bz <- xprods$bz.vx.bz
     gy.vw.gz <- xprods$gy.vw.gz
     gy.vw.gy <- xprods$gy.vw.gy
     by.vx.bz <- xprods$by.vx.bz
     by.vx.by <- xprods$by.vx.by
     
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
     
     #tau0.b<-fixef(runModelRslt$mlm0.y)[runModelRslt$treatment.gp]
     #tau1.b<-fixef(runModelRslt$mlm1.y)[runModelRslt$treatment.gp]
     
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
     #FIX:
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
     
     list(sds.y.ucm=list(sd.eps.y.ucm=sd.eps.y.ucm,sd.alpha.y.ucm=sd.alpha.y.ucm),sigma.sq.between.y.0=sigma.sq.between.y.0,sigs=sigs, bias.diffs=c(within.bias.diff,between.bias.diff,ols.bias.diff),tau.w=c(tau0.w,tau1.w),tau.b=c(tau0.b,tau1.b),tau.ols=c(runModelRslt$ols0.y$coef[2],runModelRslt$ols1.y$coef[2]),gz.vw.gz=gz.vw.gz,gy.vw.gz=gy.vw.gz,gy.vw.gy=gy.vw.gy,bz.vx.bz=bz.vx.bz,by.vx.bz=by.vx.bz,by.vx.by=by.vx.by, varY=varY,varZ=varZ,vcv.tau0=vcv.tau0,vcv.tau1=vcv.tau1,seOLS=c(se.tau0.ols, se.tau1.ols))
}



makeBnds <- function(paramObj,param='gamma') {
    
    #extract from paramObj - yields a function of tau that sets a bound on feasible eta;
    #the bounds=F option will let us look at just the one point that eta will be if bias is low or non-existent.

    gz.vw.gz <- paramObj$gz.vw.gz
    gy.vw.gy <- paramObj$gy.vw.gy
    #    sigma.sq.between.y.0 <- paramObj$sigma.sq.between.y.0
    c.b.y <- paramObj$sds.y.ucm$sd.alpha.y.ucm
    bz.vx.bz <- paramObj$bz.vx.bz
    by.vx.by <- paramObj$by.vx.by
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
    
    #set up the matrix of the 4 linear equations.
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

recover <- function(paramObj,varyParm="gamma",bnd.f,nParms=201,tau.max=2,gpSize=Inf,debug=F) {
    
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
    #needs error trapping
    if (class(try(condNum<-kappa(m,exact=T)))=="try-error") print("unable to compute condition number\n")
    mInv <- solve(m)
    for (i in 1:len) {
        delts[4] <- parmRange[i] #assume we know the omitted param.
        r[i,] <- mInv%*%delts #solve for 3 other params, given components and omitted param.
    }
    if (debug) browser()
    dimnames(r)[[2]] <- parms
    #zetaDeltaMat are the zetas and deltas consistent with the data; parmRange gives the open param for the given row of zds
    #ALT cond num:
    ##m <- makeLinEqMat(sigs[1,1],sigs[2,1],sigs[1,2],sigs[2,2],idx=3,gpSize)
    #if (class(try(condNum<-kappa(m,exact=T)))=="try-error") print("uanble to compute condition number\n")
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


##ACTION: change size of labels in this plot... is CEX= enough???
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
    
    if (autoAdjDelta && !is.null(confPts)) {
        qtl_5_95 <- quantile(confPts[,1],prob=autoAdjProbs,na.rm=T)
        delta <- seq(qtl_5_95[1],qtl_5_95[2],length=N)
    } else {
        delta <- seq(min(delta1),max(delta1),length=N) #default
    }
    if (!is.null(deltaRange)) { #override everything
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

rm.na <- function(x) x[!is.na(x)]

bsConfParmsSim <- function(modelRun,nSims=100,nTarg=201,conf=.95,bnd.f,varyParm="gamma",gpSize,tau.max=1,debug=F,allowVarBool=rep(T,7),lines=T) {
    
    #use s.e. and vcvs in pObj to resample key confouding params
    if (class(modelRun$mlm0.z)!="lmerMod") stop ("Rerun runModels using lmeSwitch=F\n")
    #extract varcomps:
    mlm1.z.sim<-sim(modelRun$mlm1.z,nSims)
    c_b1s <- apply(mlm1.z.sim@ranef[[1]][,,1],1,var)
    c_w1s <- mlm1.z.sim@sigma^2
    #extract bias diffs:
    mlm1.y.sim<-sim(modelRun$mlm1.y,nSims)
    outcome.name <- attr(terms(modelRun$outcome),"term.labels")
    treat.gpname <- modelRun$treatment.gp
    treat.varname <- attr(terms(modelRun$treatment),"term.labels")
    tau1.w <- mlm1.y.sim@fixef[,treat.varname]
    tau1.b <- mlm1.y.sim@fixef[,treat.gpname]
    #now ols:
    lm1.y.sim<-sim(modelRun$ols1.y,nSims)
    tau1.ols <- lm1.y.sim@coef[,treat.varname]
    
    #assuming indep across models yields too liberal a bound -- resample from '1' models for '0' model fits.
    data <- modelRun$newData
    ncols <- dim(data)[2]
    nrecs <- dim(data)[1]

    #setup for tau0.w and tau0.b, tau0.o
    xwCoefs <- mlm1.y.sim@fixef
    #drop first 2, and last cols of data, which are outcome and identifiers, resp.; add col. of ones, proceed.
    fePred <- xwCoefs%*%rbind(1,t(data[,-c(1:2,ncols)]))
    newId <- data[,ncols] #seems like specialized knowledge at this point...
    randEffs <- mlm1.y.sim@ranef[[1]]
    id.re <- as.numeric(dimnames(randEffs)[[2]]) #effectively, we require numeric ids..
    idx <- match(newId,id.re)
    rePred <- randEffs[,idx,1] # asms random intercept model (dim: nsims X nrecs)
    sigs <- mlm1.y.sim@sigma
    #build synthetic data from model 1 as null:
    eps <- matrix(rnorm(nSims*nrecs,mean=0,sd=rep(sigs,each=nrecs)),byrow=T,nrow=nSims)
    newYs <- fePred + rePred + eps # synthetic data, all to get a tau.
    tmpdata.y <- data[,-2] #drop id (old)

    #setup for c_w0 c_b0
    xwCoefs <- mlm1.z.sim@fixef
    #use names to select cols:
    selCols <- rm.na(match(dimnames(xwCoefs)[[2]],dimnames(data)[[2]]))
    fePred <- xwCoefs%*%rbind(1,t(data[,selCols]))
    randEffs <- mlm1.z.sim@ranef[[1]]
    id.re <- as.numeric(dimnames(randEffs)[[2]]) #effectively, we require numeric ids..
    idx <- match(newId,id.re)
    rePred <- randEffs[,idx,1] # asms random intercept model (dim: nsims X nrecs)
    sigs <- mlm1.z.sim@sigma
    #build synthetic data from model 1 as null:
    eps <- matrix(rnorm(nSims*nrecs,mean=0,sd=rep(sigs,each=nrecs)),byrow=T,nrow=nSims)
    newZs <- fePred + rePred + eps # synthetic data, all to get a tau.
    tmpdata.z <- data[,c(treat.varname,"idNew")] # 1st col is placeholder for newZs

    ncols.tmp <- dim(tmpdata.y)[[2]]
    tau0.w <- tau0.b <- tau0.ols <- sd.alpha.y.ucm <- sd.eps.y.ucm <- rep(NA,nSims)
    c_w0s <- c_b0s <- rep(NA,nSims)
    altNsims <- nSims
    if (!lines) altNsims <- 2 # partial speedup
    for (i in 1:nSims) {
        tmpdata.y[,1] <- newYs[i,] #substitute new outcomes
        #y~
        fmlaStr <- as.character(modelRun$Yfmla)
        yFmla <- as.formula(paste(fmlaStr[2],fmlaStr[1],attr(terms(modelRun$treatment),"term.labels"),"+",modelRun$treatment.gp,substring(fmlaStr[3],2),collapse=""))
        fit.y <- lmer(yFmla,data=tmpdata.y) # need kludge to get past y~. + ... formula.  Need to explicitly list z, z.mn instead
        taus <- lmeExtract(fit.y,attr(terms(modelRun$treatment),"term.labels"),modelRun$treatment.gp)$tau
        tau0.w[i] <- taus[1]
        tau0.b[i] <- taus[2]
        yFmla.ols <- as.formula(paste(attr(terms(modelRun$outcome),"term.labels"),paste(as.character(modelRun$treatment),collapse=""),sep=""))
        #add back group means
        tmpdata.y.ols <- tmpdata.y[,1:2]
        tmpdata.y.ols[,2] <- tmpdata.y.ols[,2] + tmpdata.y[,3]
        fit.ols <- lm(yFmla.ols,data=tmpdata.y.ols)
        ##browser()
        tau0.ols[i] <- lmExtract(fit.ols,attr(terms(modelRun$treatment),"term.labels"))$tau
        #z~
        ##different data needed here, and not GMCed.
        tmpdata.z[,1] <- newZs[i,] #substitute new outcomes
        fit.z <- lmer(modelRun$Zfmla,data=tmpdata.z)
        vcs <- lmeExtract.varcomp(fit.z)$vcomps
        c_w0s[i] <- vcs[2]
        c_b0s[i] <- vcs[1]
        #ucm
        fit.ucm.y <- lmer(modelRun$Yfmla,data=tmpdata.y[,c(outcome.name,"idNew")])
        ucm.vc <- lmeExtract.varcomp(fit.ucm.y)
        sd.alpha.y.ucm[i] <- sqrt(ucm.vc$vcomps[1])
        sd.eps.y.ucm[i] <- sqrt(ucm.vc$vcomps[2])
    }

    bd <-cbind(tau0.w-tau1.w,tau0.b-tau1.b,tau0.ols-tau1.ols)
    
    pObj <- extractParams(modelRun)
    fixedBdMat <- matrix(pObj$bias.diffs,ncol=3,nrow=nSims,byrow=T) #replacement with consts
    bd[,!allowVarBool[1:3]] <- fixedBdMat[,!allowVarBool[1:3]]
    
    #get gamm.z'V(w)gamm.z from model 1 (for Z, Y)
    ##****fix as.data.frame to work in call using dimnames, not names.  Fix in other call as well...
    xprods <- makeSpecProds(modelRun$Wfmla,as.data.frame(mlm1.z.sim@fixef),as.data.frame(mlm1.y.sim@fixef),modelRun$varW,modelRun$varX)
    varW <- xprods$varW #these are reordered to reflect proper variable ordering.
    varX <- xprods$varX
    gz.vw.gz <- xprods$gz.vw.gz
    bz.vx.bz <- xprods$bz.vx.bz
    gy.vw.gz <- xprods$gy.vw.gz
    gy.vw.gy <- xprods$gy.vw.gy
    #
    # check dimension -- will need to take diag when return is from a sim vector...
    #
    by.vx.bz <- xprods$by.vx.bz
    by.vx.by <- xprods$by.vx.by
    
    
    #need bnd.fnc to vary by sim, which means need new ucm model each run.
    ####*************
    
    #scale between difference by a factor (test theory on why points 'blow up':
    ###bd[,2] <- bd[,2]*betScaleFact
    sigArray <- cbind(c_w0s,c_w1s,c_b0s,c_b1s)
    fixedSigMat <- matrix(as.vector(pObj$sigs),nrow=nSims,ncol=4,byrow=T)
    sigArray[,!allowVarBool[4:7]] <- fixedSigMat[,!allowVarBool[4:7]]
    rslt <- array(NA,c(nSims,nTarg,4)) #there are 4 'output' conf parms
    if (debug) {
       pObj <- extractParams(modelRun)
       if (T) {
         curdev <- dev.cur()
         quartz()
         par(mfrow=c(3,4))
         plot(density(tau0.w));abline(v=pObj$tau.w[1],col=8)
         plot(density(tau1.w));abline(v=pObj$tau.w[2],col=8)
         plot(density(tau0.b));abline(v=pObj$tau.b[1],col=8)
         plot(density(tau1.b));abline(v=pObj$tau.b[2],col=8)
         plot(density(tau0.ols));abline(v=pObj$tau.ols[1],col=8)
         plot(density(tau1.ols));abline(v=pObj$tau.ols[2],col=8)
         plot(density(c_w0s));abline(v=pObj$sigs[1,1],col=8)
         plot(density(c_w1s));abline(v=pObj$sigs[2,1],col=8)
         plot(density(c_b0s));abline(v=pObj$sigs[1,2],col=8)
         plot(density(c_b1s));abline(v=pObj$sigs[2,2],col=8)
         quartz()
         pairs(cbind(tau0.w,tau1.w,tau0.b,tau1.b,tau0.ols,tau1.ols,c_w0s,c_w1s,c_b0s,c_b1s),pch='.')
         dev.set(curdev)
       }
       coefVar <- function(muTarg,muSamp) {sigma<-apply(muSamp,2,sd);sigma/muTarg}
       cvs <- coefVar(c(pObj$tau.w[1],pObj$tau.w[2],pObj$tau.b[1],pObj$tau.b[2],pObj$tau.ols[1],pObj$tau.ols[2],as.vector(pObj$sigs)),cbind(tau0.w,tau1.w,tau0.b,tau1.b,tau0.ols,tau1.ols,c_w0s,c_w1s,c_b0s,c_b1s))

       print(cvs[1:6])
       print(cvs[7:10])
    }

  
    #add to fakepObj
    
    for (i in 1:nSims) {
        sigs0 <- matrix(sigArray[i,],2,2,byrow=F)
        fakepObj <- pObj
        fakepObj$sigs=sigs0
        fakepObj$bias.diffs=bd[i,]
        fakepObj$gz.vw.gz= gz.vw.gz[i,i]
        fakepObj$bz.vx.bz=bz.vx.bz[i,i]   #implicitly taking diag()
        fakepObj$gy.vw.gz=gy.vw.gz[i,i]
        fakepObj$gy.vw.gy=gy.vw.gy[i,i]
        fakepObj$by.vw.bz=by.vx.bz[i,i]
        fakepObj$by.vx.by=by.vx.by[i,i]
        fakepObj$sds.y.ucm$sd.alpha.y.ucm <- sd.alpha.y.ucm[i]
        fakepObj$sds.y.ucm$sd.eps.y.ucm <-  sd.eps.y.ucm[i]
        
        bnd.f <- makeBnds(fakepObj,varyParm)
        tryrslt <- try(recovParms<-recover(fakepObj,varyParm=varyParm,bnd.f,nParms=nTarg,tau.max=tau.max,gpSize=gpSize),silent=F)
        if (!isErrClass(tryrslt)) rslt[i,,] <- recovParms[[1]]
    }
    #browser()
    
    #remove outer 1-conf % of curves (on zeta dimension):
    cutPt1 <- quantile((mins<-unlist(apply(rslt[,,1],1,min))),prob=(1-conf)/2)
    cutPt2 <- quantile((maxs<-unlist(apply(rslt[,,1],1,max))),prob=conf+(1-conf)/2)
    b <- mins>=cutPt1 & maxs<=cutPt2
    list(lines=rslt[b,,],by.vx.bz=diag(by.vx.bz))
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
    
    winBias0 <- (betaZ*varX*betaY+zetaYZ)/(betaZ^2*varX+c_W1)  #corrected 5Feb15; orig bZ(vX)bY - wrong.
    betBias0 <- (gammaZ*varW*gammaY+deltaYZ)/(gammaZ^2*varW+c_B1) #fixed gammaZ*varW*gammaY to proper denom as per text
    olsBias0 <- (betaZ*varX*betaY+gammaZ*varW*gammaY+zetaYZ+deltaYZ)/(betaZ^2*varX+gammaZ^2*varW+c_W1+c_B1) #same correction
    
    winBias1 <- zetaYZ/c_W1
    betBias1 <- deltaYZ/c_B1
    olsBias1 <- (zetaYZ+deltaYZ)/(c_W1+c_B1)
    sigs <- matrix(c(c_W0,c_W1,c_B0,c_B1),2,2,byrow=F)
    list(sigs=sigs,biasDiffs=c(winBias0-winBias1,betBias0-betBias1,olsBias0-olsBias1),winBias=c(winBias0,winBias1),betBias=c(betBias0,betBias1),olsBias=c(olsBias0,olsBias1), sd.y.ucm=c(sd.alpha.y.ucm,sd.eps.y.ucm))
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


if (F) {
 vcs <- lmeExtract.varcomp(mlm1.z)
 cat("vcomps=",vcs$vcomps,"\n")
 omega <- vcs$vcomps[1]*outer(idNew,idNew,"==")
 diag(omega) <- diag(omega)+vcs$vcomps[2]
 omegaInv <- solve(omega)
 
 ZZpInvZp <- Znew[,1]%*%omegaInv/as.numeric(Znew[,1]%*%omegaInv%*%Znew[,1])
 ZZpInvZp.mn <- Znew[,2]%*%omegaInv/as.numeric(Znew[,2]%*%omegaInv%*%Znew[,2])  #assumes Z-Zbar, Zbar in design...
 ZZpInvZp.ols <- as.numeric(Z/as.numeric(t(Z)%*%Z))
 
 #ZZpInvZp=ZZpInvZp,ZZpInvZp.mn=ZZpInvZp.mn,ZZpInvZp.ols=ZZpInvZp.ols
 
 mlm0.z.sim<-sim(modelRun$mlm0.z,nSims)
 c_b0s <- apply(mlm0.z.sim@ranef[[1]][,,1],1,var)
 c_w0s <- mlm0.z.sim@sigma^2
 
 mlm0.y.sim<-sim(modelRun$mlm0.y,nSims)
 tau0.w <- mlm0.y.sim@fixef[,treat.varname]
 tau0.b <- mlm0.y.sim@fixef[,treat.gpname]

 lm0.y.sim<-sim(modelRun$ols0.y,nSims)
 tau0.ols <- lm0.y.sim@coef[,treat.varname]
 
 #fix biases too:
 cat(dim(modelRun$ZZpInvZp),'\n')
 cat(dim(modelRun$ZZpInvZp.mn),'\n')
 cat(dim(wCoefs),'\n')
 mm.w<-model.matrix(modelRun$Xfmla,modelRun$newData)[,-1]
 mm.b<-model.matrix(modelRun$Wfmla,modelRun$newData)[,-1]
 tau0.w <- tau1.w + as.numeric(modelRun$ZZpInvZp%*%mm.w%*%t(xCoefs))
 tau0.b <- tau1.b + as.numeric(modelRun$ZZpInvZp.mn%*%mm.b%*%t(wCoefs))
 tau0.ols <- tau1.ols + as.numeric(modelRun$ZZpInvZp.ols%*%cbind(mm.w,mm.b)%*%t(cbind(xCoefs,wCoefs)))
 
}

