# family backgroud (# books in home) on reading comprehension score #

# Y: reading comprehension test raw score (608)
# G: different schools (3)
# Z: number of books in home (178)
# X: sex (165),hours of homework (167), word knowledge raw score (247)
# W: Z.groupmean,X.groupmean,size of the school (25,26),type of community served (18)

# depends:
#  biasAmpR_v2r3.R
#
# data source:
#  IEA.RData (preprocessed dataset)

library(foreign)
library(dplyr)
library(plyr)
library(reshape)
library(nlme)
library(lme4)
library(arm)
library(lattice)

#read support fns
if (!exists('runModels')) source("biasAmpR_v2r3.R")
#read data
if (!exists('scotland')) load("ObsStudies_IEA.RData")

# run lmer and lm models

#choice of country names
cnames <-c("scotland","sweden","italy")


set.seed(12341)
for (ii in 1:length(cnames)) {
    
  country <- get(cnames[ii])

  country$pop <- log(country$pop)
  country$pop[is.infinite(country$pop)] <- NA
  country$type <- factor(country$type)
  #remove nas:
  cdat <- model.frame(score~num_books + sex + word_knowl + homework + type + pop + school,data=country) #should remove NAs

  #fix type factors:
  cdat$type0 <- 1*(cdat$type==0)
  cdat$type1 <- 1*(cdat$type==1)
  cdat$type2 <- 1*(cdat$type==2)
  cdat$type3 <- 1*(cdat$type==3)
  cdat$type4 <- 1*(cdat$type==4)
  cdat$type5 <- 1*(cdat$type==5)
  cdat$type6 <- 1*(cdat$type==6)

  nPerSchool <- mean(table(cdat$school))
  nPerSchool <- prod(table(cdat$school))^(1/length(unique(cdat$school))) # works better with geo-mean

  #indiv preds
  fmlaX <- ~sex + word_knowl + homework
  #group var preds

  fmlaW.str <- "~pop"
  first <- TRUE # flag
  for (i in 0:6) {# add only types with variation
      tstr <- paste("type",as.character(i),collapse="",sep="")
      len <- length(table(cdat[,tstr]))
      if (len>1) {
           if (!first) fmlaW.str <- paste(fmlaW.str,tstr,collapse="",sep="+")
           first <- FALSE # ensures ref category
      }
  }
  fmlaW <- as.formula(fmlaW.str)


  #set treatment
  fmlaZ <- ~num_books
  fmlaY <- ~score

  #fit models
  mdl.fit <- runModels(outcome=fmlaY, treatment=fmlaZ, level1.pred = fmlaX, level2.pred = fmlaW, group = ~school, data=cdat)
  #extract variance comps & some bias diffs
  pObj <- extractParams(mdl.fit)
  #for ObsStudies paper:
  print(cnames[ii])
  print(summary(mdl.fit$mlm1.y))
  print(summary(mdl.fit$ols1.y))
  print(paste("Table 1 for: ",cnames[ii]))
  rsltTab1 <- rbind(cbind(pObj$tau.w,pObj$tau.b,pObj$tau.ols),pObj$bias.diffs)
  dimnames(rsltTab1) <- list(c("tau0","tau1","diff"),c("Within","Between","OLS"))
  print(round(rsltTab1,2))
  print(round(as.vector(t(cbind(pObj$sigs,pObj$sigs[,2]/apply(pObj$sigs,1,sum)))),2))
  print(round(cbind(pObj$sigs,pObj$sigs[,2]/apply(pObj$sigs,1,sum)),2))
  print(paste("gy,vy,t.gz ",paste(round(sqrt(c(pObj$gy.vw.gy,pObj$sds.y.ucm$sd.alpha.y.ucm^2,2^2*pObj$gz.vw.gz)),3),collapse=", "),sep=" "))

  parms <- c("zeta","delta","beta","gamma")
  paramIdx <- 4
  paramCh <- parms[paramIdx]
  rescaleParms <- c(1,1)
  bnd.f <- makeBnds(pObj,param=paramCh)

  recovParms4<-recovParms<-recover(pObj,varyParm=paramCh,bnd.f,nParm=201,tau.max=1,gpSize=nPerSchool)

  cat("cond num=",recovParms$condNum,"\n")

  #FIND THE RIGHT 4th Param from the 3rd param beta-based value for eta...
  recovParms3<-recover(pObj,varyParm="beta",makeBnds(pObj,param="beta"),nParm=201,tau.max=1,gpSize=nPerSchool)
  
  gamma.at.beta <- recovParms3$zetaDeltaMat[which.min(abs(recovParms3$zetaDeltaMat[,3]-pObj$by.vx.bz)),4]
  
  #get other 'target' values:
  nu.at.zeta.0 <- recovParms$zetaDeltaMat[which.min(abs(recovParms$zetaDeltaMat[,1])),paramIdx]
  nu.at.delta.0 <- recovParms$zetaDeltaMat[which.min(abs(recovParms$zetaDeltaMat[,2])),paramIdx]
  
  #reduce range if info from Beta warrants it
  nu.g.min <- min(recovParms3$zetaDeltaMat[,4])
  nu.g.max <- max(recovParms3$zetaDeltaMat[,4])
  #fix range here
  b.in.range <- recovParms4$zetaDeltaMat[,4] >= nu.g.min & recovParms4$zetaDeltaMat[,4] <= nu.g.max
  recovParms$zetaDeltaMat <- recovParms$zetaDeltaMat[b.in.range,]
  recovParms$parmRange <- recovParms$parmRange[b.in.range]
  recovParms$plausible.taus <- recovParms$plausible.taus[b.in.range]


  #plot params set for .png inlcuded in LaTeX file; adjust as nec.
  lcex <- 3
  pcex <- 2
  png(paste(cnames[ii],"png",sep=".",collapse=""),width=pcex*480,height=pcex*480)
  tvals <- c(gamma.at.beta,pObj$gy.vw.gz,nu.at.zeta.0,nu.at.delta.0)
  tpch <- c(0,4,1,3)
  taus <- list(ols=pObj$tau.ols[2],win=pObj$tau.w[2])
  confPts <- NULL #default
  if (!is.null(mdl.fit$bsSamps)) {
      if (dim(mdl.fit$bsSamps)[1]!=0) {
         rslt0 <- bsConfParms(mdl.fit$bsSamps,201,bnd.f,nPerSchool,tau.max=1,trim=0.025)$parms
         confPts <- cbind(as.vector(rslt0[,,2]),as.vector(rslt0[,,1]))
      }
  }
  
  rslt0 <- bsConfParmsSim(mdl.fit,nSims=250,nTarg=201,bnd.f=bnd.f,gpSize=nPerSchool,tau.max=1,debug=F,allowVarBool=c(T,T,T,rep(T,4)),line=F)
  confPts <- cbind(as.vector(rslt0$lines[,,2]),as.vector(rslt0$lines[,,1]))
  
  gamma.at.beta025 <- recovParms3$zetaDeltaMat[which.min(abs(recovParms3$zetaDeltaMat[,3]-quantile(rslt0$by.vx.bz,probs=c(0.025)))),4]
  gamma.at.beta975 <- recovParms3$zetaDeltaMat[which.min(abs(recovParms3$zetaDeltaMat[,3]-quantile(rslt0$by.vx.bz,probs=c(0.975)))),4]

  print("target for nu(by.vx.bz): ")
  print(nu.at.zeta.0)
  print("range for nu(by.vx.bz): ")
  print(c(gamma.at.beta025, gamma.at.beta,gamma.at.beta975),digits=3)
  
  plot(zdPlot(recovParms$zetaDeltaMat[,1],recovParms$zetaDeltaMat[,2],recovParms$parmRange,rescaleParms=rescaleParms,confPts=NULL,confPtsCol="darkslategrey",targetVals=tvals,targetPch=tpch,taus=taus,cW=pObj$sigs[2,1],cB=pObj$sigs[2,2],cex=lcex))
  dev.off()
}
