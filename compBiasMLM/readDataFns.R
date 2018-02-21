
makeCdat <- function(cname,stdz=F) {

#hard coded data info embedded here.
  country <- get(cname)
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
  #
  nPerSchool <- mean(table(cdat$school))
  #  nPerSchool <- prod(table(cdat$school))^(1/length(unique(cdat$school))) # geo-mean
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
  
  #treatment
  fmlaZ <- ~num_books
  fmlaY <- ~score
  #
  return(list(cdat=cdat,fmlaY=fmlaY,fmlaZ=fmlaZ,fmlaX=fmlaX,fmlaW=fmlaW,group=formula(~school),nPer=nPerSchool))
}



mlmSens <- function(outcome=fmlaY, treatment=fmlaZ, level1.pred = fmlaX, level2.pred = fmlaW, group = ~school, nPer=Inf, data=cdat) {
   
   mdl.fit <- runModels(outcome,treatment,level1.pred,level2.pred,group,data)  #by default does the gmc.
    #extract variance comps & some bias diffs
    pObj <- extractParams(mdl.fit)
    #for paper:
    print(summary(mdl.fit$mlm1.y))
    print(summary(mdl.fit$ols1.y))
    print(round(rbind(cbind(pObj$tau.w,pObj$tau.b,pObj$tau.ols),pObj$bias.diffs),2))
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
    
    ##print(nu.at.zeta.0,digits=2)
    
    #reduce range if info from Beta warrants it
    nu.g.min <- min(recovParms3$zetaDeltaMat[,4])
    nu.g.max <- max(recovParms3$zetaDeltaMat[,4])
    #fix range here
    b.in.range <- recovParms4$zetaDeltaMat[,4] >= nu.g.min & recovParms4$zetaDeltaMat[,4] <= nu.g.max
    recovParms$zetaDeltaMat <- recovParms$zetaDeltaMat[b.in.range,]
    recovParms$parmRange <- recovParms$parmRange[b.in.range]
    return(list(recovParms=recovParms,pObj=pObj,mdl.fit=mdl.fit))
}