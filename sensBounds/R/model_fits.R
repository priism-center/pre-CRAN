# Add documentation



#---------HELPER FUNCTIONS------

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


#helper function to get s.e. of param ests.
seExtract <- function(lmerObj) coef(summary(lmerObj))[,"Std. Error"]

scaleAll <- function(x) {

  # try to scale factors whenever possible

  if (is.numeric(x)) return(scale(x))
  if (is.factor(x)) return(as.numeric(levels(x)[x]))
  return(x) #else it's probably character data
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
  x <- as.array(model.matrix(level1.pred,data)[,-1, drop=F]) # drop const
  y <- data[,outcomeName]
  z.ctrd <- data[,treatName]
  z.mn.name <- paste0(treatName,".mn")
  dataNames <- names(data)
  if (!is.na(match(z.mn.name,dataNames))) z.mn <- data[,z.mn.name] else z.mn <- 0 #should work for the special case of not using runModels to build dat0
  w <- as.array(model.matrix(level2.pred,data)[,-1, drop=F]) #drop const.
  p <- dim(x)[2]
  q <- dim(w)[2]
  id <- renumber(data[,groupName])
  M <- length(unique(id))
  gpSize <- geoMeanGroupSize(id)
  #fix z; current data version is centered, and needs to not be, going forward
  list(N=N,M=M,x=x,y=y,z=z.ctrd+z.mn,w=w,p=p,q=q,id=id,gpSize=gpSize)
}


combine.formula <- function(fmla1,fmla2) {
  #works on RHS formulas, e.g., ~x1+x2...
  terms1 <- attr(terms.formula(fmla1),"term.labels")
  terms2 <- attr(terms.formula(fmla2),"term.labels")
  as.formula(paste0("~",paste(c(terms1,terms2),collapse="+")))
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


#------- Where is this used? ---------
rescaleSetVal <- function(extParms) {
  #set values for rescaling based on fitted models
  c(sqrt(extParms$varY),sqrt(extParms$varZ))
}






#--------PRINT RESULTS---------
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





