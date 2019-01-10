# # Unused packages: multcomp arm foreign reshape
#
#
# #---------FROM model_fits---------
# rescaleSetVal <- function(extParms) {
#   #set values for rescaling based on fitted models
#   c(sqrt(extParms$varY),sqrt(extParms$varZ))
# }
#
# #---------FROM plotting--------
# lwDelete <- function(fmla,data) {
#   opt0 <- options()
#   options(na.action="na.pass")
#   mm <- model.matrix(fmla,data)
#   bb<-apply(is.na(mm),1,sum)==0
#   options(opt0)
#   data[bb,]
# }
# idealizedParms <- function(varX,varW,tau,betaY,betaZ,gammaY,gammaZ,zetaY,zetaZ,deltaY,deltaZ,varAlphaY,varAlphaZ,varEpsilonY,varEpsilonZ) {
#
#   #scalar x,w only, for now.
#   #Betw. Group Size -> infty for now
#   c_W0 <- betaZ^2*varX + zetaZ^2 + varEpsilonZ
#   c_B0 <- gammaZ^2*varW + deltaZ^2 + varAlphaZ
#   c_W1 <- c_W0 - betaZ^2*varX
#   c_B1 <- c_B0 - gammaZ^2*varW
#
#   covWinYZ <- betaZ*varX*betaY + zetaZ*zetaY
#   covBetYZ <- gammaZ*varW*gammaY + deltaZ*deltaY
#
#   sd.eps.y.ucm <- sqrt(tau^2*(c_W0) + betaY^2*varX + zetaY^2 + varEpsilonY + 2*tau*covWinYZ)
#   sd.alpha.y.ucm <- sqrt(tau^2*(c_B0) + gammaY^2*varW + deltaY^2 + varAlphaZ + 2*tau*covBetYZ)
#
#   zetaYZ <- zetaY*zetaZ
#   deltaYZ <- deltaY*deltaZ
#
#   ##replace w/ function calls?
#
#   winBias0 <- (betaZ*varX*betaY+zetaYZ)/(betaZ^2*varX+c_W1)  #corrected 5Feb15; orig bZ(vX)bY - wrong.
#   betBias0 <- (gammaZ*varW*gammaY+deltaYZ)/(gammaZ^2*varW+c_B1) #fixed gammaZ*varW*gammaY to proper denom as per text  -- group size -> Infty
#   olsBias0 <- (betaZ*varX*betaY+gammaZ*varW*gammaY+zetaYZ+deltaYZ)/(betaZ^2*varX+gammaZ^2*varW+c_W1+c_B1) #same correction
#
#   winBias1 <- zetaYZ/c_W1
#   betBias1 <- deltaYZ/c_B1
#   olsBias1 <- (zetaYZ+deltaYZ)/(c_W1+c_B1)
#   sigs <- matrix(c(c_W0,c_W1,c_B0,c_B1),2,2,byrow=F)
#   list(sigs=sigs,biasDiffs=c(winBias0-winBias1,betBias0-betBias1,olsBias0-olsBias1),winBias=c(winBias0,winBias1),betBias=c(betBias0,betBias1),olsBias=c(olsBias0,olsBias1), sd.y.ucm=c(sd.alpha.y.ucm,sd.eps.y.ucm))
# }
#
#
#
