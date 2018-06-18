require(lattice)

#used to generate figure 1 in the Observational Studies paper

#set up a grid
N <- 201
lbd <- -4; ubd <- 4
#set grid values
zeta <- seq(lbd,ubd,length=N)
delta <- seq(lbd,ubd,length=N)

#functions that compute absolute bias differences for different estimators
betMwin <- function(zeta,delta,cB=1,cW=1) {abs(delta/cB)-abs(zeta/cW)}
olsMwin <- function(zeta,delta,cB=1,cW=1) {abs((zeta+delta)/(cW+cB))-abs(zeta/cW)}
glsMwin <- function(zeta,delta,cB=1,cW=1,lambda=.5) {abs((zeta+lambda*delta)/(cW+lambda*cB))-abs(zeta/cW)}

winBias <- function(zeta,delta,cB=1,cW=1) {zeta/cW}
betBias <- function(zeta,delta,cB=1,cW=1) {delta/cB}
olsBias <- function(zeta,delta,cB=1,cW=1) {(zeta+delta)/(cW+cB)}
glsBias <- function(zeta,delta,cB=1,cW=1,lambda=.5) {(zeta+lambda*delta)/(cW+lambda*cB)}

#choose the comparison here: betMwin or olsMwin or glsMwin
#one can even plot abs bias for a single estimator, of course, as these are also defined above.
compFUN <- "olsMwin"

bd.1.1 <-outer(zeta,delta,FUN=compFUN,cB=2,cW=2)
bd.1.3 <-outer(zeta,delta,FUN=compFUN,cB=1,cW=3)
bd.3.1 <-outer(zeta,delta,FUN=compFUN,cB=3,cW=1)

#helper for lattice levelplot: x & y axis lines only
pfunct = function(x, y, z,...) {
  panel.levelplot(x, y, z,...)
  panel.abline(h=0,col.line=8)
  panel.abline(v=0,col.line=8)
}

#same helper, with y=x line added.
pfunct11 = function(x, y, z,...) {
    panel.levelplot(x, y, z,...)
    panel.abline(h=0,col.line=8)
    panel.abline(v=0,col.line=8)
    panel.abline(a=0,b=1,col.line=1)
}

#set up for legend
at.pts <- seq(-4,4,length=9)
spec.lbl <- as.character(at.pts)
spec.lbl[1] <- "OLS\nbetter"
spec.lbl[length(spec.lbl)] <- "Within\nbetter"

lcex <- 2.5
png("dz_1_1.png",width=480*2,height=480*2)
plot(levelplot(t(bd.1.1), at=seq(-4,4,length=17),colorkey=list(at=at.pts,labels=list(labels=spec.lbl,cex=lcex)),scales=list(x=list(cex=lcex),y=list(cex=lcex)),panel=pfunct11,row.values=zeta,column.values=delta,xlab=list(expression(delta^{yz}),cex=lcex),ylab=list(expression(zeta^{yz}),cex=lcex)))
dev.off()
png("dz_1_3.png",width=480*2,height=480*2)
plot(levelplot(t(bd.1.3), at=seq(-4,4,length=17),colorkey=list(at=at.pts,labels=list(labels=spec.lbl,cex=lcex)),scales=list(x=list(cex=lcex),y=list(cex=lcex)),panel=pfunct,row.values=zeta,column.values=delta,xlab=list(expression(delta^{yz}),cex=lcex),ylab=list(expression(zeta^{yz}),cex=lcex)))
dev.off()
png("dz_3_1.png",width=480*2,height=480*2)
plot(levelplot(t(bd.3.1), at=seq(-4,4,length=17),colorkey=list(at=at.pts,labels=list(labels=spec.lbl,cex=lcex)),scales=list(x=list(cex=lcex),y=list(cex=lcex)),panel=pfunct,row.values=zeta,column.values=delta,xlab=list(expression(delta^{yz}),cex=lcex),ylab=list(expression(zeta^{yz}),cex=lcex)))
dev.off()

