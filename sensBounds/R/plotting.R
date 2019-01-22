# Add documentation
#' @importFrom plyr .


#---------HELPER FUNCTIONS-----------

makeSpecProds <- function(level2.pred, fe_mlm1.z, fe_mlm1.y, varW, varX) {
    # helper function to extractParams to generate terms needed eventually for bounds get gamm.z'V(w)gamm.z
    # from model 1 (for Z) add sd*se to gamma or beta [ref]

    wNames <- attr(terms(level2.pred), "term.labels")
    nmsZ <- names(fe_mlm1.z)
    if (is.null(nmsZ))
        nmsZ <- dimnames(fe_mlm1.z)[[2]]  #handles 2 cases
    idxW <- match(wNames, nmsZ)
    xNames <- names(fe_mlm1.z)[-idxW][-1]  #drop 'intercept'
    idxX <- match(xNames, nmsZ)
    if (class(fe_mlm1.z) == "data.frame") {
        gamm.z <- as.matrix(fe_mlm1.z[, idxW])
        beta.z <- as.matrix(fe_mlm1.z[, idxX])
    } else {
        # assume single row
        gamm.z <- t(as.matrix(fe_mlm1.z[idxW, drop = F]))
        beta.z <- t(as.matrix(fe_mlm1.z[idxX, drop = F]))
    }

    # next bit is based on biased parms, but worth doing for comparison purposes (and bounds):

    nmsY <- names(fe_mlm1.y)
    if (is.null(nmsY))
        nmsY <- dimnames(fe_mlm1.y)[[2]]  #handles 2 cases
    idxW <- match(wNames, nmsY)
    idxX <- match(xNames, nmsY)  #need to exclude treatment from these.
    if (class(fe_mlm1.y) == "data.frame") {
        gamm.y <- as.matrix(fe_mlm1.y[, idxW])
        beta.y <- as.matrix(fe_mlm1.y[, idxX])
    } else {
        # single row
        gamm.y <- t(as.matrix(fe_mlm1.y[idxW, drop = F]))
        beta.y <- t(as.matrix(fe_mlm1.y[idxX, drop = F]))
    }
    # fix cols of var mat to align with gamma & beta
    idxW <- match(wNames, dimnames(varW)[[1]])
    idxX <- match(xNames, dimnames(varX)[[1]])
    varW <- as.matrix(varW[idxW, idxW])
    varX <- as.matrix(varX[idxX, idxX])

    bndProdList <- list(gz.vw.gz = gamm.z %*% varW %*% t(gamm.z), bz.vx.bz = beta.z %*% varX %*% t(beta.z),
        gy.vw.gz = gamm.y %*% varW %*% t(gamm.z), gy.vw.gy = gamm.y %*% varW %*% t(gamm.y), by.vx.bz = beta.y %*%
            varX %*% t(beta.z), by.vx.by = beta.y %*% varX %*% t(beta.y))

    list(varX = varX, varW = varW, bndProdList = bndProdList)

}

extractParams <- function(runModelRslt) {

    # extract certain ests of var comps & bias diffs & taus (extraneous) for use in sensitivity plots (and
    # diagnostics)

    data <- runModelRslt$data  #some vars are added to this data.frame, so best to use it instead of original data
    level2.pred <- runModelRslt$level2.pred
    treatment <- runModelRslt$treatment
    outcome <- runModelRslt$outcome
    # get gamm.z'V(w)gamm.z from model 1 (for Z)
    xprods <- makeSpecProds(level2.pred, lme4::fixef(runModelRslt$mlm1.z), lme4::fixef(runModelRslt$mlm1.y),
        runModelRslt$varW, runModelRslt$varX)
    varW <- xprods$varW  #these are reordered to reflect proper variable ordering.
    varX <- xprods$varX

    outcome.varname <- attr(terms(outcome), "term.labels")
    varY <- var(data[, outcome.varname])
    # note: treatment.gp is already extracted from the formula form... could be changed to be consistent w/
    # other passing protocol.
    treatment.varname <- attr(terms(treatment), "term.labels")
    treatment.gpname <- paste0(treatment.varname, ".mn")
    varZ <- var(data[, treatment.varname])

    # WITHIN: there need to be two extract pgms.  ONE does fixefs and se, other does VC and lvcv... BETWEEN:
    # most effects are extracted.  just need between est of tau (more elegant would be to reuse lmeExtract, but
    # would be redundant and add comp time):

    # add s.e. as well
    extr0 <- lmeExtract(runModelRslt$mlm0.y, treatment.varname, treatment.gpname)
    tau0.w <- extr0$tau[1]
    tau0.b <- extr0$tau[2]
    #
    extr1 <- lmeExtract(runModelRslt$mlm1.y, treatment.varname, treatment.gpname)
    tau1.w <- extr1$tau[1]
    tau1.b <- extr1$tau[2]

    within.bias.diff <- tau0.w - tau1.w  # based on lme, fixef
    between.bias.diff <- tau0.b - tau1.b  #based on lme

    # OLS more complicated to find correct coef to extract.  need to find the `I(treat+treat.mn)` entry.
    treatConstrName <- attr(terms(runModelRslt$ols0.y), "term.labels")  #there is only one predictor in this model & it's the treatment.

    extr <- lmExtract(runModelRslt$ols0.y, treatConstrName)
    tau0.ols <- extr$tau
    extr <- lmExtract(runModelRslt$ols1.y, treatConstrName)
    tau1.ols <- extr$tau
    #
    ols.bias.diff <- tau0.ols - tau1.ols

    ##### varcorr:
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

    sigs <- matrix(c(sigma.sq.within.0, sigma.sq.within.1, sigma.sq.between.0, sigma.sq.between.1), 2, 2, byrow = F)

    # these are used for bounds, so sampling variability of these is less crucial:
    ucm.vc <- lmeExtract.varcomp(runModelRslt$mlm.ucm.y)
    sd.alpha.y.ucm <- sqrt(ucm.vc$vcomps[1])
    sd.eps.y.ucm <- sqrt(ucm.vc$vcomps[2])

    list(sds.y.ucm = list(sd.eps.y.ucm = sd.eps.y.ucm, sd.alpha.y.ucm = sd.alpha.y.ucm), sigma.sq.between.y.0 = sigma.sq.between.y.0,
        sigs = sigs, bias.diffs = c(within.bias.diff, between.bias.diff, ols.bias.diff), tau.w = c(tau0.w,
            tau1.w), tau.b = c(tau0.b, tau1.b), tau.ols = c(tau0.ols, tau1.ols), bndProdList = xprods$bndProdList,
        varY = varY, varZ = varZ)
}

makeBnds <- function(paramObj, param = "gamma") {

    # extract from paramObj - yields a function of tau that sets a bound on feasible eta; param allows one to
    # switch the 'open' param (see paper) and this yields a different set of bounds.

    # for convenience, extract these here:

    gz.vw.gz <- paramObj$bndProdList$gz.vw.gz
    gy.vw.gy <- paramObj$bndProdList$gy.vw.gy
    c.b.y <- paramObj$sds.y.ucm$sd.alpha.y.ucm
    bz.vx.bz <- paramObj$bndProdList$bz.vx.bz
    by.vx.by <- paramObj$bndProdList$by.vx.by
    c.w.y <- paramObj$sds.y.ucm$sd.eps.y.ucm

    if (param == "gamma") {
        f <- function(tau.max, bound = T) {
            if (bound) {
                return(as.numeric(c(sqrt(gz.vw.gz) * max(abs(tau.max) * sqrt(gz.vw.gz), sqrt(gy.vw.gy), c.b.y))))
            } else {
                return(as.numeric(sqrt(gy.vw.gy * gz.vw.gz)))
            }
        }
    } else {
        f <- function(tau.max, bound = T) {
            if (bound) {
                return(as.numeric(c(sqrt(bz.vx.bz) * max(abs(tau.max) * sqrt(bz.vx.bz), sqrt(by.vx.by), c.w.y))))
            } else {
                return(as.numeric(sqrt(by.vx.by * bz.vx.bz)))
            }
        }
    }
    return(f)
}

makeLinEqMat <- function(c_w0, c_w1, c_b0, c_b1, idx = 4, gpSize = Inf) {

    # set up the matrix of the 4 linear equations.  See paper.
    m <- matrix(0, 4, 4)
    m[1, 1] <- 1/c_w0 - 1/c_w1
    m[1, 3] <- 1/c_w0

    if (is.finite(gpSize)) {
        m[2, 1] <- (1/gpSize)/(c_w0/gpSize + c_b0) - (1/gpSize)/(c_w1/gpSize + c_b1)
        m[2, 2] <- 1/(c_w0/gpSize + c_b0) - 1/(c_w1/gpSize + c_b1)
        m[2, 3] <- (1/gpSize)/(c_w0/gpSize + c_b0)
        m[2, 4] <- 1/(c_w0/gpSize + c_b0)
    } else {
        m[2, 2] <- 1/c_b0 - 1/c_b1
        m[2, 4] <- 1/c_b0
    }
    m[3, 1] <- m[3, 2] <- 1/(c_w0 + c_b0) - 1/(c_w1 + c_b1)
    m[3, 3] <- m[3, 4] <- 1/(c_w0 + c_b0)
    m[4, idx] <- 1
    m
}

recover.n <- function(paramObj, varyParm = "gamma", bnd.f, nParms = 201, tau.max = 2, gpSize = Inf) {

    # solves a system of equations for 3 out of 4 of (zeta_1, delta_1,byxbz,gywgz) in terms of the 4th, based
    # on values for sigs and differences in bias for 3 estimation methods sigs: 2x2 matrix rows are (c.w0,c.b0)
    # (c.w1,c.b1) deltaBiases: bias(w,b,ols) or c(tau0.w-tau1.w,tau0.b-tau1.b,tau0.ols-tau1.ols)

    # one param must be 'open' - zeta_1, delta_1,etc.  varyParm choose it by naming it from parms IMPT: the
    # BOUNDS are only derived for gamma (gywgz) as the omitted param.

    # extract from paramObj:
    sigs <- paramObj$sigs
    deltaBiases <- paramObj$bias.diffs

    parms <- c("zeta", "delta", "beta", "gamma")

    idx <- match(varyParm, parms)
    m <- makeLinEqMat(sigs[1, 1], sigs[2, 1], sigs[1, 2], sigs[2, 2], idx, gpSize)

    # generate parmRange
    scaleFactor <- seq(-1, 1, length = nParms)
    parmRange <- scaleFactor * as.numeric(bnd.f(tau.max))

    # iterate through the 4th (omitted) param
    len <- length(parmRange)
    r <- matrix(NA, len, 4)  # result
    delts <- c(deltaBiases, 0)  #placeholder for 3 deltas and one open param
    # error trapping
    if (class(try(condNum <- kappa(m, exact = T))) == "try-error")
        print("unable to compute condition number\n")
    mInv <- solve(m)
    for (i in 1:len) {
        delts[4] <- parmRange[i]  #assume we know the omitted param.
        r[i, ] <- mInv %*% delts  #solve for 3 other params, given components and omitted param.
    }
    dimnames(r)[[2]] <- parms
    # zetaDeltaMat are the zetas and deltas consistent with the data; parmRange gives the open param for the
    # given row of zds ALT cond num: m <- makeLinEqMat(sigs[1,1],sigs[2,1],sigs[1,2],sigs[2,2],idx=3,gpSize) if
    # (class(try(condNum<-kappa(m,exact=T)))=='try-error') print('unable to compute condition number\n')
    list(zetaDeltaMat = r, plausible.taus = scaleFactor * tau.max, parmRange = parmRange, condNum = condNum)
}






# helper functs (mostly for zdPlot; useful in other contexts

# functions that compute asymptotic absolute bias differences for different estimators
betMwin <- function(zeta, delta, cB = 1, cW = 1) {
    abs(delta/cB) - abs(zeta/cW)
}
olsMwin <- function(zeta, delta, cB = 1, cW = 1) {
    abs((zeta + delta)/(cW + cB)) - abs(zeta/cW)
}
glsMwin <- function(zeta, delta, cB = 1, cW = 1, lambda = 0.5) {
    abs((zeta + lambda * delta)/(cW + lambda * cB)) - abs(zeta/cW)
}

winBias <- function(zeta, delta, cB = 1, cW = 1) {
    zeta/cW
}
betBias <- function(zeta, delta, cB = 1, cW = 1, gpSize = Inf) {
    delta/(cW/gpSize + cB)  #relies on 1/Inf == 0
}
olsBias <- function(zeta, delta, cB = 1, cW = 1) {
    (zeta + delta)/(cW + cB)
}
glsBias <- function(zeta, delta, cB = 1, cW = 1, lambda = 0.5) {
    (zeta + lambda * delta)/(cW + lambda * cB)
}

correctedTau.o <- function(tau.o, zeta, delta, cW, cB) {
    return(tau.o - olsBias(zeta, delta, cW, cB))
}



prePlotParams <- function(mdlFit, nGridPoints = 201, tau.max = 1, gpSize = Inf) {
    # Bounds calcs params needed for plots cond numb
    # nGridPoints: number of points to build the confounder line.
    # tau.max: abs(tau) used in the bounds calculation (plausible value supplied by user)
    # gpSize: for balanced designs, the number of subjects/group.  For others, try geom. mean or Inf.

    # FIND THE RIGHT open Param from the gamma- or beta-based value for eta.
    pObj <- extractParams(mdlFit)
    # corresponds to the bound when gamma is the open param:
    recovParmsGammaBased <- recover.n(pObj, varyParm = "gamma", makeBnds(pObj, param = "gamma"), nParm = nGridPoints,
        tau.max = tau.max, gpSize = gpSize)
    # beta-based equiv.:
    recovParmsBetaBased <- recover.n(pObj, varyParm = "beta", makeBnds(pObj, param = "beta"), nParm = nGridPoints,
        tau.max = tau.max, gpSize = gpSize)
    # get the gamma corresp. to this beta by searching through the list of evaluated points and finding closest
    # corresponding. this is one of the target values to report
    gammaYZ.at.beta <- recovParmsBetaBased$zetaDeltaMat[which.min(abs(recovParmsBetaBased$zetaDeltaMat[, "beta"] -
        as.numeric(pObj$bndProdList$by.vx.bz))), "gamma"]
    # get other 'target' values:
    gammaYZ.at.zeta.0 <- recovParmsGammaBased$zetaDeltaMat[which.min(abs(recovParmsGammaBased$zetaDeltaMat[,
        "zeta"])), "gamma"]
    gammaYZ.at.delta.0 <- recovParmsGammaBased$zetaDeltaMat[which.min(abs(recovParmsGammaBased$zetaDeltaMat[,
        "delta"])), "gamma"]
    gammaYZ.CS <- pObj$bndProdList$gy.vw.gz  #bound derived from C-S Ineq.

    # reduce range if info from Beta warrants it
    nu.g.min <- min(recovParmsBetaBased$zetaDeltaMat[, "gamma"])
    nu.g.max <- max(recovParmsBetaBased$zetaDeltaMat[, "gamma"])
    # adjust range here
    b.in.range <- recovParmsGammaBased$zetaDeltaMat[, "gamma"] >= nu.g.min & recovParmsGammaBased$zetaDeltaMat[,
        "gamma"] <= nu.g.max
    # by this point, it might be better not to call it 'gamma based'
    recovParmsGammaBased$zetaDeltaMat <- recovParmsGammaBased$zetaDeltaMat[b.in.range, ]
    recovParmsGammaBased$parmRange <- recovParmsGammaBased$parmRange[b.in.range]
    recovParmsGammaBased$plausible.taus <- recovParmsGammaBased$plausible.taus[b.in.range]

    bndVals <- c(gammaYZ.at.beta, gammaYZ.CS, gammaYZ.at.zeta.0, gammaYZ.at.delta.0)
    names(bndVals) <- c("gYZbeta", "gYZcs", "gYZzeta", "gYZdelta")

    list(zetaDeltaMat = recovParmsGammaBased$zetaDeltaMat, parmRange = recovParmsGammaBased$parmRange, plausible.taus = recovParmsGammaBased$plausible.taus,
        condNum = recovParmsGammaBased$condNum, bndVals = bndVals)
}





#---------MAIN FUNCTION---------

## ACTION: change size of labels in this plot... (is CEX= enough?)
zdPlot <- function(zeta1, delta1, parmRange, rescaleParms = c(1, 1), targetVals = c(0), targetPch = c(0), taus = NULL,
    offset = 5, cW = sqrt(5), cB = sqrt(5), n.pts = 9, N = 201, autoAdjZeta = F, zetaRange = NULL, deltaRange = NULL,
    cex = 1, zInflator = 1, debug = F, ...) {
    # the function that plots 'danger zones' use targetVals to show where 'upper bd' might be or where eta wd
    # be if we had unbiased Y eqn. recale sd parms (y.w,y.b,z.w,z.b) defaults to no rescale.  O/w feed s.d.s
    # based on '0' models for z,y(ucm) NOTE: since between ests aren't used in the plots, the gpSize is never
    # used by this function and so is not a param.

    pfunct = function(x, y, z, cols0, target, targetPch, taus, cex, ...) {
        lattice::panel.levelplot(x, y, z, ...)
        lattice::panel.abline(h = 0, col.line = 1)
        lattice::panel.abline(v = 0, col.line = 1)
        lattice::panel.points(x = delta1, y = zeta1, pch = 16, cex = 0.75 * cex, col.symbol = cols0)
        lattice::panel.points(x = target[, 1], y = target[, 2], pch = targetPch, lwd = 3, cex = 1.4 * cex,
            col = 1)
        if (!is.null(taus)) {
            len <- length(taus$taus[, 1])
            for (i in 1:len) {
                if (!is.na(taus$switch[i])) {
                  if (taus$switch[i]) {
                    lattice::panel.text(x = target[i, 1] + taus$offset, y = target[i, 2], label = bquote(tilde(tau)[o] ==
                      .(round(taus$taus[, 1], 2)[i])), adj = 0, cex = 0.75 * cex)
                  } else lattice::panel.text(x = target[i, 1] - taus$offset, y = target[i, 2], label = bquote(tilde(tau)[w] ==
                    .(round(taus$taus[, 2], 2)[i])), adj = 1, cex = 0.75 * cex)
                }
            }
        }
    }

    zeta <- seq(min(zeta1), max(zeta1), length = N)  #default
    if (!is.null(zetaRange)) {
        # override everything
        zeta <- seq(zetaRange[1], zetaRange[2], length = N)
    }

    delta <- seq(min(delta1), max(delta1), length = N)  #default
    if (!is.null(deltaRange)) {
        # override
        delta <- seq(deltaRange[1], deltaRange[2], length = N)
    }
    bdiff <- outer(zeta, delta, FUN = "olsMwin", cW = cW, cB = cB) * zInflator

    rng.bd <- max(abs(min(bdiff)), abs(max(bdiff)))
    at.pts <- round(seq(-rng.bd, +rng.bd, length = n.pts), 2)
    at.pts2 <- seq(-rng.bd, +rng.bd, length = 1 + 2 * (n.pts - 1))

    spec.lbl <- as.character(at.pts)
    spec.lbl[1] <- "OLS\nbetter"
    spec.lbl[length(spec.lbl)] <- "Within\nbetter"

    # use sigma for y,z, b & w, to scale zeta & delta product params
    zeta <- zeta/prod(rescaleParms[1])
    delta <- delta/prod(rescaleParms[2])
    # need to rescale zeta1 & delta1 too....
    zeta1 <- zeta1/prod(rescaleParms[1])
    delta1 <- delta1/prod(rescaleParms[2])
    # locate target point (for beta or gamma):
    nLocs <- length(targetVals)
    locs <- rep(NA, nLocs)  #initialize
    gap <- parmRange[2] - parmRange[1]
    for (i in 1:nLocs) {
        locs[i] <- which.min(abs(parmRange - targetVals[i]))
        if (locs[i] == 1 && targetVals[i] < parmRange[1] - gap) {
            locs[i] <- NA  #boundary cond - drop
            cat("Below minimum range, eta=", targetVals[i], "point beyond", x[1], y[1])
        }
        if (locs[i] == N && targetVals[i] > parmRange[1] + gap) {
            locs[i] <- NA  #boundary cond - drop
            cat("Above minimum range, eta=", targetVals[i], "point beyond", x[N], y[N])
        }
    }
    targetPts <- matrix(cbind(delta1[locs], zeta1[locs]), nLocs, 2, byrow = F)  #byrow=F for when you pass more than one point...
    # here is where we get vals for plausible taus for the plot - evaluate at targetVals these are ests. of tau
    # based on zeta,delta
    if (!is.null(taus)) {
        corr.tau.o <- taus$ols - olsBias(zeta[locs], delta[locs], cW = cW, cB = cB)
        corr.tau.w <- taus$win - winBias(zeta[locs], cW = cW)
        tau.switch <- sign(olsMwin(zeta[locs], delta[locs], cW = cW, cB = cB)) <= 0
        if (debug) {
            cat("Bias-corrected (model-based) tau evaluated at specified points (least abs bias indicated on plot):\n")
            dstr <- cbind(zeta[locs], delta[locs], corr.tau.o, corr.tau.w, tau.switch)
            dimnames(dstr) <- list(NULL, c("delta", "zeta", "OLS-based", "Within-based", "OLS-better?"))
            print(round(dstr, digits = 2))
        }
        tau.to.plot <- list(taus = cbind(corr.tau.o, corr.tau.w), offset = offset * (delta[2] - delta[1]),
            switch = tau.switch)
    } else {
        tau.to.plot <- NULL  #for passing to levelplot pfunct
    }
    df <- data.frame(x = rep(delta, N), y = rep(zeta, each = N), z = c(t(bdiff)))

    # alt version of cols0=grey(1-abs(parmRange)/max(abs(parmRange)))
    # add argument to change colors:
    # col.regions = grDevices::colorRampPalette(c("dark orange", "white", "dark blue"), space = "Lab")
    # add grDevices to DESCRIPTION to use
    lattice::levelplot(z ~ x * y, data = df, at = at.pts2, colorkey = list(at = at.pts, labels = list(at = at.pts,
        label = spec.lbl, cex = 0.75 * cex)), panel = pfunct, row.values = delta, column.values = zeta, xlab = list(label = expression(delta^{
        yz
    }), cex = 0.85 * cex), ylab = list(label = expression(zeta^{
        yz
    }), cex = 0.85 * cex), zeta1 = zeta1, delta1 = delta1, cols0 = 8, target = targetPts, targetPch = targetPch,
        taus = tau.to.plot, cex = cex, scales = list(x = list(cex = 0.85 * cex), y = list(cex = 0.85 * cex),
            cex = cex), ...)
}









