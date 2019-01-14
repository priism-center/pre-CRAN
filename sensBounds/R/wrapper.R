

# country$pop[country$pop<=0] <- NA

#---------MAIN FUNCTION-------

#' Function to get model fits and plotting parameters
#' @export

sensBounds <- function(formula = score ~ num_books + factor(sex) + word_knowl + homework + factor(type) + log(pop),
    grouplevel = ~school, data = italy, condition_max = 1000) {
    # Add center and standardize options? These are set to TRUE for makeCdat

    # Prep data and formulas
    dd <- makeCdat(fmla = formula, data = data, groupFmla = grouplevel, cntrTreat = TRUE, stdz = TRUE, tol = 0.01)

    # Core call to run the sensitivity analysis
    mdl.fit <- runModels(outcome = dd$fmlaY, treatment = dd$fmlaZ, level1.pred = dd$fmlaX, level2.pred = dd$fmlaW,
        group = dd$group, data = dd$cdat)


    # Plot parameters: first gather pre-plot info:
    gpSize <- geoMeanGroupSize(dd$cdat$school)
    ppParm <- prePlotParams(mdl.fit, nGridPoints = 201, tau.max = 1, gpSize = gpSize)

    # report stability of calcs used in the determination of confounding line
    cat("condition number (for matrix used to identify line in confounding space):", round(ppParm$condNum,
        2))
    # ADD WARNING HERE?




    return(list(model_fit = mdl.fit, plot_parameters = ppParm, model_name = deparse(substitute(data))))

}



#---------PRINT RESULTS---------

#' Print results
#' @export

printResults <- function(sbobject, digits = 3, debug = FALSE) {
    # helper function to give the output from the models

    # extract variance comps & some bias diffs
    pObj <- extractParams(sbobject$model_fit)
    # used in ObsStudies paper:
    print(paste("Intermediary Model Fits for: ", sbobject$model_name))
    print("Multilevel model fit:")
    print(summary(sbobject$model_fit$mlm1.y))
    print("OLS regression model fit:")
    print(summary(sbobject$model_fit$ols1.y))
    print(paste("Table 1 for: ", sbobject$model_name))
    rsltTab1 <- rbind(cbind(pObj$tau.w, pObj$tau.b, pObj$tau.ols), pObj$bias.diffs)
    dimnames(rsltTab1) <- list(c("tau0", "tau1", "diff"), c("Within", "Between", "OLS"))
    print(round(rsltTab1, digits = digits))
    print("ICCs for models:")
    print(round(cbind(pObj$sigs, pObj$sigs[, 2]/apply(pObj$sigs, 1, sum)), digits = digits))
    if (debug)
        print(paste("gy,vy,t.gz ", paste(round(sqrt(c(pObj$bndProdList$gy.vw.gy, pObj$sds.y.ucm$sd.alpha.y.ucm^2,
            2^2 * pObj$bndProdList$gz.vw.gz)), digits = digits), collapse = ", "), sep = " "))
}





#---------PLOT RESULTS---------

#' Plots the boundaries
#' @export

sensPlot <- function(sbobject, cex = 2) {
    # plot params set for .png inlcuded in LaTeX file; adjust as nec.
    lcex <- cex
    tpch <- c(0, 4, 1, 3)
    pObj <- extractParams(sbobject$model_fit)  # to get taus from two model fits.
    taus <- list(ols = pObj$tau.ols[2], win = pObj$tau.w[2])  #index 2 catches the CWC version of treatment Z.

    plot(zdPlot(sbobject$plot_parameters$zetaDeltaMat[, "zeta"], sbobject$plot_parameters$zetaDeltaMat[, "delta"],
        sbobject$plot_parameters$parmRange, rescaleParms = c(1, 1), targetVals = sbobject$plot_parameters$bndVals,
        targetPch = tpch, taus = taus, cW = pObj$sigs[2, 1], cB = pObj$sigs[2, 2], cex = lcex))
}





