

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
    cat("Condition number (for matrix used to identify line in confounding space):", round(ppParm$condNum,
        2), "\n")
    # Warning about high condition number
    if (ppParm$condNum > condition_max) {
      warning("Condition number exceeds maximum limit. Results may not be stable.")
    }




    return(list(model_fit = mdl.fit, plot_parameters = ppParm, model_name = deparse(substitute(data))))

}



#---------PRINT RESULTS---------

#' Print results
#' @export

printResults <- function(sbobject, digits = 3, debug = FALSE) {
    # extract variance comps & some bias diffs
    pObj <- extractParams(sbobject$model_fit)

    # MLM Fit
    cat("Intermediary Model Fits for", sbobject$model_name, "Data \n \n")
    cat("Multilevel Model Fit \n")
    print(summary(sbobject$model_fit$mlm1.y))

    # OLS Fit
    cat("\n")
    cat("OLS Regression Model Fit")
    print(summary(sbobject$model_fit$ols1.y))

    # Taus Table
    cat("Table 1: Multilevel model-based estimates of treatment effect \n \n")
    rsltTab1 <- rbind(cbind(pObj$tau.w, pObj$tau.b, pObj$tau.ols), pObj$bias.diffs)
    dimnames(rsltTab1) <- list(c("Tau (without predictors)",
                                 "Tau (with predictors)",
                                 "Difference (change in bias)"),
                               c("Within", "Between", "OLS"))
    print(round(rsltTab1, digits = digits))

    # ICC Table
    cat("\n")
    cat("Table 2: Estimates of multilevel model-based within and between group-level variance components \n and intra-class correlation (ICC) of treatment \n \n")
    rsltTab2 <- cbind(pObj$sigs, pObj$sigs[, 2]/apply(pObj$sigs, 1, sum))
    dimnames(rsltTab2) <- list(c("Model without predictors",
                                 "Model with predictors"),
                               c("Within", "Between", "ICC"))
    print(round(rsltTab2, digits = digits))
    if (debug)
        print(paste("gy,vy,t.gz ", paste(round(sqrt(c(pObj$bndProdList$gy.vw.gy, pObj$sds.y.ucm$sd.alpha.y.ucm^2,
            2^2 * pObj$bndProdList$gz.vw.gz)), digits = digits), collapse = ", "), sep = " "))
}





#---------PLOT RESULTS---------

#' Plots the boundaries
#' @export

sensPlot <- function(sbobject, cex = 2) {
  #lattice::lattice.options(default.theme = lattice::standard.theme(color = FALSE))
    # plot params set for .png inlcuded in LaTeX file; adjust as nec.
    lcex <- cex
    tpch <- c(0, 4, 1, 3)
    pObj <- extractParams(sbobject$model_fit)  # to get taus from two model fits.
    taus <- list(ols = pObj$tau.ols[2], win = pObj$tau.w[2])  #index 2 catches the CWC version of treatment Z.

    plot(zdPlot(sbobject$plot_parameters$zetaDeltaMat[, "zeta"], sbobject$plot_parameters$zetaDeltaMat[, "delta"],
        sbobject$plot_parameters$parmRange, rescaleParms = c(1, 1), targetVals = sbobject$plot_parameters$bndVals,
        targetPch = tpch, taus = taus, cW = pObj$sigs[2, 1], cB = pObj$sigs[2, 2], cex = lcex))
}





