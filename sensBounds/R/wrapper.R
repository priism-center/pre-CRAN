#' Function to get model fits and plotting parameters
#' @export

#country$pop[country$pop<=0] <- NA

# Rename to something different than package name?
sensBounds <- function(
  formula = score ~ num_books + factor(sex) + word_knowl + homework + factor(type) + log(pop),
  grouplevel = ~school,
  data = italy,
  plot = TRUE) {
  # Add center and standardize options? These are set to TRUE for makeCdat

  # Prep data and formulas
  dd <- makeCdat(fmla = formula,
                 data = data,
                 groupFmla = grouplevel,
                 cntrTreat = TRUE,
                 stdz = TRUE,
                 tol = 0.01)

  # Core call to run the sensitivity analysis
  mdl.fit <- runModels(outcome = dd$fmlaY,
                       treatment = dd$fmlaZ,
                       level1.pred = dd$fmlaX,
                       level2.pred = dd$fmlaW,
                       group = dd$group,
                       data = dd$cdat)


  # Plot Results:
  #first gather pre-plot info:
  gpSize <- geoMeanGroupSize(dd$cdat$school)
  ppParm <- prePlotParams(mdl.fit, nGridPoints = 201, tau.max = 1, gpSize = gpSize)

  #report stability of calcs used in the determination of confounding line
  cat("condition number (for matrix used to identify line in confounding space):",
      round(ppParm$condNum, 2))
  # ADD WARNING HERE?




  return(list(model_fit = mdl.fit,
              plot_parameters = ppParm,
              model_name = deparse(substitute(data))))

}



# # Print Results:
# printResults(mdlFit = model_fit, mdlName = model_name, digits = 3, debug = FALSE)
# # From paper
# printResults(mdl.fit, cnames[ii], digits = 3)
#
#
#
#
#
# # Plot Results:
#
# #plot params set for .png inlcuded in LaTeX file; adjust as nec.
#   lcex <- 3
#   pcex <- 2
#   png(paste(data,"png", sep = ".", collapse = ""), width = pcex*480, height = pcex*480)
#   tpch <- c(0, 4, 1, 3)
#   pObj <- extractParams(mdl.fit) # to get taus from two model fits.
#   taus <- list(ols = pObj$tau.ols[2], win = pObj$tau.w[2])  #index 2 catches the CWC version of treatment Z.
#
# plot(zdPlot(ppParm$zetaDeltaMat[,"zeta"],
#                             ppParm$zetaDeltaMat[,"delta"],
#                             ppParm$parmRange,
#                             rescaleParms = c(1,1),
#                             targetVals = ppParm$bndVals,
#                             targetPch = tpch,
#                             taus = taus,
#                             cW = pObj$sigs[2,1],
#                             cB = pObj$sigs[2,2],
#                             cex = lcex))
#   dev.off()
#
#
#
