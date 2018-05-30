run lmer and lm models
======================

This is a script that produces model results for the ObsStudies paper.
Plots follow model runs.

    #chose country names
    #each name corresponds to a dataframe object in ObsStudies_IEA.RData 

    cnames <-c("scotland","sweden","italy")

Data Setup first
----------------

Common formulas set up
----------------------

Then run models, all in a loop one for each country selected.
-------------------------------------------------------------

    set.seed(12341)
    for (ii in 1:length(cnames)) {
        
      country <- get(cnames[ii])

      #BETWEEN ##### could be a class to a 'helper' function called IEApreprocess that takes logs, removes NAs, etc.
      ################################
      country$pop <- log(country$pop)  #transform towards symmetry
      country$pop[is.infinite(country$pop)] <- NA
      
      country$type <- factor(country$type)
      
      #remove NAs:
      cdat <- model.frame(score~num_books + sex + word_knowl + homework + type + pop + school,data=country, na.action=na.omit) 

      #set up indicator vars for factors (preferred over +factor(type) as this would require more complicated re-use of data with commands like model.matrix:
      
      factForm <- buildFactorFormula(cdat,"type")
      cdat <- factForm$dsn
      
      nPerSchool <- prod(table(cdat$school))^(1/length(unique(cdat$school))) # use geo-mean (looking for a common denominator in the bias formulas with unbalanced design; geomean an approx. proxy)

      
      ## Prepare the formulas for the model fit calls.  
      #indiv-lwvel preds
      fmlaX <- ~sex + word_knowl + homework
      #group-level preds (pop + type, but type is "separate indicators")
      fmlaW <- as.formula(paste("~pop",factForm$fmlaString,sep="",collapse=""))
      #set treatment & outcome
      fmlaZ <- ~num_books
      fmlaY <- ~score
      ####################################
      
      #fit models
      mdl.fit <- runModels(outcome=fmlaY, treatment=fmlaZ, level1.pred = fmlaX, level2.pred = fmlaW, group = ~school, data=cdat)
      printResults(mdl.fit,cnames[ii],digits=2)
      


      # Bounds calcs -

      ppParm <- prePlotParams(mdl.fit,nGridPoints=201,tau.max=1,gpSize=nPerSchool)
      
      #report stability of calcs used in the determination of confounding line
      cat("condition number (for matrix used to identify line in confounding space):",round(ppParm$condNum,2))

      #plot params set for .png inlcuded in LaTeX file; adjust as nec.
      lcex <- 3
      pcex <- 2
      png(paste(cnames[ii],"png",sep=".",collapse=""),width=pcex*480,height=pcex*480)
      tpch <- c(0,4,1,3)
      pObj <- extractParams(mdl.fit) # to get taus from two model fits.
      taus <- list(ols=pObj$tau.ols[2],win=pObj$tau.w[2])  #index 2 catches the CWC version of treatment Z.
      
      plot(zdPlot(ppParm$zetaDeltaMat[,"zeta"],ppParm$zetaDeltaMat[,"delta"],ppParm$parmRange,rescaleParms=c(1,1),targetVals=ppParm$bndVals,targetPch=tpch,taus=taus,cW=pObj$sigs[2,1],cB=pObj$sigs[2,2],cex=lcex))
      dev.off()
    }

    ## [1] "Intermediary Model Fits for:  scotland"
    ## [1] "Multilevel model fit:"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## score ~ num_books + num_books.mn + sex0 + word_knowl + homework +  
    ##     sex0.mn + word_knowl.mn + homework.mn + pop + type1 + type2 +  
    ##     type3 + type4 + (1 | idNew)
    ##    Data: 
    ## as.data.frame(cbind(cbind(YZdata[, c(Yname, IDname)], Znew, cbind(Xnew,  
    ##     W)), idNew))
    ## 
    ## REML criterion at convergence: 4402.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.6957 -0.5823  0.0549  0.6087  8.1771 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  idNew    (Intercept) 0.04188  0.2047  
    ##  Residual             0.54263  0.7366  
    ## Number of obs: 1913, groups:  idNew, 106
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error t value
    ## (Intercept)   -0.001821   0.027802   -0.07
    ## num_books      0.077116   0.019656    3.92
    ## num_books.mn   0.523222   0.065159    8.03
    ## sex0           0.059338   0.017471    3.40
    ## word_knowl     0.624077   0.018964   32.91
    ## homework      -0.026672   0.017982   -1.48
    ## sex0.mn        0.045266   0.097523    0.46
    ## word_knowl.mn  0.132407   0.059534    2.22
    ## homework.mn    0.005156   0.076896    0.07
    ## pop            0.024711   0.026698    0.93
    ## type1          0.047655   0.035129    1.36
    ## type2          0.061497   0.037415    1.64
    ## type3         -0.016690   0.026658   -0.63
    ## type4          0.048688   0.027656    1.76

    ## 
    ## Correlation matrix not shown by default, as p = 14 > 12.
    ## Use print(summary(mdl.fit$mlm1.y), correlation=TRUE)  or
    ##   vcov(summary(mdl.fit$mlm1.y))   if you need it

    ## [1] "OLS regression model fit:"
    ## 
    ## Call:
    ## lm(formula = formula(terms(Yfmla, data = cbind(Z, Xnew, W), allowDotAsName = T)), 
    ##     data = as.data.frame(cbind(YZdata, Xnew, W)))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -6.6861 -0.4713  0.0560  0.4759  6.1461 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)    9.964e-17  1.764e-02   0.000  1.00000    
    ## num_books      1.463e-01  1.870e-02   7.824 8.44e-15 ***
    ## sex0           5.958e-02  1.830e-02   3.256  0.00115 ** 
    ## word_knowl     6.114e-01  1.980e-02  30.883  < 2e-16 ***
    ## homework      -2.543e-02  1.883e-02  -1.351  0.17701    
    ## sex0.mn        1.498e-01  6.959e-02   2.152  0.03150 *  
    ## word_knowl.mn  4.315e-01  4.391e-02   9.827  < 2e-16 ***
    ## homework.mn    8.384e-02  5.230e-02   1.603  0.10911    
    ## pop           -1.085e-02  1.924e-02  -0.564  0.57283    
    ## type1          1.598e-02  2.177e-02   0.734  0.46281    
    ## type2          1.298e-02  2.264e-02   0.573  0.56647    
    ## type3          5.355e-03  1.800e-02   0.298  0.76607    
    ## type4          1.800e-02  1.835e-02   0.981  0.32673    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7715 on 1900 degrees of freedom
    ## Multiple R-squared:  0.4086, Adjusted R-squared:  0.4048 
    ## F-statistic: 109.4 on 12 and 1900 DF,  p-value: < 2.2e-16
    ## 
    ## [1] "Table 1 for:  scotland"
    ##      Within Between  OLS
    ## tau0   0.20    0.52 0.28
    ## tau1   0.08    0.52 0.15
    ## diff   0.12    0.00 0.14
    ## [1] "ICCs for models:"
    ##      [,1] [,2] [,3]
    ## [1,] 0.81 0.19 0.19
    ## [2,] 0.78 0.14 0.15
    ## condition number (for matrix used to identify line in confounding space): 136.43

    ## [1] "Intermediary Model Fits for:  sweden"
    ## [1] "Multilevel model fit:"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## score ~ num_books + num_books.mn + sex0 + word_knowl + homework +  
    ##     sex0.mn + word_knowl.mn + homework.mn + pop + type1 + type2 +  
    ##     type3 + (1 | idNew)
    ##    Data: 
    ## as.data.frame(cbind(cbind(YZdata[, c(Yname, IDname)], Znew, cbind(Xnew,  
    ##     W)), idNew))
    ## 
    ## REML criterion at convergence: 4585
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.7315 -0.7036  0.1034  0.7359  3.5117 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  idNew    (Intercept) 0.02397  0.1548  
    ##  Residual             0.71972  0.8484  
    ## Number of obs: 1790, groups:  idNew, 96
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error t value
    ## (Intercept)   -0.0021206  0.0262874  -0.081
    ## num_books      0.0771670  0.0216164   3.570
    ## num_books.mn   0.2773272  0.0871697   3.181
    ## sex0           0.0960053  0.0207225   4.633
    ## word_knowl     0.4513929  0.0211106  21.382
    ## homework      -0.0163107  0.0216661  -0.753
    ## sex0.mn        0.1024270  0.1065520   0.961
    ## word_knowl.mn  0.7255855  0.0993624   7.302
    ## homework.mn    0.0472699  0.0721234   0.655
    ## pop            0.0086300  0.0276483   0.312
    ## type1          0.0179358  0.0366207   0.490
    ## type2          0.0003057  0.0415866   0.007
    ## type3         -0.0016143  0.0317164  -0.051

    ## 
    ## Correlation matrix not shown by default, as p = 13 > 12.
    ## Use print(summary(mdl.fit$mlm1.y), correlation=TRUE)  or
    ##   vcov(summary(mdl.fit$mlm1.y))   if you need it

    ## [1] "OLS regression model fit:"
    ## 
    ## Call:
    ## lm(formula = formula(terms(Yfmla, data = cbind(Z, Xnew, W), allowDotAsName = T)), 
    ##     data = as.data.frame(cbind(YZdata, Xnew, W)))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1872 -0.5927  0.0854  0.6379  2.9291 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -5.724e-16  2.037e-02   0.000    1.000    
    ## num_books      9.238e-02  2.096e-02   4.407 1.11e-05 ***
    ## sex0           9.613e-02  2.105e-02   4.566 5.31e-06 ***
    ## word_knowl     4.494e-01  2.143e-02  20.971  < 2e-16 ***
    ## homework      -1.647e-02  2.201e-02  -0.748    0.454    
    ## sex0.mn        1.167e-01  8.662e-02   1.348    0.178    
    ## word_knowl.mn  7.881e-01  7.962e-02   9.898  < 2e-16 ***
    ## homework.mn    1.375e-02  5.673e-02   0.242    0.808    
    ## pop            1.808e-02  2.168e-02   0.834    0.404    
    ## type1          3.964e-02  2.821e-02   1.405    0.160    
    ## type2          4.259e-02  2.969e-02   1.434    0.152    
    ## type3          1.681e-02  2.476e-02   0.679    0.497    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8619 on 1778 degrees of freedom
    ## Multiple R-squared:  0.2617, Adjusted R-squared:  0.2571 
    ## F-statistic: 57.29 on 11 and 1778 DF,  p-value: < 2.2e-16
    ## 
    ## [1] "Table 1 for:  sweden"
    ##      Within Between  OLS
    ## tau0   0.14    0.40 0.17
    ## tau1   0.08    0.28 0.09
    ## diff   0.06    0.12 0.07
    ## [1] "ICCs for models:"
    ##      [,1] [,2] [,3]
    ## [1,] 0.93 0.08 0.08
    ## [2,] 0.91 0.04 0.04
    ## condition number (for matrix used to identify line in confounding space): 382.56

    ## [1] "Intermediary Model Fits for:  italy"
    ## [1] "Multilevel model fit:"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## score ~ num_books + num_books.mn + sex0 + word_knowl + homework +  
    ##     sex0.mn + word_knowl.mn + homework.mn + pop + type1 + type2 +  
    ##     type3 + type4 + type5 + type6 + (1 | idNew)
    ##    Data: 
    ## as.data.frame(cbind(cbind(YZdata[, c(Yname, IDname)], Znew, cbind(Xnew,  
    ##     W)), idNew))
    ## 
    ## REML criterion at convergence: 7483.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.6299 -0.5635  0.0821  0.6071  3.7840 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  idNew    (Intercept) 0.2342   0.4839  
    ##  Residual             0.4801   0.6929  
    ## Number of obs: 3313, groups:  idNew, 256
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)    0.01967    0.03753   0.524
    ## num_books      0.05937    0.01464   4.055
    ## num_books.mn   0.10874    0.06087   1.786
    ## sex0           0.01979    0.01284   1.541
    ## word_knowl     0.40951    0.01533  26.721
    ## homework       0.00898    0.01378   0.652
    ## sex0.mn        0.07442    0.07995   0.931
    ## word_knowl.mn  0.70364    0.05476  12.849
    ## homework.mn    0.01598    0.05136   0.311
    ## pop           -0.02853    0.03895  -0.732
    ## type1          0.05213    0.04058   1.285
    ## type2          0.03424    0.06240   0.549
    ## type3          0.04487    0.03678   1.220
    ## type4          0.02678    0.04168   0.642
    ## type5         -0.02294    0.04090  -0.561
    ## type6          0.01339    0.04160   0.322

    ## 
    ## Correlation matrix not shown by default, as p = 16 > 12.
    ## Use print(summary(mdl.fit$mlm1.y), correlation=TRUE)  or
    ##   vcov(summary(mdl.fit$mlm1.y))   if you need it

    ## [1] "OLS regression model fit:"
    ## 
    ## Call:
    ## lm(formula = formula(terms(Yfmla, data = cbind(Z, Xnew, W), allowDotAsName = T)), 
    ##     data = as.data.frame(cbind(YZdata, Xnew, W)))
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6548 -0.5078  0.0626  0.5564  3.7109 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -3.930e-16  1.418e-02   0.000   1.0000    
    ## num_books      8.123e-02  1.483e-02   5.478 4.63e-08 ***
    ## sex0           2.023e-02  1.513e-02   1.338   0.1811    
    ## word_knowl     4.061e-01  1.800e-02  22.561  < 2e-16 ***
    ## homework       8.388e-03  1.623e-02   0.517   0.6053    
    ## sex0.mn        6.487e-02  4.119e-02   1.575   0.1154    
    ## word_knowl.mn  7.196e-01  2.593e-02  27.752  < 2e-16 ***
    ## homework.mn    1.341e-02  3.035e-02   0.442   0.6586    
    ## pop           -3.125e-02  1.866e-02  -1.675   0.0940 .  
    ## type1          1.481e-02  1.913e-02   0.774   0.4388    
    ## type2          5.785e-02  2.902e-02   1.993   0.0463 *  
    ## type3          2.223e-02  1.848e-02   1.203   0.2292    
    ## type4         -6.594e-03  1.953e-02  -0.338   0.7357    
    ## type5         -1.637e-02  1.781e-02  -0.919   0.3583    
    ## type6          1.250e-02  1.775e-02   0.704   0.4813    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8162 on 3298 degrees of freedom
    ## Multiple R-squared:  0.3366, Adjusted R-squared:  0.3338 
    ## F-statistic: 119.5 on 14 and 3298 DF,  p-value: < 2.2e-16
    ## 
    ## [1] "Table 1 for:  italy"
    ##      Within Between  OLS
    ## tau0   0.12    0.33 0.22
    ## tau1   0.06    0.11 0.08
    ## diff   0.06    0.23 0.14
    ## [1] "ICCs for models:"
    ##      [,1] [,2] [,3]
    ## [1,] 0.75 0.28 0.27
    ## [2,] 0.74 0.22 0.23
    ## condition number (for matrix used to identify line in confounding space): 80.14

Replicate Danger Zones Plots from Obs. Studies (2018) paper
-----------------------------------------------------------

    require(lattice)
    #used to generate figure 1 in the Observational Studies paper

    #set up a grid
    N <- 201
    lbd <- -4; ubd <- 4
    #set grid values
    zeta <- seq(lbd,ubd,length=N)
    delta <- seq(lbd,ubd,length=N)

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

    ## quartz_off_screen 
    ##                 2

    png("dz_1_3.png",width=480*2,height=480*2)
    plot(levelplot(t(bd.1.3), at=seq(-4,4,length=17),colorkey=list(at=at.pts,labels=list(labels=spec.lbl,cex=lcex)),scales=list(x=list(cex=lcex),y=list(cex=lcex)),panel=pfunct,row.values=zeta,column.values=delta,xlab=list(expression(delta^{yz}),cex=lcex),ylab=list(expression(zeta^{yz}),cex=lcex)))
    dev.off()

    ## quartz_off_screen 
    ##                 2

    png("dz_3_1.png",width=480*2,height=480*2)
    plot(levelplot(t(bd.3.1), at=seq(-4,4,length=17),colorkey=list(at=at.pts,labels=list(labels=spec.lbl,cex=lcex)),scales=list(x=list(cex=lcex),y=list(cex=lcex)),panel=pfunct,row.values=zeta,column.values=delta,xlab=list(expression(delta^{yz}),cex=lcex),ylab=list(expression(zeta^{yz}),cex=lcex)))
    dev.off()

    ## quartz_off_screen 
    ##                 2

Done.
