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
      country$pop[country$pop<=0] <- NA  #can't use a school with 0 popn
      dd <- makeCdat(score~num_books + factor(sex) + word_knowl + homework + factor(type) + log(pop), data=country, groupFmla=~school,stdz=T)
      #core call to run the sensitivity analysis. 
      mdl.fit<-runModels(outcome=dd$fmlaY, treatment=dd$fmlaZ, level1.pred = dd$fmlaX, level2.pred = dd$fmlaW, group = dd$group, data=dd$cdat)

      # Print Results: 
      printResults(mdl.fit,cnames[ii],digits=3)
      
      # Plot Results: 
      #first gather pre-plot info:
      gpSize <- geoMeanGroupSize(dd$cdat$school)
      ppParm <- prePlotParams(mdl.fit,nGridPoints=201,tau.max=1,gpSize=gpSize)
      
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
    ## score ~ `factor(sex)0` + word_knowl + homework + `factor(sex)0.mn` +  
    ##     word_knowl.mn + homework.mn + `factor(type)1` + `factor(type)0` +  
    ##     `factor(type)4` + `factor(type)3` + `log(pop)` + num_books +  
    ##     num_books.mn + (1 | school)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4402.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.6957 -0.5823  0.0549  0.6087  8.1771 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  school   (Intercept) 0.04188  0.2047  
    ##  Residual             0.54263  0.7366  
    ## Number of obs: 1913, groups:  school, 106
    ## 
    ## Fixed effects:
    ##                    Estimate Std. Error t value
    ## (Intercept)       -0.001821   0.027802  -0.066
    ## `factor(sex)0`     0.059338   0.017471   3.396
    ## word_knowl         0.624077   0.018964  32.909
    ## homework          -0.026672   0.017982  -1.483
    ## `factor(sex)0.mn`  0.045266   0.097523   0.464
    ## word_knowl.mn      0.132407   0.059534   2.224
    ## homework.mn        0.005156   0.076897   0.067
    ## `factor(type)1`   -0.003820   0.031389  -0.122
    ## `factor(type)0`   -0.057647   0.035073  -1.644
    ## `factor(type)4`    0.028296   0.026821   1.055
    ## `factor(type)3`   -0.028685   0.026877  -1.067
    ## `log(pop)`         0.024711   0.026698   0.926
    ## num_books          0.077116   0.019656   3.923
    ## num_books.mn       0.523222   0.065159   8.030

    ## 
    ## Correlation matrix not shown by default, as p = 14 > 12.
    ## Use print(summary(mdlFit$mlm1.y), correlation=TRUE)  or
    ##   vcov(summary(mdlFit$mlm1.y))    if you need it

    ## [1] "OLS regression model fit:"
    ## 
    ## Call:
    ## lm(formula = Yfmla1.lm, data = data)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -6.6861 -0.4713  0.0560  0.4759  6.1461 
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                  9.110e-17  1.764e-02   0.000  1.00000    
    ## `factor(sex)0`               5.958e-02  1.830e-02   3.256  0.00115 ** 
    ## word_knowl                   6.114e-01  1.980e-02  30.883  < 2e-16 ***
    ## homework                    -2.543e-02  1.883e-02  -1.351  0.17701    
    ## `factor(sex)0.mn`            1.498e-01  6.959e-02   2.152  0.03150 *  
    ## word_knowl.mn                4.315e-01  4.391e-02   9.827  < 2e-16 ***
    ## homework.mn                  8.384e-02  5.230e-02   1.603  0.10911    
    ## `factor(type)1`              5.121e-03  1.942e-02   0.264  0.79207    
    ## `factor(type)0`             -1.217e-02  2.122e-02  -0.573  0.56647    
    ## `factor(type)4`              1.369e-02  1.800e-02   0.761  0.44680    
    ## `factor(type)3`              2.823e-03  1.802e-02   0.157  0.87551    
    ## `log(pop)`                  -1.085e-02  1.924e-02  -0.564  0.57283    
    ## I(num_books + num_books.mn)  1.463e-01  1.870e-02   7.824 8.44e-15 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7715 on 1900 degrees of freedom
    ## Multiple R-squared:  0.4086, Adjusted R-squared:  0.4048 
    ## F-statistic: 109.4 on 12 and 1900 DF,  p-value: < 2.2e-16
    ## 
    ## [1] "Table 1 for:  scotland"
    ##      Within Between   OLS
    ## tau0  0.201   0.522 0.282
    ## tau1  0.077   0.523 0.146
    ## diff  0.124  -0.002 0.136
    ## [1] "ICCs for models:"
    ##       [,1]  [,2]  [,3]
    ## [1,] 0.806 0.189 0.190
    ## [2,] 0.780 0.136 0.149
    ## condition number (for matrix used to identify line in confounding space): 136.43

    ## [1] "Intermediary Model Fits for:  sweden"
    ## [1] "Multilevel model fit:"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## score ~ `factor(sex)0` + word_knowl + homework + `factor(sex)0.mn` +  
    ##     word_knowl.mn + homework.mn + `factor(type)1` + `factor(type)0` +  
    ##     `factor(type)3` + `log(pop)` + num_books + num_books.mn +  
    ##     (1 | school)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 4585.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.7315 -0.7036  0.1034  0.7359  3.5117 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  school   (Intercept) 0.02397  0.1548  
    ##  Residual             0.71972  0.8484  
    ## Number of obs: 1790, groups:  school, 96
    ## 
    ## Fixed effects:
    ##                     Estimate Std. Error t value
    ## (Intercept)       -0.0021206  0.0262874  -0.081
    ## `factor(sex)0`     0.0960053  0.0207225   4.633
    ## word_knowl         0.4513929  0.0211106  21.382
    ## homework          -0.0163107  0.0216661  -0.753
    ## `factor(sex)0.mn`  0.1024270  0.1065520   0.961
    ## word_knowl.mn      0.7255855  0.0993624   7.302
    ## homework.mn        0.0472699  0.0721234   0.655
    ## `factor(type)1`    0.0176449  0.0317568   0.556
    ## `factor(type)0`   -0.0002535  0.0344856  -0.007
    ## `factor(type)3`   -0.0018161  0.0296256  -0.061
    ## `log(pop)`         0.0086300  0.0276483   0.312
    ## num_books          0.0771670  0.0216164   3.570
    ## num_books.mn       0.2773272  0.0871697   3.181

    ## 
    ## Correlation matrix not shown by default, as p = 13 > 12.
    ## Use print(summary(mdlFit$mlm1.y), correlation=TRUE)  or
    ##   vcov(summary(mdlFit$mlm1.y))    if you need it

    ## [1] "OLS regression model fit:"
    ## 
    ## Call:
    ## lm(formula = Yfmla1.lm, data = data)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1872 -0.5927  0.0854  0.6379  2.9291 
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                 -8.073e-16  2.037e-02   0.000    1.000    
    ## `factor(sex)0`               9.613e-02  2.105e-02   4.566 5.31e-06 ***
    ## word_knowl                   4.494e-01  2.143e-02  20.971  < 2e-16 ***
    ## homework                    -1.647e-02  2.201e-02  -0.748    0.454    
    ## `factor(sex)0.mn`            1.167e-01  8.662e-02   1.348    0.178    
    ## word_knowl.mn                7.881e-01  7.962e-02   9.898  < 2e-16 ***
    ## homework.mn                  1.375e-02  5.673e-02   0.242    0.808    
    ## `factor(type)1`             -8.885e-04  2.338e-02  -0.038    0.970    
    ## `factor(type)0`             -3.532e-02  2.462e-02  -1.434    0.152    
    ## `factor(type)3`             -1.131e-02  2.211e-02  -0.512    0.609    
    ## `log(pop)`                   1.808e-02  2.168e-02   0.834    0.404    
    ## I(num_books + num_books.mn)  9.238e-02  2.096e-02   4.407 1.11e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8619 on 1778 degrees of freedom
    ## Multiple R-squared:  0.2617, Adjusted R-squared:  0.2571 
    ## F-statistic: 57.29 on 11 and 1778 DF,  p-value: < 2.2e-16
    ## 
    ## [1] "Table 1 for:  sweden"
    ##      Within Between   OLS
    ## tau0  0.137   0.399 0.167
    ## tau1  0.077   0.277 0.092
    ## diff  0.060   0.121 0.075
    ## [1] "ICCs for models:"
    ##       [,1]  [,2]  [,3]
    ## [1,] 0.927 0.076 0.076
    ## [2,] 0.914 0.040 0.042
    ## condition number (for matrix used to identify line in confounding space): 382.56

    ## [1] "Intermediary Model Fits for:  italy"
    ## [1] "Multilevel model fit:"
    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## score ~ `factor(sex)0` + word_knowl + homework + `factor(sex)0.mn` +  
    ##     word_knowl.mn + homework.mn + `factor(type)1` + `factor(type)0` +  
    ##     `factor(type)5` + `factor(type)4` + `factor(type)3` + `factor(type)6` +  
    ##     `log(pop)` + num_books + num_books.mn + (1 | school)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 7484.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.6299 -0.5635  0.0821  0.6071  3.7840 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  school   (Intercept) 0.2342   0.4839  
    ##  Residual             0.4801   0.6929  
    ## Number of obs: 3313, groups:  school, 256
    ## 
    ## Fixed effects:
    ##                    Estimate Std. Error t value
    ## (Intercept)        0.019669   0.037528   0.524
    ## `factor(sex)0`     0.019787   0.012842   1.541
    ## word_knowl         0.409508   0.015326  26.721
    ## homework           0.008980   0.013780   0.652
    ## `factor(sex)0.mn`  0.074419   0.079954   0.931
    ## word_knowl.mn      0.703641   0.054762  12.849
    ## homework.mn        0.015981   0.051363   0.311
    ## `factor(type)1`    0.030847   0.040062   0.770
    ## `factor(type)0`   -0.023150   0.042196  -0.549
    ## `factor(type)5`   -0.039849   0.038835  -1.026
    ## `factor(type)4`    0.005024   0.040604   0.124
    ## `factor(type)3`    0.023431   0.040381   0.580
    ## `factor(type)6`   -0.002742   0.039349  -0.070
    ## `log(pop)`        -0.028529   0.038952  -0.732
    ## num_books          0.059368   0.014642   4.055
    ## num_books.mn       0.108739   0.060869   1.786

    ## 
    ## Correlation matrix not shown by default, as p = 16 > 12.
    ## Use print(summary(mdlFit$mlm1.y), correlation=TRUE)  or
    ##   vcov(summary(mdlFit$mlm1.y))    if you need it

    ## [1] "OLS regression model fit:"
    ## 
    ## Call:
    ## lm(formula = Yfmla1.lm, data = data)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6548 -0.5078  0.0626  0.5564  3.7109 
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                 -7.575e-16  1.418e-02   0.000  1.00000    
    ## `factor(sex)0`               2.023e-02  1.513e-02   1.338  0.18109    
    ## word_knowl                   4.061e-01  1.800e-02  22.561  < 2e-16 ***
    ## homework                     8.388e-03  1.623e-02   0.517  0.60532    
    ## `factor(sex)0.mn`            6.487e-02  4.119e-02   1.575  0.11537    
    ## word_knowl.mn                7.196e-01  2.593e-02  27.752  < 2e-16 ***
    ## homework.mn                  1.341e-02  3.035e-02   0.442  0.65856    
    ## `factor(type)1`             -2.115e-02  1.641e-02  -1.289  0.19755    
    ## `factor(type)0`             -3.912e-02  1.963e-02  -1.993  0.04630 *  
    ## `factor(type)5`             -4.493e-02  1.507e-02  -2.981  0.00290 ** 
    ## `factor(type)4`             -4.335e-02  1.599e-02  -2.711  0.00674 ** 
    ## `factor(type)3`             -1.400e-02  1.750e-02  -0.800  0.42378    
    ## `factor(type)6`             -1.476e-02  1.502e-02  -0.983  0.32572    
    ## `log(pop)`                  -3.125e-02  1.866e-02  -1.675  0.09399 .  
    ## I(num_books + num_books.mn)  8.123e-02  1.483e-02   5.478 4.63e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8162 on 3298 degrees of freedom
    ## Multiple R-squared:  0.3366, Adjusted R-squared:  0.3338 
    ## F-statistic: 119.5 on 14 and 3298 DF,  p-value: < 2.2e-16
    ## 
    ## [1] "Table 1 for:  italy"
    ##      Within Between   OLS
    ## tau0  0.118   0.334 0.219
    ## tau1  0.059   0.109 0.081
    ## diff  0.059   0.225 0.138
    ## [1] "ICCs for models:"
    ##       [,1]  [,2]  [,3]
    ## [1,] 0.752 0.282 0.273
    ## [2,] 0.735 0.216 0.227
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
