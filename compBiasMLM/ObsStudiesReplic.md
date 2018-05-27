run lmer and lm models
======================

Data Setup first
----------------

Then run models, all in a loop one for each country selected.
-------------------------------------------------------------

    #choice of country names
    cnames <-c("scotland","sweden","italy")

    set.seed(12341)
    for (ii in 1:length(cnames)) {
        
      country <- get(cnames[ii])

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

      nPerSchool <- mean(table(cdat$school))
      nPerSchool <- prod(table(cdat$school))^(1/length(unique(cdat$school))) # geo-mean

      #model spec. 
      
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

      #set treatment
      fmlaZ <- ~num_books
      fmlaY <- ~score

      #fit models
      mdl.fit <- runModels(outcome=fmlaY, treatment=fmlaZ, level1.pred = fmlaX, level2.pred = fmlaW, group = ~school, data=cdat)
      #extract variance comps & some bias diffs
      pObj <- extractParams(mdl.fit)
      #for ObsStudies paper:
      print(cnames[ii])
      print(summary(mdl.fit$mlm1.y))
      print(summary(mdl.fit$ols1.y))
      print(paste("Table 1 for: ",cnames[ii]))
      rsltTab1 <- rbind(cbind(pObj$tau.w,pObj$tau.b,pObj$tau.ols),pObj$bias.diffs)
      dimnames(rsltTab1) <- list(c("tau0","tau1","diff"),c("Within","Between","OLS"))
      print(round(rsltTab1,2))
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
      
      #reduce range if info from Beta warrants it
      nu.g.min <- min(recovParms3$zetaDeltaMat[,4])
      nu.g.max <- max(recovParms3$zetaDeltaMat[,4])
      #fix range here
      b.in.range <- recovParms4$zetaDeltaMat[,4] >= nu.g.min & recovParms4$zetaDeltaMat[,4] <= nu.g.max
      recovParms$zetaDeltaMat <- recovParms$zetaDeltaMat[b.in.range,]
      recovParms$parmRange <- recovParms$parmRange[b.in.range]
      recovParms$plausible.taus <- recovParms$plausible.taus[b.in.range]


      #plot params set for .png inlcuded in LaTeX file; adjust as nec.
      lcex <- 3
      pcex <- 2
      png(paste(cnames[ii],"png",sep=".",collapse=""),width=pcex*480,height=pcex*480)
      tvals <- c(gamma.at.beta,pObj$gy.vw.gz,nu.at.zeta.0,nu.at.delta.0)
      tpch <- c(0,4,1,3)
      taus <- list(ols=pObj$tau.ols[2],win=pObj$tau.w[2])
      confPts <- NULL #default
      if (!is.null(mdl.fit$bsSamps)) {
          if (dim(mdl.fit$bsSamps)[1]!=0) {
             rslt0 <- bsConfParms(mdl.fit$bsSamps,201,bnd.f,nPerSchool,tau.max=1,trim=0.025)$parms
             confPts <- cbind(as.vector(rslt0[,,2]),as.vector(rslt0[,,1]))
          }
      }
      
      rslt0 <- bsConfParmsSim(mdl.fit,nSims=250,nTarg=201,bnd.f=bnd.f,gpSize=nPerSchool,tau.max=1,debug=F,allowVarBool=c(T,T,T,rep(T,4)),line=F)
      confPts <- cbind(as.vector(rslt0$lines[,,2]),as.vector(rslt0$lines[,,1]))
      
      gamma.at.beta025 <- recovParms3$zetaDeltaMat[which.min(abs(recovParms3$zetaDeltaMat[,3]-quantile(rslt0$by.vx.bz,probs=c(0.025)))),4]
      gamma.at.beta975 <- recovParms3$zetaDeltaMat[which.min(abs(recovParms3$zetaDeltaMat[,3]-quantile(rslt0$by.vx.bz,probs=c(0.975)))),4]

      print("target for nu(by.vx.bz): ")
      print(nu.at.zeta.0)
      print("range for nu(by.vx.bz): ")
      print(c(gamma.at.beta025, gamma.at.beta,gamma.at.beta975),digits=3)
      
      plot(zdPlot(recovParms$zetaDeltaMat[,1],recovParms$zetaDeltaMat[,2],recovParms$parmRange,rescaleParms=rescaleParms,confPts=NULL,confPtsCol="darkslategrey",targetVals=tvals,targetPch=tpch,taus=taus,cW=pObj$sigs[2,1],cB=pObj$sigs[2,2],cex=lcex))
      dev.off()
    }

    ## [1] "scotland"
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
    ## Correlation of Fixed Effects:
    ##             (Intr) nm_bks nm_bk. sex0   wrd_kn homwrk sx0.mn wrd_k. hmwrk.
    ## num_books    0.000                                                        
    ## num_boks.mn  0.018  0.000                                                 
    ## sex0         0.000  0.004  0.000                                          
    ## word_knowl   0.000 -0.190  0.000  0.047                                   
    ## homework     0.000  0.020  0.000 -0.036  0.043                            
    ## sex0.mn     -0.010  0.000 -0.019  0.000  0.000  0.000                     
    ## wrd_knwl.mn -0.002  0.000 -0.374  0.000  0.000  0.000 -0.022              
    ## homework.mn -0.006  0.000 -0.128  0.000  0.000  0.000 -0.069 -0.018       
    ## pop          0.156  0.000  0.114  0.000  0.000  0.000  0.062 -0.063 -0.048
    ## type1        0.007  0.000  0.133  0.000  0.000  0.000 -0.137 -0.143  0.032
    ## type2       -0.001  0.000  0.250  0.000  0.000  0.000 -0.131 -0.092  0.105
    ## type3       -0.011  0.000 -0.094  0.000  0.000  0.000  0.005  0.123 -0.047
    ## type4       -0.019  0.000  0.169  0.000  0.000  0.000 -0.006 -0.085  0.026
    ##             pop    type1  type2  type3 
    ## num_books                              
    ## num_boks.mn                            
    ## sex0                                   
    ## word_knowl                             
    ## homework                               
    ## sex0.mn                                
    ## wrd_knwl.mn                            
    ## homework.mn                            
    ## pop                                    
    ## type1       -0.371                     
    ## type2       -0.437  0.559              
    ## type3       -0.090  0.096  0.107       
    ## type4       -0.147  0.238  0.291  0.038
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
    ## [1] 0.81 0.19 0.19 0.78 0.14 0.15
    ##      [,1] [,2] [,3]
    ## [1,] 0.81 0.19 0.19
    ## [2,] 0.78 0.14 0.15
    ## [1] "gy,vy,t.gz  0.104, 0.32, 0.413"
    ## cond num= 136.4313 
    ## [1] "target for nu(by.vx.bz): "
    ##     gamma 
    ## 0.0534951 
    ## [1] "range for nu(by.vx.bz): "
    ##    gamma    gamma    gamma 
    ## -0.00751  0.04139  0.08703 
    ## 0.1397986 1.047607 -0.08404467 1.482858 
    ## 0.3040622 1.211871 0.0802189 1.647122

    ## [1] "sweden"
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
    ## Correlation of Fixed Effects:
    ##             (Intr) nm_bks nm_bk. sex0   wrd_kn homwrk sx0.mn wrd_k. hmwrk.
    ## num_books    0.000                                                        
    ## num_boks.mn  0.004  0.000                                                 
    ## sex0         0.000  0.008  0.000                                          
    ## word_knowl   0.000 -0.132  0.000  0.074                                   
    ## homework     0.000 -0.010  0.000 -0.039 -0.050                            
    ## sex0.mn      0.008  0.000 -0.166  0.000  0.000  0.000                     
    ## wrd_knwl.mn  0.007  0.000 -0.120  0.000  0.000  0.000  0.090              
    ## homework.mn  0.031  0.000 -0.008  0.000  0.000  0.000 -0.013 -0.176       
    ## pop          0.027  0.000 -0.061  0.000  0.000  0.000  0.039  0.169  0.054
    ## type1        0.048  0.000 -0.265  0.000  0.000  0.000  0.092 -0.147 -0.108
    ## type2        0.036  0.000 -0.469  0.000  0.000  0.000  0.050 -0.114 -0.106
    ## type3        0.024  0.000 -0.150  0.000  0.000  0.000  0.033 -0.054 -0.124
    ##             pop    type1  type2 
    ## num_books                       
    ## num_boks.mn                     
    ## sex0                            
    ## word_knowl                      
    ## homework                        
    ## sex0.mn                         
    ## wrd_knwl.mn                     
    ## homework.mn                     
    ## pop                             
    ## type1       -0.192              
    ## type2       -0.272  0.655       
    ## type3       -0.219  0.486  0.506
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
    ## [1] 0.93 0.08 0.08 0.91 0.04 0.04
    ##      [,1] [,2] [,3]
    ## [1,] 0.93 0.08 0.08
    ## [2,] 0.91 0.04 0.04
    ## [1] "gy,vy,t.gz  0.203, 0.264, 0.389"
    ## cond num= 382.5636 
    ## [1] "target for nu(by.vx.bz): "
    ##      gamma 
    ## 0.02048167 
    ## [1] "range for nu(by.vx.bz): "
    ##   gamma   gamma   gamma 
    ## -0.0437  0.0115  0.0635 
    ## 0.3275606 0.2708298 0.07227225 0.5828489 
    ## 0.3335862 0.2768555 0.07829791 0.5888745

    ## [1] "italy"
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
    ## Correlation of Fixed Effects:
    ##             (Intr) nm_bks nm_bk. sex0   wrd_kn homwrk sx0.mn wrd_k. hmwrk.
    ## num_books    0.000                                                        
    ## num_boks.mn  0.050  0.000                                                 
    ## sex0         0.000  0.023  0.000                                          
    ## word_knowl   0.000 -0.148  0.000  0.014                                   
    ## homework     0.000 -0.029  0.000  0.009 -0.054                            
    ## sex0.mn     -0.032  0.000 -0.025  0.000  0.000  0.000                     
    ## wrd_knwl.mn  0.004  0.000 -0.261  0.000  0.000  0.000  0.009              
    ## homework.mn -0.010  0.000 -0.070  0.000  0.000  0.000  0.130 -0.065       
    ## pop          0.242  0.000 -0.036  0.000  0.000  0.000  0.016  0.117 -0.087
    ## type1        0.018  0.000 -0.158  0.000  0.000  0.000  0.031  0.060  0.076
    ## type2        0.073  0.000 -0.167  0.000  0.000  0.000  0.026 -0.197  0.069
    ## type3        0.034  0.000 -0.039  0.000  0.000  0.000  0.006 -0.073  0.126
    ## type4        0.028  0.000  0.056  0.000  0.000  0.000  0.039 -0.096  0.028
    ## type5        0.031  0.000 -0.111  0.000  0.000  0.000  0.000 -0.073  0.061
    ## type6        0.025  0.000  0.009  0.000  0.000  0.000  0.018 -0.061  0.040
    ##             pop    type1  type2  type3  type4  type5 
    ## num_books                                            
    ## num_boks.mn                                          
    ## sex0                                                 
    ## word_knowl                                           
    ## homework                                             
    ## sex0.mn                                              
    ## wrd_knwl.mn                                          
    ## homework.mn                                          
    ## pop                                                  
    ## type1       -0.291                                   
    ## type2       -0.611  0.491                            
    ## type3       -0.151  0.361  0.435                     
    ## type4       -0.346  0.369  0.502  0.360              
    ## type5       -0.290  0.321  0.442  0.300  0.316       
    ## type6       -0.340  0.305  0.428  0.282  0.324  0.263
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
    ## [1] 0.75 0.28 0.27 0.74 0.22 0.23
    ##      [,1] [,2] [,3]
    ## [1,] 0.75 0.28 0.27
    ## [2,] 0.74 0.22 0.23
    ## [1] "gy,vy,t.gz  0.423, 0.681, 0.535"
    ## cond num= 80.14383 
    ## [1] "target for nu(by.vx.bz): "
    ##    gamma 
    ## 0.111183 
    ## [1] "range for nu(by.vx.bz): "
    ##  gamma  gamma  gamma 
    ## 0.0581 0.0989 0.1356 
    ## 0.1207769 0.4215799 -0.07064312 0.6129999 
    ## 0.2532973 0.5541002 0.06187724 0.7455203
