#--------Read data and create new variables for testing-----------

classroom <- read.csv("/Users/andrea/Documents/summer/Classroom_data_tests/classroom.csv")

### Creation of different types of variables ###
#check that every school has a sample
range(unique(classroom$schoolid))
unique(classroom$schoolid)==c(1:107)


#create extra school-level variables by averaging classroom vars
classroom$avgyears <- ave(classroom$yearstea, classroom$schoolid)
classroom$avgknow <- ave(classroom$mathknow, classroom$schoolid)
classroom$avgprep <- ave(classroom$mathprep, classroom$schoolid)

#create factor school-level vars
classroom$yearsfact <- 0
classroom$yearsfact[classroom$avgyears <= 10] <- 1
classroom$yearsfact[classroom$avgyears > 10 & classroom$avgyears <= 20] <- 2
classroom$yearsfact[classroom$avgyears > 20] <- 3

classroom$knowfact <- 0
classroom$knowfact[classroom$avgknow <= -1] <- 1
classroom$knowfact[classroom$avgknow > -1 & classroom$avgknow <= 1] <- 2
classroom$knowfact[classroom$avgknow > 1] <- 3
classroom$knowfact[is.na(classroom$avgknow)] <- NA

classroom$prepfact <- 0
classroom$prepfact[classroom$avgprep > 2.75] <- 1

#create declared factor var X
classroom$sex.fact <- factor(classroom$sex)

classroom$mathkind.fact <- 0
classroom$mathkind.fact[classroom$mathkind <= 400] <- 1
classroom$mathkind.fact[classroom$mathkind > 400 & classroom$mathkind <= 500] <- 2
classroom$mathkind.fact[classroom$mathkind > 500] <- 3
classroom$mathkind.fact <- factor(classroom$mathkind.fact)



#create factor var X with characters (non-numeric)
classroom$sex.non <- 0
classroom$sex.non[classroom$sex==0] <- "M"
classroom$sex.non[classroom$sex==1] <- "F"

#create factor var W with characters (non-numeric)
classroom$yearsfact.non <- 0
classroom$yearsfact.non[classroom$yearsfact==1] <- "A"
classroom$yearsfact.non[classroom$yearsfact==2] <- "B"
classroom$yearsfact.non[classroom$yearsfact==3] <- "C"



#create X vars that are already GMC
classroom$mathkind.gmc <- 0
mathkind.gmc <- tapply(classroom$mathkind, INDEX = classroom$schoolid, FUN = function(x) scale(x, scale=F))

for (i in 1:length(mathkind.gmc)) {
  classroom$mathkind.gmc[classroom$schoolid == i] <- as.numeric(mathkind.gmc[[i]])
}

classroom$minority.gmc <- 0
minority.gmc <- tapply(classroom$minority, INDEX = classroom$schoolid, FUN = function(x) scale(x, scale=F))

for (i in 1:length(minority.gmc)) {
  classroom$minority.gmc[classroom$schoolid == i] <- as.numeric(minority.gmc[[i]])
}


#check which vars have NAs
apply(classroom, 2, function(x) sum(is.na(x)))
#only NAs in mathknow, avgknow, and knowfact



#---------Former setup--------



set.seed(12341)
#clean with makeCdat
dd <- makeCdat(mathgain~ses + , data=classroom, groupFmla=~schoolid,stdz=T)

#Stop to check formulas
View(dd)



#run models
mdl.fit<-runModels(outcome=dd$fmlaY, treatment=dd$fmlaZ, level1.pred = dd$fmlaX, level2.pred = dd$fmlaW, group = dd$group, data=dd$cdat)

# Print Results: 
printResults(mdl.fit,"classroom data",digits=3)

# Plot Results: 
#first gather pre-plot info:
gpSize <- geoMeanGroupSize(dd$cdat$schoolid)
ppParm <- prePlotParams(mdl.fit,nGridPoints=201,tau.max=1,gpSize=gpSize)

#report stability of calcs used in the determination of confounding line
cat("condition number (for matrix used to identify line in confounding space):",round(ppParm$condNum,2))

#plot params set for .png inlcuded in LaTeX file; adjust as nec.
lcex <- 3
pcex <- 2
png(paste("classroom_data","png",sep=".",collapse=""),width=pcex*480,height=pcex*480)
tpch <- c(0,4,1,3)
pObj <- extractParams(mdl.fit) # to get taus from two model fits.
taus <- list(ols=pObj$tau.ols[2],win=pObj$tau.w[2])  #index 2 catches the CWC version of treatment Z.

plot(zdPlot(ppParm$zetaDeltaMat[,"zeta"],ppParm$zetaDeltaMat[,"delta"],ppParm$parmRange,rescaleParms=c(1,1),targetVals=ppParm$bndVals,targetPch=tpch,taus=taus,cW=pObj$sigs[2,1],cB=pObj$sigs[2,2],cex=lcex))
dev.off()


#----------New setup----------------

# Run case
test_case <- sensBounds(formula = test_formula, grouplevel = ~schoolid, data = classroom)
printResults(sbobject = test_case, digits = 3, debug = FALSE)
sensPlot(sbobject = test_case, cex = 2)

# Checks


#---------Test Cases---------
# Outcome variable: mathgain
# Treatment: ses
# X = individual level variable (student); W = group level variable (school)


#----------Additive relationship, all continuous--------
# 1 continuous X, 0 W
test_formula <- mathgain ~ ses + mathkind

# 1 continuous X and 1 continous W
test_formula <- mathgain ~ ses + mathkind + housepov

# 0 X, 1 continuous W 
test_formula <- mathgain ~ ses + housepov


# Error in parse(text = x, keep.source = FALSE) : 
#   <text>:2:0: unexpected end of input
# 1: ~
#   ^ 

# Warning message to add to makeCdat:
# if (b.varWithinGroup = NULL) stop("Warning: X NULL")
# Or change makeCdat so that on condition of b.varWithinGroup = NULL, formula is ~X.NULL and warning either given here or in biasAmp?


#----------Additive relationship, with factors--------

# 1 factor X, 0 W
test_formula <- mathgain ~ ses + factor(sex)

# 1 continuous X and 1 factor W
test_formula <- mathgain ~ ses + mathkind + factor(yearsfact)

# 1 factor X (classed, but undeclared in formula) and 0 W
test_formula <- mathgain ~ ses + sex.fact
test_formula <- mathgain ~ ses + mathkind.fact


# 1 factor X (character) and 0 W 
test_formula <- mathgain ~ ses + factor(sex.non)


# 1 factor X (character, undeclared in formula) and 0 W 
test_formula <- mathgain ~ ses + sex.non

# 1 cont. X and 1 factor W (character, undeclared in formula) 
test_formula <- mathgain ~ ses + mathkind + yearsfact.non




#----------Interaction terms-------
### BROKEN ###
# Error or Y var gets added to end
# None work unless otherwise noted

# Continuous*continuous
test_formula <- mathgain ~ ses + mathkind + mathkind:sex
### WORKS ###
test_formula <- mathgain ~ ses + mathkind + I(mathkind*sex)

# Continuous*factor
test_formula <- mathgain ~ ses + mathkind + mathkind:sex.fact
test_formula <- mathgain ~ ses + mathkind + mathkind*sex.fact
test_formula <- mathgain ~ ses + mathkind + mathkind*factor(sex)
test_formula <- mathgain ~ ses + mathkind + I(mathkind*sex.fact)

# Factor*factor
test_formula <- mathgain ~ ses + mathkind + factor(minority):factor(sex)
test_formula <- mathgain ~ ses + mathkind + factor(minority)*factor(sex)






#----------Polynomial terms-------
test_formula <- mathgain ~ ses + mathkind + I(mathkind)^2

test_formula <- mathgain ~ ses + mathkind + I(mathkind^3)

test_formula <- mathgain ~ ses + mathkind + I(mathkind^2) + I(mathkind^3)


#----------Mix-------



#----------Other cases------------

# 3 X, no within group variation (same as only W variables)
# Same error as above in additive relationships
test_formula <- mathgain ~ ses + avgprep + avgyears + housepov

# Error in parse(text = x, keep.source = FALSE) : 
#   <text>:2:0: unexpected end of input
# 1: ~
#   ^ 


# X already GMC 
test_formula <- mathgain ~ ses + mathkind.gmc
test_formula <- mathgain ~ ses + minority.gmc

### After running model ###
# Warning messages:
#   1: Some predictor variables are on very different scales: consider rescaling 
# 2: Some predictor variables are on very different scales: consider rescaling 

# Coefficient on the mean W (sex.gmc.mn or minority.gmc.mn) is very large (e+14)





