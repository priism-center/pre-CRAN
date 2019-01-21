#--------Read data and create new variables for testing-----------

classroom <- read.csv("classroom.csv")

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



#----------Test setup----------------

# Run case
test_case <- sensBounds(formula = test_formula, 
                        grouplevel = ~schoolid, 
                        data = classroom, 
                        condition_max = 1000)

# Checks
# Formulas (part of sensBounds object)

# ~'outcome var'
test_case[["model_fit"]][["outcome"]]

# ~'treatment var'
test_case[["model_fit"]][["treatment"]]

# ~'all group-level predictors (X.mn and W vars, with every level of factor vars except ref group)'
test_case[["model_fit"]][["level2.pred"]]


# Model outputs and ICCs from printResults
printResults(sbobject = test_case, digits = 3, debug = FALSE)
# check the model formulas correct?



# Plot prints properly
sensPlot(sbobject = test_case, cex = 2)



#---------Test Cases---------
# Outcome variable: mathgain
# Treatment: ses
# X = individual level variable (student); W = group level variable (school)


#----------Additive relationship, all continuous--------
# 1 continuous X, 0 W
test_formula <- mathgain ~ ses + mathkind

# 1 continuous X and 1 continous W
test_formula <- mathgain ~ ses + mathkind + housepov


### BROKEN ###
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
### WORKS ### factor treated as continuous and multiplied beforehand
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





