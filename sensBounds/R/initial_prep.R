# Contents of readDataFnsV4R0.R



# Add function documentation


#####
# The following line comes before the cleaning function. It gets rid of schools with population of zero. Should this happen before saving the datasets? Raw dataset source and cleaning code can be documented in data-raw/ and included in .Rbuildignore as ^data-raw$ so it isn't included in the built package.

#country$pop[country$pop<=0] <- NA  #can't use a school with 0 popn


#####
# implement rhs, lhs locally if library: formula.tools not avail.

# rhs <- function(x) { as.formula(paste0('~',strsplit(as.character(x),'~')[[1]][2])) }
# lhs <- function(x) { strsplit(as.character(x),'~')[[1]][1] } modified funct to be less IEA specific...

makeCdat <- function(fmla, data, groupFmla, cntrTreat = TRUE, stdz = TRUE, tol = 0.01) {

    # remove NAs & process the formula:
    groupName <- as.character(formula.tools::rhs(groupFmla))
    respName <- as.character(formula.tools::lhs(fmla))
    treatName <- attr(terms(fmla, response = F), "term.labels")[1]  # first predictor is treatment (standard approach)
    addResp <- as.formula(paste0("~.+", respName))
    fmlaInclResp <- as.formula(paste0(paste0("~", strsplit(as.character(fmla), "~")[[1]][2]), "+",
        respName))  #add response to end of fmla to extract from df
    addGroup <- as.formula(paste0("~.+", groupName))
    fmlaAll <- update(fmla, addGroup)
    cdat <- model.frame(fmlaAll, data = data, na.action = na.omit)
    ids <- model.frame(groupFmla, data = cdat)
    cdat[, groupName] <- NULL  # drop from data.frame now, so that IDs don't get weird
    cdat <- as.data.frame(model.matrix(fmlaInclResp, cdat)[, -1])  #forces all named factors to be indicators, going fwd; drops intercept, as usual

    # stdz first as scale is more obvious for eval'g tol (between var)
    if (stdz)
        cdat <- as.data.frame(sapply(cdat, scaleAll))  #still use scaleAll as there may be lurking factors (not a problem if they can be treated as numeric values - e.g., sex 0/1).
    meanByGroup <- aggregate(cdat, by = list(ids[, 1]), mean)  # for CWC
    varByGroup <- aggregate(cdat, by = list(ids[, 1]), var)
    b.varWithinGroup <- apply(varByGroup[, -1] > tol, 2, any, na.rm = T)  #these have at least some within variation.
    # might need to trap condition of no variation within group anywhere, e.g. at this point, the
    # treatment is first and the response is last in the df.  Exploit this structure.
    dimnames(meanByGroup)[[2]][1] <- groupName  # to faciliate merge.
    id.groupMeans <- meanByGroup[, 1]  #save these
    df.groupMeans <- meanByGroup[, -1][, b.varWithinGroup]  #only want those that vary (need to append .mn to name))  ## and need to drop id from this - it's first col.
    ctrdVarNames <- dimnames(df.groupMeans)[[2]]
    # at this point, last entries are outcome var, so prob. want to remove from fmla soon.
    groupVarNames <- paste0(ctrdVarNames, ".mn")
    df.groupMeans <- cbind(id.groupMeans, df.groupMeans)
    dimnames(df.groupMeans)[[2]] <- c(groupName, groupVarNames)
    df.centeredVars <- cbind(cdat[, b.varWithinGroup], ids)  #these will be centered
    df.X <- merge(df.centeredVars, df.groupMeans, by = groupName)
    df.Worig <- cdat[, !b.varWithinGroup, drop = F]  #original group level vars
    # hold out the response var (need to put back):
    respVar <- df.X[, respName]
    # hold out treatment var (in case need to put back):
    treatVar <- df.X[, treatName]
    # center Xs:
    df.X[, ctrdVarNames] <- df.X[, ctrdVarNames] - df.X[, groupVarNames]
    # put back response var - it shouldn't be de-meaned
    df.X[, respName] <- respVar
    if (!cntrTreat) {
        # put back mean if response variable, and flag set:
        df.X[, treatName] <- treatVar
    }
    datNew <- cbind(df.X, df.Worig)
    ## Prepare the formulas for the model fit calls. indiv-level preds need to use `` due to factor
    ## vars and logged vars...  should work, though clumsy
    XvarNames <- paste0(paste0("`", ctrdVarNames[-1]), "`")  #[-1] drops treatment var
    lenXvarNames <- length(XvarNames)
    fmlaX <- as.formula(paste0("~", paste(XvarNames[-lenXvarNames], collapse = "+")))  #[-len..] drops outcome var
    # group-level preds (pop + type, but type is 'separate indicators')
    lenGroupVarNames <- length(groupVarNames)
    WvarNames <- paste0(paste0("`", c(groupVarNames[-c(1, lenGroupVarNames)], dimnames(df.Worig)[[2]])),
        "`")  # -c(1,len) to drop outcome and treatment.
    fmlaW <- as.formula(paste0("~", paste(WvarNames, collapse = "+")))  #same drop of treatment (mean) var
    # set treatment & outcome
    fmlaZ <- as.formula(paste0("~", treatName))
    fmlaZ.mn <- as.formula(paste0("~", groupVarNames[1]))  # rather than reconstruct the .mn suffix.
    fmlaY <- as.formula(paste0("~", respName))

    ####################################
    return(list(cdat = datNew, fmlaY = fmlaY, fmlaZ = fmlaZ, fmlaZ.mn = fmlaZ.mn, fmlaX = fmlaX,
        fmlaW = fmlaW, group = groupFmla))
}

