#
dR2 <- function(NullModel, FullModel, B=1000, seed=NULL, weights=FALSE)
{
  # Given two nested models (NullModel << FullModel)
  # computes the increase in R^2 from the NullModel to the Full
  # and its p-value under a permutation procedure
  # weights: logical value indicating whether WLS shold be applied
  
  if(length(NullModel$residuals) != length(FullModel$residuals))
  {
    stop("Both models are not built on the same units. Check whether there are some missing values")
  }
  if(!is.null(seed))
  {
    set.seed(seed)
  }
  aux <- calcR2perm(NullModel, FullModel, B, weights)
  dR2.obs <- aux$dR2.obs
  dR2.perm <- aux$dR2.perm
  pvalue <- sum(dR2.perm > dR2.obs)/B
  out <- list("dR2"= dR2.obs, "pvalue"= pvalue)
  return(out)
}

dD <- function(NullModel, FullModel, B=1000, seed=NULL)
{
  # Given two nested models (NullModel << FullModel)
  # computes the increase in D from the NullModel to the Full
  # and its p-value under a permutation procedure
  
  if(length(NullModel$residuals) != length(NullModel$residuals))
  {
    stop("Both models are not built on the same units. Check whether there are some missing values")
  }
  if(!is.null(seed))
  {
    set.seed(seed)
  }
  # Comprobar que pasa con los NA
  aux <- calcDperm(NullModel, FullModel, B)
  dD.obs <- aux$dD.obs
  dD.perm <- aux$dD.perm
  pvalue <- sum(dD.perm > dD.obs)/B
  out <- list("dD"= dD.obs, "pvalue"= pvalue)
  return(out)
}

orderR2 <- function(data, yname, prsname = "PRS.")
{
  # Input: 
  #      data: data.frame with Trait and all PRS
  #      yname: character, name of the Trait variable
  #      prsname: character, common names for PRSs
  # Output:
  #      PRSs ordered by sum R^2
  
  pT <- dim(data)[2]
  indexy <- which(names(data)==yname)
  y <- data[, indexy]
  PRSindexes <- grep(prsname, colnames(data))
  nPRS <- length(PRSindexes)
  PRSnames <- names(data)[PRSindexes]
  indexrest <- (1:pT)[-c(indexy, PRSindexes)]
  xrest <- data[, indexrest]
  indexfactor <- which(sapply(xrest, is.factor))
  nF <- length(indexfactor)
  if(nF>2)
  {
    nF <- 2
    warning("Only interaction with respect to the two first factors will be considered")
  }
  nmodels <- switch(nF+1, 1, 2, 4)
  R <- matrix(NA, nrow=nPRS, ncol=nmodels,
              dimnames = list(paste0(prsname, 1:nPRS), paste0("Model", 1:nmodels)))
  
  for (i in 1:nPRS)
  {
    xi <- data.frame(y,  data[, PRSindexes[i]], xrest)
    xi <- na.omit(xi)
    names(xi) <- c(yname, PRSnames[i], names(xrest) )
    p <- dim(xi)[2]
    formula <- substitute(y ~ . , list(y=as.name(yname)))
    Mi <- lm(formula, data=xi)
    ximodel <- model.matrix(Mi)
    R[i, 1] <- summary(Mi)$r.squared
    if (nF==1) # There is an only factor
    {
      f <- indexfactor[1]
      geninteraction <- grep(names(xrest)[f], colnames(ximodel))
      ximodelextended <- cbind(ximodel, ximodel[, 2]*ximodel[, geninteraction])
      R[i, 2] <- R2NW(xi[, 1], ximodelextended)
    }else
    {
      if(nF==2)
      {
        f <- indexfactor[1:nF]
        geninteraction <- grep(names(xrest)[f[1]], colnames(ximodel))
        ximodelextended <- cbind(ximodel, ximodel[, 2]*ximodel[, geninteraction])
        last <- ximodelextended[, dim(ximodelextended)[2]]
        geninteraction <- grep(names(xrest)[f[2]], colnames(ximodel))
        R[i, 2] <- R2NW(xi[, 1], ximodelextended)
        ximodelextended <- cbind(ximodel, ximodel[, 2]*ximodel[, geninteraction])
        R[i, 3] <- R2NW(xi[, 1], ximodelextended)
        ximodelextended <- cbind(ximodelextended, last)
        R[i, 4] <- R2NW(xi[, 1], ximodelextended)
      }
    }
  }
  if(nF>=1)
  {
    total <- apply(R, 1, sum)
    R <- cbind(R, total)
    colnames(R)[nmodels+1] <- "Sum"
  }else
  {
    total <- R
  }
  
  o <- order(total, decreasing=TRUE)
  return(R[o,])
}

orderBin <- function(data, yname, prsname = "PRS.", statistic="D")
{
  # Input: 
  #      data: data.frame with binary Trait and all PRS
  #      yname: character, name of the Trait variable
  #      prsname: character, common names for PRSs
  #      statistic: "D" for Tjur's coefficient of dicrimination; "PseudoR2" for Nagelkerke's pseudo R^2
  # Output:
  #      PRSs ordered by sum of the selected statistic.
  
  
  fstatistic <- switch(statistic,
                       "D" = Dcalc,
                       "PseudoR2"= PseudoR2calc)
  pT <- dim(data)[2]
  indexy <- which(names(data)==yname)
  y <- data[, indexy]
  PRSindexes <- grep(prsname, colnames(data))
  nPRS <- length(PRSindexes)
  PRSnames <- names(data)[PRSindexes]
  indexrest <- (1:pT)[-c(indexy, PRSindexes)]
  xrest <- data[, indexrest]
  indexfactor <- which(sapply(xrest, is.factor))
  nF <- length(indexfactor)
  if(nF>2)
  {
    nF <- 2
    warning("Only interaction with respect to the two first factors will be considered")
  }
  nmodels <- switch(nF+1, 1, 2, 4)
  D <- matrix(NA, nrow=nPRS, ncol=nmodels,
              dimnames = list(paste0(prsname, 1:nPRS), paste0("Model", 1:nmodels)))
  
  for (i in 1:nPRS)
  {
    xi <- data.frame(y,  data[, PRSindexes[i]], xrest)
    xi <- na.omit(xi)
    yi <- xi[, 1]
    names(xi) <- c(yname, PRSnames[i], names(xrest) )
    p <- dim(xi)[2]
    formula <- substitute(y ~ . , list(y=as.name(yname)))
    Mi <- glm(formula, data=xi, family=binomial())
    ximodel <- model.matrix(Mi)
   # predMi <- Mi$fitted.values
    # aux <- tapply(predMi, y, mean)
    D[i, 1] <- fstatistic(yi, ximodel) #aux[2] - aux[1] # Coefficient of discrimination
    if (nF==1) # There is an only factor
    {
      f <- indexfactor[1]
      geninteraction <- grep(names(xrest)[f], colnames(ximodel))
      ximodelextended <- cbind(ximodel, ximodel[, 2]*ximodel[, geninteraction])
      D[i, 2] <- Dcalc(xi[, 1], ximodelextended)
    }else
    {
      if(nF==2)
      {
        f <- indexfactor[1:nF]
        geninteraction <- grep(names(xrest)[f[1]], colnames(ximodel))
        ximodelextended <- cbind(ximodel, ximodel[, 2]*ximodel[, geninteraction])
        last <- ximodelextended[, dim(ximodelextended)[2]]
        geninteraction <- grep(names(xrest)[f[2]], colnames(ximodel))
        D[i, 2] <- fstatistic(xi[, 1], ximodel) #Dcalc(xi[, 1], ximodelextended)
        ximodelextended <- cbind(ximodel, ximodel[, 2]*ximodel[, geninteraction])
        D[i, 3] <- fstatistic(xi[, 1], ximodelextended) #Dcalc(xi[, 1], ximodelextended)
        ximodelextended <- cbind(ximodelextended, last)
        D[i, 4] <-  fstatistic(xi[, 1], ximodelextended) #Dcalc(xi[, 1], ximodelextended)
      }
    }
  }
  if(nF>=1)
  {
    total <- apply(D, 1, sum)
    D <- cbind(D, total)
    colnames(D)[nmodels+1] <- "Sum"
  }else
  {
    total <- D
  }
  
  o <- order(total, decreasing=TRUE)
  return(D[o, ])
}


#----------------------------------------------------------------
# Auxiliar functions
R2W <- function(y, X, w){
  # X: design matrix
  # y: response variable
  # Output, coefdet: coefficient of determination
  aux <- stats::lm.wfit(X, y, w=w)
  SSreg <- sum((aux$fitted.values - mean(aux$fitted.values))^2)
  RSS <- sum(aux$residuals^2)
  coefdet <- SSreg/(SSreg+ RSS)
  return(coefdet)
}

R2NW <- function(y, X){
  # X: design matrix
  # y: response variable
  # Output, coefdet: coefficient of determination
  aux <- stats::lm.fit(X, y) 
  SSreg <- sum((aux$fitted.values - mean(aux$fitted.values))^2)
  RSS <- sum(aux$residuals^2)
  coefdet <- SSreg/(SSreg+ RSS)
  return(coefdet)
}

Dcalc <- function(y, X){
  # X: design matrix
  # y: response variable
  # Output, coefdet: discrimination coefficent
  M <- stats::glm.fit(X, y, family=binomial()) 
  pred <- M$fitted.values
  aux <- tapply(pred, y, mean)
  discrcoef <- aux[2] - aux[1]
  return(discrcoef)
}

PseudoR2calc <- function(y, X){
  # X: design matrix
  # y: response variable
  # Output, coefdet: Nagelkerke's pseudo R2
  M <- stats::glm.fit(X, y, family=binomial()) 
  n <- length(y)
  DevM <- M$deviance
  Dev0 <- M$null.deviance
  L1 <- exp(-DevM/2)
  L0 <- exp(-Dev0/2)
  coxsnell <- 1 - (L0/L1)^(2/n)
  pseudoR2 <- coxsnell/(1 - L0^(2/n))
  return(pseudoR2)
}

CoxSnellcalc <- function(y, X){
  # X: design matrix
  # y: response variable
  # Output, coefdet: Cox-Snell's pseudo R2
  M <- stats::glm.fit(X, y, family=binomial()) 
  n <- length(y)
  DevM <- M$deviance
  Dev0 <- M$null.deviance
  L1 <- exp(-DevM/2)
  L0 <- exp(-Dev0/2)
  pseudoR2 <- 1 - (L0/L1)^(2/n)
  return(pseudoR2)
}

Dmodel <- function(model){
  # model: logistic model
  # Output, coefdet: discrimination coefficent
  y <- model$y
  pred <- model$fitted.values
  aux <- tapply(pred, y, mean)
  discrcoef <- aux[2] - aux[1]
  return(discrcoef)
}


# Third step in Freedman & Lane
calcR2perm <- function(NullModel, FullModel, B, weights=weights)
{
  if (weights)
  {
    out <- calcR2permW(NullModel=NullModel, FullModel=FullModel, B=B)
  }else
  {
    out <- calcR2permNW(NullModel=NullModel, FullModel=FullModel, B=B)
  }
  return(out)
}

calcR2permNW <- function(NullModel, FullModel, B)
{
  Prednull <- NullModel$fitted.values
  Resnull <- NullModel$residuals
  n <- length(Prednull)
  Xbig <- model.matrix(FullModel)
  Xsmall <- model.matrix(NullModel)
  Y <- Prednull + Resnull
  
  dR2.obs <- R2NW(Y, Xbig) - R2NW(Y, Xsmall)
  dR2.perm <- rep(0, B)
  for (b in 1:B)
  {
    perm <- sample(n, replace=FALSE)
    Yperm <- Resnull[perm] + Prednull # Step 3 in Freedman & Lane
    dR2.perm[b] <- R2NW(Yperm, Xbig) - R2NW(Yperm, Xsmall)
  }
 
  out <- list(dR2.obs=dR2.obs, dR2.perm=dR2.perm)
  return(out)
}

calcR2permW <- function(NullModel, FullModel, B)
{
  # Weights for the Null and the Full Model are not the same.
    Pred <- FullModel$fitted.values
    Res <- FullModel$residuals
    ResAbs <- abs(Res)
    Y <- Pred + Res
    WF <- (lm(ResAbs ~ Pred))$fitted.values^2
    WF <- 1/WF
    Pred <- NullModel$fitted.values
    Res <- NullModel$residuals
    ResAbs <- abs(Res)
    WN <- (lm(ResAbs ~ Pred))$fitted.values^2
    WN <- 1/WN
    Xbig <- model.matrix(FullModel)
    Xsmall <- model.matrix(NullModel)
    dR2.obs <- R2W(Y, Xbig, w=WF) - R2W(Y, Xsmall, w=WN)
    
    M <- sqrt(WN)*Xsmall
    Z <- sqrt(WN)*Y
    auxNull <- lm.fit(Xsmall, Y)
    Prednull <- auxNull$fitted.values
    Resnull <- auxNull$residuals # Exchangeable residuals
    n <- length(Prednull)
    
    dR2.perm <- rep(0, B)
    for (b in 1:B)
    {
       perm <- sample(n, replace=FALSE)
       Zperm <- Resnull[perm] + Prednull # Step 3 in Freedman & Lane
       Yperm <- 1/sqrt(WN)*Zperm
       dR2.perm[b] <- R2W(Yperm, Xbig, w=WF) - R2W(Yperm, Xsmall, w=WN)
    }
    out <- list(dR2.obs=dR2.obs, dR2.perm=dR2.perm)
    return(out)
}


calcDperm <- function(NullModel, FullModel, B)
{
  dD.obs <- Dmodel(FullModel) - Dmodel(NullModel)
  Xbig <- model.matrix(FullModel)
  Xsmall <- model.matrix(NullModel)
  Y <- FullModel$model[[1]] 
  n <- length(Y)
  dD.perm <- rep(0, B)
  for (b in 1:B)
  {
    perm <- sample(n, replace=FALSE)
    Yperm <- Y[perm] # Step 3 in Freedman & Lane
    NullModelperm <- stats::glm.fit(Xsmall, Yperm, family=binomial()) 
    FullModelperm <- stats::glm.fit(Xbig, Yperm, family=binomial()) 
    dD.perm[b] <- Dmodel(FullModelperm) - Dmodel(NullModelperm)
  }
  
  out <- list(dD.obs=dD.obs, dD.perm=dD.perm)
  return(out)
}

