# clear 
# rm(list=ls())
graphics.off()
library(tidyverse)

#----------------------------------------------------------#
# objective function to be minimized (for parameter fit)
#----------------------------------------------------------#

func_minimize <- function(modelfunc, param, data, prior)
{
  ret <- modelfunc(param, data, prior)
  
  # return negative log-likelihood
  return(ret$negll)
}

#----------------------------------------------------------#
# Parameter fit ML
#----------------------------------------------------------#

paramfitML <- function(modelfunctions, data, nParam)
{
  
  nModel <- length(modelfunctions)
  aic <- numeric(nModel)
  bic <- numeric(nModel)
  negll <- numeric(nModel)
  paramlist <- list()
  
  for (idxm in 1:nModel) {
    fvalmin <- Inf
    print(sprintf("Model %d:", idxm))
    
    for (idx in 1:10) {
      
      # set initial value
      initparam <- runif(nParam[idxm], 0, 1.0)
      
      res <- optim(initparam, func_minimize,
                   hessian = TRUE, modelfunc = modelfunctions[[idxm]], 
                   data=data, prior = NULL)
      
      if (res$value < fvalmin) {
        paramest <- res$par
        fvalmin <- res$value
      }
    }
    aic[idxm] <- 2*fvalmin + 2*2
    bic[idxm] <- 2*fvalmin + nParam * log(T)
    negll[idxm] <- fvalmin
    paramlist[[idxm]] <- paramest
    
    print(sprintf("Estimated value: %.2f", paramest))
    print(sprintf("log-likelihood: %.2f, AIC: %.2f, -BIC/2: %.2f", 
                  negll[idxm], 
                  aic[idxm],
                  -bic[idxm]/2))
  }
  return(list(negll = negll, aic = aic, bic = bic, paramlist = paramlist))
}

#----------------------------------------------------------#
# Parameter fit MAP
#----------------------------------------------------------#

paramfitMAP <- function(modelfunctions, data,  nParam, prior)
{
  
  nModel <- length(modelfunctions)
  lml <- numeric(nModel)
  negll <- numeric(nModel)
  paramlist <- list()
  
  for (idxm in 1:nModel) {
    fvalmin <- Inf
    print(sprintf("Model %d:", idxm))
    
    for (idx in 1:10) {
      
      # set initial value
      initparam <- runif(nParam[idxm], 0, 1.0)
      
      res <- optim(initparam, func_minimize,
                   hessian = TRUE, modelfunc = modelfunctions[[idxm]], 
                   data=data, prior = prior[[idxm]])
      
      print(res$hessian)
      
      if (res$value < fvalmin) {
        paramest <- res$par
        fvalmin <- res$value
        lp <- -fvalmin
        H <- res$hessian
      }
    }
    negll[idxm] <- fvalmin
    paramlist[[idxm]] <- paramest
    
    # log marginal likelihood (Laplace)
    lml[idxm] <- lp + nParam[idxm]/2 * log(2*pi) - 0.5 * log(det(H)) 
    
    print(sprintf("Estimated value: %.2f", paramest))
    print(sprintf("log marginal likelihood: %.2f", lml[idxm]))
    
  }
  return(list(negll = fvalmin, lml = lml, paramlist = paramlist))
}
