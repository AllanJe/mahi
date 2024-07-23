# called in function loss (f)
lossgrad <- function(w,treatment,mediators,outcome,covariables=NULL,weights){

  binarymed <- is.binary(mediators)
  if(!is.null(covariables)){
    mhat <- treatment%*%w$alpha + covariables %*% w$xi
    outcomehat = treatment %*% w$gamma + mediators %*% w$beta + covariables %*%w$psi
  }
  else{
    mhat <- treatment%*%w$alpha
    outcomehat = treatment %*% w$gamma + mediators %*% w$beta
  }
  
  for (k in 1:dim(mediators)[2]){
    if (is.binary(as.matrix(mediators[,k]))){
      mhat[,k] <- 1/(1+exp(-mhat[,k]))
    }
  }

  if (is.2(length(unique(outcome)))){
    outcomehat <- 1/(1+exp(-outcomehat))
  }

  
  P=dim(treatment)[2]
  if (is.null(weights)){
    weights$med=rep(1/P,dim(mediators)[2]) #default value for the weights is 1/P
    weights$out=max(1,dim(mediators)[2]/100)}

  wgrad <- c()
  wgrad$alpha <- t(treatment) %*% (mhat-mediators) %*% diag(weights$med) /length(outcome)
  wgrad$beta <- weights$out*as.vector(t(mediators) %*% (outcomehat-outcome)) /length(outcome)
  wgrad$gamma <- weights$out*as.vector(t(treatment) %*% (outcomehat-outcome))  /length(outcome)
  
  if(!is.null(covariables)){
    wgrad$psi <- weights$out*as.vector(t(covariables) %*% (outcomehat-outcome))  /length(outcome)
    wgrad$xi <- t(covariables) %*% (mhat-mediators) %*% diag(weights$med) /length(outcome)
  }

  return(wgrad)
}
