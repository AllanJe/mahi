#function f
loss <- function(w,treatment,mediators,outcome,covariables=NULL,weights){
  
  K=dim(mediators)[2]
  P=dim(treatment) [2]
  binarymed <- is.binary(mediators)
  outcome<- as.matrix(outcome)
  if(is.null(covariables)){
    mhat <- treatment %*% w$alpha 
    zhat = treatment %*% w$gamma + mediators %*% w$beta 
  }
  else{
    mhat <- treatment %*% w$alpha + covariables %*% w$xi
    zhat = treatment %*% w$gamma + mediators %*% w$beta + covariables %*% w$psi
  }
  
  
  if (is.null(weights)){ 
    weights$med=rep(1/P,K) #default value for the weights is 1/P
    weights$out=max(1,K/100)}
  
  loss<- sum(variableLoss(estimation=mhat,truth=mediators,binary=binarymed)*(weights$med)) +
    weights$out*variableLoss(estimation=zhat,truth=outcome,binary=is.2(length(unique(outcome))))
  
  return(loss/length(outcome))
}
