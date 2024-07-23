#' lambdaChoice
#'
#' 'lambdaChoice' is used to select mediators trough a high number of mediators and estimate various quantities for causal mediation analysis using "multimediate", with the selected mediators.
#'
#'@param data a data.frame with the exposure, the outcome and mediators.
#'@param name.exposure a character string indicating the name of the exposure(s) in data.
#'@param name.outcome a character string indicating the name of the outcome in data.
#'@param name.mediators a character string indicating the name of mediators in data.
#'@param name.covariables a character string indicating the name of mediators in data.
#'@param weights writing in progress
#'@param groups writing in progress
#'@param lambdamax writing in progress
#'@param N1 writing in progress
#'@param selectedMin writing in progress
#'@param selectedMax writing in progress
#'@param L0 writing in progress
#'@param eta writing in progress
#'@param epsilon writing in progress
#'@return lambdaChoice returns a list which contains the lambda to use in mahi function.
#' @export
#'
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar


lambdaChoice <- function(data,name.exposure,name.outcome,name.mediators,name.covariables=NULL,weights=NULL,groups=NULL,lambdamax=50,N1=10,selectedMin=.3*length(outcome),selectedMax=.6*length(outcome),L0=1,eta=2,epsilon=.001){

  lambdamin=0
  lambda=lambdamaxx <-lambdamax

  treatment=as.matrix(data[,name.exposure])
  mediators=as.matrix(data[,name.mediators])
  outcome=data[,name.outcome]
  if(!is.null(name.covariables)){
    covariables=as.matrix(data[,name.covariables])
  }else{
    covariables=NULL
  }
  P=length(name.exposure)
  K=length(name.mediators)
  #if (is.null(weights)){ weights$med=rep(1/P,K)
  #                      weights$out=max(1,K/100)}
  #default value for the weights is 1/P
  #if (is.null(groups)){ groups=c(1:K) }  #default groups are each mediator on its own
  if(is.null(weights)){
    weights=c()
    weights$med=rep(1/P,K)
    weights$out=1:8
  }
  NwY=length(weights$out)
  lambdaout=c()

  pb <- utils::txtProgressBar(min = 0, max = NwY, style = 3)
  for (nwy in  1:NwY){
    weightsi=weights
    weightsi$out=weights$out[[nwy]]
    lambdamin=0
    lambda=lambdamaxx


  selected <- c()
  for (i in 1:N1){
    selected <- c(selected,singlehdml(treatment,mediators,outcome,covariables=covariables,groups=groups,lambda=lambda,L0=L0,eta=eta,weights=weightsi,epsilon=epsilon)$selectedNumber)
  }



  while(mean(selected)>selectedMax){
    print("choosing a greater lambdamax")
    lambda=lambdamaxx=lambda+50
    selected <- c()
    for (i in 1:N1){
      selected <- c(selected,singlehdml(treatment,mediators,outcome,covariables=covariables,groups=groups,lambda=lambda,L0=L0,eta=eta,weights=weightsi,epsilon=epsilon)$selectedNumber)
    }
  }


    i=1       # ne pas faire plus plus de 10 divisions
    lambdamax=lambda
    test=FALSE
    while((test==FALSE) & (i<10)){

      i=i+1
      lambda=mean(c(lambdamin,lambdamax))
      selected <- c()
      for (j in 1:N1){
        selected <- c(selected,singlehdml(treatment,mediators,outcome,covariables,groups=groups,lambda=lambda,L0=L0,eta=eta,weights=weightsi,epsilon=epsilon)$selectedNumber)
      }

      if(mean(selected)<selectedMin){ lambdamax <- lambda}
      else if(mean(selected)>selectedMax){ lambdamin <- lambda}
      else{test <- TRUE}
    }



  print(paste("lambda =",lambda))
  lambdaout=c(lambdaout,lambda)
  utils::setTxtProgressBar(pb, nwy)
  }
  close(pb)
  return(lambdaout)

}
