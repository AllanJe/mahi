boothdml <- function(treatment,mediators,outcome,covariables,lambda,Nboot=100,L0=.1,eta=2,weights=NULL,groups=NULL,epsilon=.001){
  
  start=Sys.time() 
  
  n <- dim(treatment)[1] 
  P <- dim(treatment)[2]
  K <- dim(mediators)[2]
  
  if (n!=dim(mediators)[1]){
    print("The number of observations between mediators and treatment does not correspond")
  }
  
  if (n!=length(outcome)){
    print("The number of observations between outcome and treatment does not correspond")
  }
  
  
  
  counts <- rep(0,K)
  wnew <- list()
  #print("Step 1: from high to low dimension")
  #pb <- txtProgressBar(min = 0, max = Nboot, style = 3)
  for(j in 1:Nboot){
    #construction des données bootstrap
    #print(paste("j=",j,"/",Nboot))
    bootsample <- sample(c(1:n),floor(n/2),replace=FALSE)
    boottreatment <- as.matrix(treatment[bootsample,])
    bootmediators <- as.matrix(mediators[bootsample,])
    if (is.null(covariables)){
      bootcovariables <- NULL
    }
    else{
      bootcovariables <- as.matrix(covariables[bootsample,])
    }
    bootoutcome <- outcome[bootsample]
    singleh = singlehdml(treatment = boottreatment, mediators = bootmediators, outcome = bootoutcome,covariables=bootcovariables, groups = groups, lambda = lambda, L0 = L0, eta = eta, weights = weights, epsilon = epsilon)
    counts <- counts + singleh$selected
    wnew[[j]] <- singleh$wnew
    #setTxtProgressBar(pb, j)
  }
  #close(pb)
  
  stop=Sys.time()
  return(list(counts = counts, time = stop-start, wnew = wnew))
}
