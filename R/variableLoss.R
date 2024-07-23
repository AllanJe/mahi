variableLoss<- function(estimation,truth,binary=FALSE){ #returns the quadratic loss for non binary variables, the logistic one for binary variables
  return((1-binary)*0.5*apply((estimation-truth)^2,2,sum)+binary*apply(-estimation*truth+log(1+exp(estimation*binary)),2,sum))
  #return((1-binary)*0.5*apply((estimation-truth)^2,2,sum)+binary*apply(truth*log(1+exp(-estimation*binary))+(1-truth)*log(1+exp(estimation*binary)),2,sum))
} 

