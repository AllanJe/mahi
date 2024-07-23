criterion <- function(w,treatment,mediators,outcome,covariables=NULL,weights,mu,groups){
  K=dim(mediators)[2]
  if (is.null(groups)){ groups=c(1:K) }  #default groups are each mediator on its own
  G <- max(groups)
  omega <- 0
  
  
  for (g in 1:G){
    #print(g)
    #print(length(groups==g))
    omega <- omega + norm(as.matrix( c(w$alpha[-1,groups==g],w$beta[groups==g]) ),type="2") # omega is only computed with alphas and beta
    # MODIFICATIONS omega <- omega + sqrt(sum(c(as.vector(w$alpha[,groups==g]),w$beta[groups==g])^2))
    # omega <- omega + norm(as.matrix( c(w$alpha[-1,groups==g],w$beta[groups==g]) ),type="2")
    }
  return(loss(w,treatment,mediators,outcome,covariables=covariables,weights=weights)+mu*omega)
}
