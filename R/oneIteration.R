oneIteration <- function(w,treatment,mediators,outcome,covariables=NULL,weights,groups,lambda,L0,eta){
  
  L <- L0
  
  u <- c()
  wgrad <- lossgrad(w,treatment=treatment,mediators=mediators,outcome=outcome,covariables=covariables,weights=weights)
  u$alpha <- w$alpha - 1/L*wgrad$alpha 
  u$beta <- w$beta - 1/L*wgrad$beta
  u$gamma <- w$gamma - 1/L*wgrad$gamma
  if(!is.null(covariables)){
    u$xi <- w$xi - 1/L*wgrad$xi
    u$psi <- w$psi - 1/L*wgrad$psi
  }
  
  wL <- proximalapprox(u,mu=lambda/L,groups=groups)
  # test=f(v)-f(v+1)+<Gradf(v),v+1 - v> +L/2 ||v+1 -v||22
  #avec wL= v+1 et w=v
  test <-  loss(w,treatment,mediators,outcome,covariables,weights=weights) - loss(wL,treatment,mediators,outcome,covariables,weights=weights) +
    sum(wgrad$alpha*(wL$alpha-w$alpha))+ # MODIFICATIONS sum(wgrad$alpha*(wL$alpha-w$alpha))+
    sum(wgrad$beta*(wL$beta-w$beta))+
    sum(wgrad$gamma*(wL$gamma-w$gamma))+
    L/2*norm((wL$alpha-w$alpha),type="2")^2 +  # MODIFICATIONS L/2*norm((wL$alpha-w$alpha),type="F")^2 +
    L/2*norm(as.matrix(wL$beta-w$beta),type="2")^2+ # MODIFICATIONS L/2*norm(as.matrix(wL$beta-w$beta),type="F")^2+
    L/2*norm(as.matrix(wL$gamma-w$gamma),type="2")^2 # MODIFICATIONS L/2*norm(as.matrix(wL$gamma-w$gamma),type="F")^2
    if(!is.null(covariables)){
      test = test + 
        sum(wgrad$xi*(wL$xi-w$xi))+
        sum(wgrad$psi*(wL$psi-w$psi))+ 
        L/2*norm((wL$xi-w$xi),type="2")^2 + 
        L/2*norm(as.matrix(wL$psi-w$psi),type="2")^2
    }
    # increase L until test<0
  while(test<0 & L<10000){
    #print(paste("L=",L))
    L <-eta*L
    
    
    wgrad <- lossgrad(w,treatment=treatment,mediators=mediators,outcome=outcome,covariables=covariables,weights=weights)
    u$alpha <- w$alpha - 1/L*wgrad$alpha 
    u$beta <- w$beta - 1/L*wgrad$beta
    u$gamma <- w$gamma - 1/L*wgrad$gamma
    if(!is.null(covariables)){
      u$xi <- w$xi - 1/L*wgrad$xi
      u$psi <- w$psi - 1/L*wgrad$psi
    }
    wL <- proximalapprox(u,mu=lambda/L,groups=groups)
    
    test <- loss(w,treatment,mediators,outcome,covariables=covariables,weights=weights)- loss(wL,treatment,mediators,outcome,covariables=covariables,weights=weights) +
      sum(wgrad$alpha*(wL$alpha-w$alpha))+ # MODIFICATIONS sum(wgrad$alpha*(wL$alpha-w$alpha))+
      sum(wgrad$beta*(wL$beta-w$beta))+
      sum(wgrad$gamma*(wL$gamma-w$gamma))+
      L/2*norm((wL$alpha-w$alpha),type="2")^2 +  # MODIFICATIONS L/2*norm((wL$alpha-w$alpha),type="F")^2 +
      L/2*norm(as.matrix(wL$beta-w$beta),type="2")^2+ # MODIFICATIONS L/2*norm(as.matrix(wL$beta-w$beta),type="F")^2+
      L/2*norm(as.matrix(wL$gamma-w$gamma),type="2")^2 # MODIFICATIONS L/2*norm(as.matrix(wL$gamma-w$gamma),type="F")^2
    if(!is.null(covariables)){
      test = test + 
        sum(wgrad$xi*(wL$xi-w$xi))+
        sum(wgrad$psi*(wL$psi-w$psi))+ 
        L/2*norm((wL$xi-w$xi),type="2")^2 + 
        L/2*norm(as.matrix(wL$psi-w$psi),type="2")^2
    }
    
  }
  
  return(wL)
  
}
