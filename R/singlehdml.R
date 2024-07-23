singlehdml <- function(treatment,mediators,outcome,covariables=NULL,lambda,L0=1,eta=2,weights=NULL,groups=NULL,epsilon=.001){

  n <- dim(treatment)[1]
  P <- dim(treatment)[2]
  K <- dim(mediators)[2]


  if (is.null(groups)){ groups=c(1:K) }  #default groups are each mediator on its own

  if (n!=dim(mediators)[1]){
    print("The number of observations between mediators and treatment does not correspond")
  }

  if (n!=length(outcome)){
    print("The number of observations between outcome and treatment does not correspond")
  }

  # ajout d'une colonne de 1 au début de T qui permet de prendre en compte les ordonnées à l'origine dans les deux combinaisons linéaires
  treatment <- cbind(rep(1,n),treatment)

  # initialisation de alpha, beta, gamma
  w <- c()
  w$alpha <- matrix(stats::runif(K*P+K,-1,1),nrow = P+1,ncol=K)   #la première ligne contient les ordonnées à l'origine pour les Mk, \alpha_{pk} correspond à alpha[p+1,k]
  w$beta <- stats::runif(K,-1,1)
  w$gamma <- stats::runif(P+1,-1,1)               #la première ligne contient l'ordonnée à l'origine pour Y, \gamma_{p} correspond à gamma[p+1]
  if(!is.null(covariables)){
    Q <- dim(covariables)[2]
    w$xi = matrix(stats::runif(K*Q,-1,1),nrow = Q,ncol=K)
    w$psi= stats::runif(Q,-1,1)
  }


  # vecteurs des groupes pour retrouver les coefficients à prendre en compte dans la pénalité
  # pour alpha, la k^e colonne est dans le groupe k sauf l'ordonnée à l'origine qui est dans 0, pour beta, k^e coord dans k, tout gamma dans 0
  #groups <- c()
  #groups$alpha <- rbind(rep(0,K),matrix(c(1:K),nrow=P,ncol=K,byrow=TRUE))
  #groups$beta <- c(1:K)
  #groups$gamma <- rep(0,P+1)

  #print("start")
  #print(criterion(w,treatment,mediators,outcome,weights=weights,mu=lambda,groups=groups))
  wnew <- oneIteration(w,treatment,mediators,outcome,covariables=covariables,groups=groups,lambda=lambda,L0=L0,eta=eta,weights=weights)


  i=1
  while(criterion(w,treatment,mediators,outcome,covariables,weights=weights,mu=lambda,groups=groups)-criterion(wnew,treatment,mediators,outcome,covariables,weights=weights,mu=lambda,groups=groups)>epsilon & i<10){
    i=i+1
    #print(paste("i=",i))
    #print(criterion(w,treatment,mediators,outcome,weights=weights,mu=lambda,groups=groups))
    w <- wnew
    wnew <- oneIteration(w,treatment,mediators,outcome,covariables,groups=groups,lambda=lambda,L0=L0,eta=eta,weights = weights)
  }

  selectedNumber <- 0
  selected <- rep(0,K)
  for(k in 1:K){
    knorm=wnew$alpha[-1,k]*wnew$beta[k]
    selected[k]=(sum(knorm==0)==0)*1
    selectedNumber <- selectedNumber + selected[k]
    # MODIFICATIONS knorm=norm(as.matrix( c(wnew$alpha[,k],wnew$beta[k])),type="F")
    # knorm=norm(as.matrix( c(wnew$alpha[-1,k],wnew$beta[k])),type="2")
    # selected[k] <- (knorm>0)*1
    #selectedNumber <- selectedNumber + (knorm>0)
  }
  #print(paste("selectedNumber =",selectedNumber))
  return(list(wnew=wnew,selectedNumber=selectedNumber,selected=selected,criterion=criterion(wnew,treatment,mediators,outcome,covariables,weights=weights,mu=lambda,groups=groups)))
}
