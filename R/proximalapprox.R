#proximal approximation for the group lasso with the group structure given by the group vector
proximalapprox <- function(u,mu,groups){
  
  G = max(groups)
  
  for(g in 1:G){
    medg <- which(groups==g)
    knorm=norm(as.matrix( c(u$alpha[-1,groups==g],u$beta[groups==g])),type="2") # norme de u|Gr the subvector of u which coordinates correspond to those of alpha|Gr and beta||Gr ,
    # MODIFICATIONS 
    coeff = ifelse((1-mu/knorm)>0, 1-mu/knorm, 0)   #coefficient par lequel multiplier les elements du groupe k
    u$alpha[,groups==g] <- u$alpha[,groups==g]*coeff
    # MODIFICATIONS u$alpha[,groups==g] <- u$alpha[,groups==g]*coeff
    u$beta[groups==g] <- u$beta[groups==g]*coeff
  }
  return(u)
}
