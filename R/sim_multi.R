#' @import stats
#' @importFrom MASS mvrnorm


sim_multi=function(n=1000, NM=3,
                   binaryOutcome = FALSE,
                   coeff, correlation=c(0.9,0.6,0.3),
                   covariables){

  Treatment=sample(0:1,n,replace=TRUE)
  data=data.frame(Intercept=1,Treatment=Treatment)


  lengthcovariables <- dim(covariables)[1]
  if(!is.null(covariables)){
  for(lc in 1:dim(covariables)[2]){
    if( length(covariables[,lc])!= n ){
      stop("Not all covariables are of length n.")
    }
    assign(names(covariables)[lc],covariables[,lc])
  }
    data=cbind(data,covariables)

  }



  nuplet=unique(data)



  ntre=c("Intercept","Treatment")
  nmed=paste("M",seq(1,NM),sep="")
  row.names(coeff)=c(nmed,"Outcome")

  if(!is.null(covariables)){
    ncov=paste("C",seq(1,dim(covariables)[2]),sep="")
    colnames(coeff)=c(ntre,ncov,nmed)
  }
  else{
    colnames(coeff)=c(ntre,nmed)
  }




Valuev=NULL
for(i in 1:dim(nuplet)[1]){
  Valuet=NULL
  for (j in 1:(dim(nuplet)[2]-1)){
    Valuet=paste(Valuet,paste(nuplet[i,-1])[j],sep="")
  }
  Valuev=c(Valuev,Valuet)
}
cf.med.names=NULL
for(nm in 1:NM){
  cf.med.names=c(cf.med.names,paste(nmed[nm],Valuev,sep="."))
}

muM=NULL
for (nm in 1:NM){
  for(i in 1:dim(nuplet)[1]){
    muM=c(muM,sum(coeff[nm,1:(ncol(coeff)-NM)]*nuplet[i,]))
  }
}

Sigma=matrix(0,length(muM),length(muM))
u=nm=1
v=dim(nuplet)[1]
for(nm in 1:NM){
  Sigma[u:v,u:v]=matrix(1,dim(nuplet)[1],dim(nuplet)[1])
  nm=nm+1
  u=v+1
  v=dim(nuplet)[1]*nm
}

a=seq(1,length(muM),by=dim(nuplet)[1])

l=1
i=1
while(i<NM){
  u=a[i]
  for(k in (i+1):(length(a))){
    v=a[k]
    w=v+dim(nuplet)[1]
    Sigma[u:(u+dim(nuplet)[1]-1),v:(w-1)]=Sigma[v:(w-1),u:(u+dim(nuplet)[1]-1)]=matrix(correlation[l],dim(nuplet)[1],dim(nuplet)[1])
    l=l+1
  }
  i=i+1
}

cf.data=MASS::mvrnorm(n=n,muM,Sigma,tol=10)
cf.data=as.data.frame(cf.data)
names(cf.data)=cf.med.names
cf.data$Treatment=Treatment



for (i in 1:dim(nuplet)[1]){
  transct=cf.data[,paste(nmed,Valuev[i],sep=".")]
  if (binaryOutcome==TRUE){
    Yerror <- stats::rlogis(n,0,1)
    cf.data[,paste("Y",Valuev[i],sep=".")]=apply(as.matrix(coeff[NM+1,])%*%t(cbind(nuplet[i,],transct)),2,sum)+Yerror>0
  }
  else{
    Yerror <- stats::rnorm(n,0,1)
    cf.data[,paste("Y",Valuev[i],sep=".")]=apply(as.matrix(coeff[NM+1,])%*%t(cbind(nuplet[i,],transct)),2,sum)+Yerror
  }

}


for(nm in 1:NM){
  data[nmed[nm]]=NA
}
data$Outcome=NA
data=data[,-1]

com=c()
for (i in 1:dim(nuplet)[1]){
  com = c(com,parse(text= paste(paste(names(nuplet)[-1],nuplet[i,-1], sep = "=="), collapse = " & ")))
}
for (i in 1:dim(nuplet)[1]){
  for(nm in 1:NM){
    data[nmed[nm]][which(eval(com[i])),]=cf.data[,paste(nmed[nm],Valuev[i],sep=".")][which(eval(com[i]))]
  }
  data$Outcome[which(eval(com[i]))]=cf.data[,paste("Y",Valuev[i],sep=".")][which(eval(com[i]))]
}

return(list(data=data,cf.data=cf.data,coeff=coeff,Sigma=Sigma))
}

