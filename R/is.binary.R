# is.binary <- function(x){
#   return(length(unique(x))==2)
# }

is.binary=function(mediators){

  if(is.null(dim(mediators)[2]) | dim(mediators)[2]==1){
    out=is.2(length(unique(mediators)))
  }
  else{
    out=unlist(lapply(lapply(lapply(X=as.data.frame(mediators),2,FUN=unique),length),is.2))
  }
  return(out)
}
