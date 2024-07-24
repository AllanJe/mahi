#' mahi
#'
#'
#' 'mahi' is used to select mediators trough a high number of mediators and estimate various quantities for causal mediation analysis using "multimediate", with the selected mediators.
#'
#'@param data a data.frame with the exposure, the outcome and mediators.
#'@param name.exposure a character string indicating the name of the exposure in data.
#'@param name.outcome a character string indicating the name of the outcome in data.
#'@param name.mediators a character string indicating the name of mediators in data.
#'@param name.covariables a character string indicating the name of mediators in data.
#'@param lambda vector or a single value of lambda use for group-lasso.
#'@param Nboot number of bootstrap for stability selection of mediators. Default is 100.
#'@param L0 writing in progress
#'@param eta writing in progress
#'@param weights writing in progress
#'@param groups writing in progress
#'@param epsilon writing in progress
#'@param Kmax maximum of mediators keept after the ranking obtained with stability selection.
#'@param p.adjust.method a character string indicating the name of the method for the correction for the multiple test. See help with p.adjust.methods.
#'@param pvalseuil the p-value for the multiple test
#'@param bin a logical value indicating if the outcome is binary. if 'TRUE' a probit regression will be use is the second step. Default is 'FALSE'.
#'@param dostep2 a logical value indicating if the a multiple test using the p-value calculated with multimediate (a multiple mediation analysis) have to be made. Default is 'FALSE'.
#'@param p.adjust.method correction method to use on the p-value calculated with multimediate among c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none"). See p.adjust help for more details. Default is "bonferroni".
#'@param pvalseuil p-value threshold for multiple testing. Default is 0.05.
#'@return mahi returns an object of class "mahi", a list that contains the components listed below.
#' @export
#'
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @import stats
#' @import multimediate
#'
#' @examples
#' # example code
#'
#' #data(data_cont)
#'
#' #lambda <- lambdaChoice(data_cont,name.exposure="Exposure",name.outcome=
#' #"Outcome",name.mediators=paste0("M", 1:500),lambdamax=30600,selectedMin=50,
#' #selectedMax=100)
#'
#' #mahifit <- mahi(data=data_cont,name.exposure="Exposure",name.outcome=
#' #"Outcome",name.mediators=paste0("M", 1:500),
#' #lambda=lambda,Nboot=30,L0=.1,eta=2,epsilon=.001, Kmax=30,bin=FALSE,
#' #dostep2=TRUE,p.adjust.method="hochberg",pvalseuil=0.05)
#'
#' ##The function mahi returns a list containing several results such as the
#' ## mediators selected in the first step.
#' #summary(mahifit)
#'
#' ##To see the mediators selected after the step 1, use :
#' #   mahifit$Kmaxrank
#'
#' ## Remember that the first 30 candidates are the real mediators. The first
#' ## 10 (1 to 10) have a strong mediated effect, the next 10 (11 to 20) have a
#' ## medium mediated effect, and the other 10 (21 to 30) have a very
#' ## weak effect.
#'
#' ##To see the mediators count for each mediator after the stability selection
#' ## (step 1), use :
#' #  mahifit$bootcount
#'
#' ##To see the ranking after the stability selection (step 1), use :
#' #mahifit$ranking
#'
#' ##To see the results of the multiple analysis with mediators selected at
#' ## the step 1, use :
#' #summary(mahifit$multimed)
#'
#' ## You can see the p-value in the previous table or use :
#' #   mahifit$pvals
#'
#' ## The corrected p-value for the multiple test obtained with method given
#' ## can be disp with :
#' #mahifit$pvalscorr
#'
#' ##It is then possible to make another correction for a multiple test with the
#' ## mediators selected in step 1.
#' ##To see the mediators selected after step 2 you can use these lines of
#' # codes.
#'
#' #mahifit$step2




mahi=function(data,name.exposure,name.outcome,name.mediators,name.covariables=NULL,lambda=NULL,Nboot=100,L0=.1,eta=2,weights=NULL,groups=NULL,epsilon=.001,Kmax=NULL,bin=FALSE,dostep2=FALSE,p.adjust.method="bonferroni",pvalseuil=0.05){
  n=dim(data)[1]
  K=length(name.mediators)
  P=length(name.exposure)
  if(is.null(Kmax)){
    Kmax=2*n/log(n)
  }

  #if (is.null(weights)){
  #  weights$med=rep(1/P,K)
  #  weights$out=max(1,K/100)
  #}
  if(is.null(weights)){
    weights=c()
    weights$med=rep(1/P,K)
    weights$out=1
  }

  if(is.null(lambda)){
    lambda = lambdaChoice(data=data,name.exposure=name.exposure,name.outcome=name.outcome,name.mediators=name.mediators,
                          name.covariables=name.covariables,
                          weights=weights,groups=groups,lambdamax=20,selectedMin=10,selectedMax=50,L0=L0,eta=eta,epsilon=epsilon)
  }
  treatment=as.matrix(data[,name.exposure])
  mediators=as.matrix(data[,name.mediators])
  outcome=data[,name.outcome]
  if(is.null(name.covariables)){
    covariables=NULL
  }
  else{
    covariables=as.matrix(data[,name.covariables])
  }
  print("Step 1: from high to low dimension")
  bootstep=ranking.id=ranking=Kmaxrank=list()
  #select_seuil=NULL
  NwY=length(weights$out)
  #seuils = Nboot*0.8#max(Nboot*0.8,qbinom((1/(100*NwY)), Nboot, selectedMax/K, lower.tail = FALSE, log.p = FALSE))
  pb <- utils::txtProgressBar(min = 0, max = NwY, style = 3)
  for (i in 1:NwY ) {
    weightsi=weights
    weightsi$out=weights$out[[i]]
    bootstep[[paste("wY =",i)]] <- boothdml(treatment,mediators,outcome,covariables,lambda=lambda[i],Nboot,L0,eta,weights=weightsi,groups=NULL,epsilon)
    #select_seuil=c(select_seuil,name.mediators[bootstep[[paste("wY =",i)]]$counts>=seuils])
    if ( i == 1){
      bootcount <-bootstep[[paste("wY =",i)]]$counts
    }
    bootcount <- rbind(bootcount,bootstep[[paste("wY =",i)]]$counts)
    # ranking.id[[paste("wY =",i)]]=order(bootstep[[paste("wY =",i)]]$counts,decreasing = TRUE)
    # ranking[[paste("wY =",i)]]=name.mediators[ranking.id[[paste("wY =",i)]]]
    # Kmaxrank[[paste("wY =",i)]]= ranking[1:Kmax]
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  #select_seuil=unique(select_seuil)

  maxcount = apply(bootcount,2,max)
  ranking.id=order(maxcount ,decreasing = TRUE)
  ranking=name.mediators[ranking.id]
  Kmaxrank= ranking[1:Kmax]

  pvaleur=stats::pbinom(maxcount, Nboot, 0.25, lower.tail = FALSE, log.p = FALSE)

  step1=data.frame(Mediators=name.mediators,Pvalue=as.numeric(pvaleur),Maxcount=maxcount)

  step1=step1[order(step1$Pvalue,decreasing = FALSE),]
  if(!dostep2){
    out=list(step1=step1,ranking=ranking,ranking.id=ranking.id,bootcount=bootcount,bootstep=bootstep,
             n=n,K=K,Kmax=Kmax,lambda=lambda)
  }
  else{ print("Step 2: Multiple test on indirect effects")
    step2=rep(FALSE,K)
    if(P>1){
      step2byt=array(FALSE,dim=c(K,P))

    }
    multimed =pvals=pvalscorr=list()
    pvalsbytreat=matrix(NA,Kmax,P)
    for(p in 1:P){
      lmodel.m=list()
      medforout=""
      for(medch in Kmaxrank){
        lmodel.m[[medch]] = stats::lm(stats::formula(paste( medch,"~",name.exposure[p])), data = data)
        medforout=paste(medforout,medch,sep=" + ")
      }
      formout1 = paste(name.outcome,"~",name.exposure[p])
      if(bin){
        model.y = stats::glm(stats::formula( paste(formout1,medforout,sep=" ")), data = data, family=stats::binomial(link = "logit"))
      }
      else{model.y = stats::lm(stats::formula( paste(formout1,medforout,sep=" ")), data = data)}


      correlated=TRUE
      if(length(Kmaxrank)==1){correlated=FALSE}
      multimed[[p]]=multimediate::multimediate(lmodel.m,correlated,model.y,treat=name.exposure[p])
      pvalsbytreat[,p]=summary(multimed[[p]],opt="avg")[seq(3,length(Kmaxrank)*2+2,2),5]
    }

    pvalscorr=apply(pvalsbytreat,2,stats::p.adjust,method = p.adjust.method)
    selmed1=pvalscorr<=pvalseuil
    if(P>1){
      for (p in 1:P){
        selmedbyt=Kmaxrank[selmed1[,p]]
        step2byt[which(name.mediators %in% selmedbyt),p]=TRUE
      }}

    selmed2=apply(selmed1,1,sum)
    selmed3=Kmaxrank[selmed2==P]
    step2[which(name.mediators %in% selmed3)]=TRUE




    if(P>1){
      out=list(step2byt=step2byt,step2=step2,
               step1=step1,ranking=ranking,ranking.id=ranking.id,bootcount=bootcount,
               multimed=multimed,pvals=pvalsbytreat,pvalscorr=pvalscorr,
               n=n,K=K,Kmax=Kmax,lambda=lambda,bin=bin,step2=step2)
    }
    else{
      out=list(step2=step2,
               step1=step1,ranking=ranking,ranking.id=ranking.id,bootcount=bootcount,
               multimed=multimed[[1]],pvals=pvalsbytreat,pvalscorr=pvalscorr,
               n=n,K=K,Kmax=Kmax,lambda=lambda,bin=bin,step2=step2)
    }}

  class(out) <- "mahi"
  return(out)
}
