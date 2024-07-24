\name{data_cont}
\alias{data_cont}
\docType{data}
\title{
Toy dataset
}

\description{A toy dataset for selection and causal mediation analysis with 30 uncausally related mediators.
}

\usage{data("data_cont")}
\format{
  A data frame with 100 observations and 500 canditates mediators, one exposure and an outcome.
The first 10 mediators have strong indirect effects (and in principle are easier to select), the next 10 have mild indirect effects (less easy to detect) and the next 10 have weak indirect effects (hard to detect). All other 470 variables are not mediators because either alpha_k= 0 or beta_k= 0
}

\examples{
#data(data_cont)
#
# lambda <- lambdaChoice(data_cont,name.exposure="Exposure",name.outcome=
#"Outcome",name.mediators=paste0("M", 1:500),lambdamax=30600,selectedMin=50,
#selectedMax=100)
#
# #mahifit <- mahi(data=data_cont,name.exposure="Exposure",name.outcome=
#"Outcome",name.mediators=paste0("M", 1:500),
# #lambda=lambda,Nboot=30,L0=.1,eta=2,epsilon=.001, Kmax=30,bin=FALSE,
# #dostep2=TRUE,p.adjust.method="hochberg",pvalseuil=0.05)
#
# ##The function mahi returns a list containing several results such as the
#mediators selected in the first step.
# summary(mahifit)
#
##To see the mediators selected after the step 1, use :
#   mahifit$Kmaxrank
#
## Remember that the first 30 candidates are the real mediators. The first 10
#(1 to 10) have a strong mediated effect, the next 10 (11 to 20) have a medium
#mediated effect, and the other 10 (21 to 30) have a very weak effect.
#
##To see the mediators count for each mediator after the stability selection
#(step 1), use :
#  mahifit$bootcount
#
##To see the ranking after the stability selection (step 1), use :
#mahifit$ranking
##To see the results of the multiple analysis with mediators selected at the
#step 1, use :
#summary(mahifit$multimed)
## You can see the p-value in the previous table or use :
#   mahifit$pvals
#
## The corrected p-value for the multiple test obtained with method given
#can be disp with :
#mahifit$pvalscorr
#
##It is then possible to make another correction for a multiple test with the
#mediators selected in step 1.
##To see the mediators selected after step 2 you can use these lines of codes.
#
#mahifit$step2
}
\keyword{datasets}
