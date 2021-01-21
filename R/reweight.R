reweight<-function(y,delta,tau){
  # calculate weights for step I in cqpcorr
  # y is the survival time
  # delta is the censoring status
  # tau is the quantile level of interest
  kmtau<-survival::survfit(survival::Surv(y,delta)~1,type="kaplan-meier")
  fhat<-(1-kmtau$surv)[order(order(y))]
  n=length(delta)
  w<-rep(1,n)
  index<-which(delta==0)
  if(length(index)>0){
    for(i in 1:length(index)){
      if(fhat[index[i]]<tau) w[index[i]]<-(tau-fhat[index[i]])/(1-fhat[index[i]])
    }}
  return(w)
}
