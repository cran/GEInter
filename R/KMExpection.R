KMExpection<-function(y,delta){
  fit<-survival::survfit(survival::Surv(time=y,event=delta) ~ 1)
  St<-fit$surv
  step_point<-fit$time
  n=length(y)
  y_estimate<-y
  for (i in 1:n){
    if (delta[i]==0){
      temp1=sum((step_point>=y[i])*(St-c(St[-1],0)))+1e-10
      temp2=sum(step_point*(step_point>=y[i])*(St-c(St[-1],0)))
      y_estimate[i]<-temp2/temp1
    }
  }

  return(y_estimate=y_estimate)

}
