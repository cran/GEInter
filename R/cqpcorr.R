cqpcorr<-function(y,delta,x,z,tau,w=NULL){
  # note that x and z are univariates
  # w is the weights for step I
  # x and z are representing environmental and genetic factors respectively
  n=length(delta)
  if (is.null(w)){
    w<-reweight(y,delta,tau)
  }

  inter<-x*z

  ### step I
  #construct pseudo observations for y^\infty
  index<-which(w!=1)
  y.pse<-rep(max(y)+9999,length(index))
  x.pse<-x[index]
  z.pse<-z[index]

  yy<-c(y,y.pse)
  xx<-c(x,x.pse)
  zz<-c(z,z.pse)
  ww<-c(w,1-w[index])

  part2<-quantreg::rq(yy~xx+zz+1,weights=ww,tau)

  ### step II
  part1<-stats::lm(inter~x+z+1)

  ### step III: quantile partial corr
  ind<-as.numeric(part2$residuals[1:n]<0)
  var1<-tau-w*ind
  var2<-part1$residuals
  yxz<-(mean(var1*var2)-mean(var1)*mean(var2))/sqrt((mean(w^2)*tau-mean(w)^2*tau^2)*stats::var(var2))

  return(yxz)
}
