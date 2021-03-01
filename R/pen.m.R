pen.m<- function(y,X,r,kmweight=NULL)
{
  n=length(y)
  rr=MASS::psi.huber(r,k=1.345*stats::mad(r))
  W=rr
  if (!is.null(kmweight)){
    W=W*kmweight
  }

  estimates=matrix(0,dim(X)[2],1)

  Xtemp=X*W
  dddd=t(Xtemp[,colSums(Xtemp)!=0])%*%X[,colSums(Xtemp)!=0]
  dde=dim(dddd)[1]
  ctemp= solve(dddd+diag(1e-5,dde,dde))%*%(t(Xtemp[,colSums(Xtemp)!=0]))
  mbeta <-ctemp%*% y
  y_hat=X[,colSums(Xtemp)!=0]%*%mbeta
  estimates[colSums(Xtemp)!=0]=mbeta

  return(list(y_hat=y_hat,estimates=estimates,scale=scale))

}
