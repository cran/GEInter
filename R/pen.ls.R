pen.ls <- function(y, X,kmweight=NULL)
{

  if (is.null(kmweight)) {
    beta.ls=matrix(0,dim(X)[2],1)
    ttt=t(X[,colSums(X)!=0]) %*% X[,colSums(X)!=0]
    dd=dim(ttt)[1]
    temp=solve( t(X[,colSums(X)!=0]) %*% X[,colSums(X)!=0]+diag(1e-5,dd,dd))
    beta.ls[colSums(X)!=0] <- as.vector(temp %*% t(X[,colSums(X)!=0]) %*% y )
  } else {
    beta.ls=matrix(0,dim(X)[2],1)
    Xtemp=X*kmweight
    ttt=t(Xtemp[,colSums(Xtemp)!=0])%*%X[,colSums(Xtemp)!=0]
    dd=dim(ttt)[1]
    ctemp= solve(t(Xtemp[,colSums(Xtemp)!=0])%*%X[,colSums(Xtemp)!=0]+diag(1e-5,dd,dd))%*%(t(Xtemp[,colSums(Xtemp)!=0]))
    beta.ls[colSums(Xtemp)!=0] <-ctemp%*% y
  }
  y_hat=X%*%beta.ls


  return(list(beta=beta.ls,y_hat=y_hat))
}
