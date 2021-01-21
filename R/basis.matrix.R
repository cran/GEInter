basis.matrix<-function(x,num.knots,knots=NULL,Boundary.knots=NULL,degree=3,NorM){

  n=length(x)

  if (is.null(knots)){
    X = splines::bs(x, df=num.knots+degree+1, intercept=TRUE, degree=degree)
    knots=attr(X,'knots')
    Boundary.knots=attr(X,'Boundary.knots')
  } else {
    X = splines::bs(x, knots=knots, intercept=TRUE, degree=degree,Boundary.knots = Boundary.knots)
    knots=attr(X,'knots')
    Boundary.knots=attr(X,'Boundary.knots')
  }

  X=X[,-1]
  X=X-(matrix(1,n,1)%*%NorM)

  X=cbind(1,X)

  #X = cbind(1,X[,-1])

  return(list(X=X,knots=knots,Boundary.knots=Boundary.knots))
}
