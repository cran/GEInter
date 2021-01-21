robust.spline<-function(designmat,y,r,kmweight=NULL){



  X=designmat$X
  knots=designmat$knots
  Boundary.knots=designmat$Boundary.knots
  degree=designmat$degree

  fit<-pen.m(y, X,r,kmweight)
  return(list(estimates=fit$estimates,y_hat=fit$y_hat,knots=knots,Boundary.knots=Boundary.knots,degree=degree,X=X))
}