ls.spline<-function(designmat,y,kmweight=NULL){


  X=designmat$X
  knots=designmat$knots
  Boundary.knots=designmat$Boundary.knots
  degree=designmat$degree

  fit<-pen.ls(y, X,kmweight)

  return(list(estimates=fit$beta,y_hat=fit$y_hat,knots=knots,Boundary.knots=Boundary.knots,degree=degree,X=X))
}