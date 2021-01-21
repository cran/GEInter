predict_boosting2<-function(fit,x_new,v_type){

  y_predict=0
  id=fit$id
  v=fit$v
  for (t in 1:id){
    knots=fit$spline_result[[t]]$knots
    coef=fit$spline_result[[t]]$estimates
    degree=fit$spline_result[[t]]$degree
    Boundary.knots=fit$spline_result[[t]]$Boundary.knots
    if (is.list(x_new)){
      x_temp=x_new[[fit$variable[t]]]
    } else {
      x_temp=x_new[,fit$variable[t]]
    }
    N=design.matrix(x_temp,knots=knots,Boundary.knots=Boundary.knots,degree=degree,v_type=v_type[fit$variable[t]],model_type='nonlinear',NorM=fit$NorM)$X
    temp=N%*%coef
    y_predict=y_predict+v*temp
  }
  return(y_predict)
}