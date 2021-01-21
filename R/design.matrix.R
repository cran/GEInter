design.matrix<-function(x,num.knots,knots,Boundary.knots,degree,v_type,model_type='nonlinear',NorM=NULL){
  ##v_type  EC: continuous E factor, ED: decrete E factor, G: G factor
  ## x EC,ED,G: n*1, G-EC, G-ED: n*2  x[,1]: E factor x[,2]: G factor


  if (model_type=='linear'){

    if ((v_type=='EC') | (v_type=='G') | (v_type=='ED')){
      n=length(x)
      x=matrix(x,n,1)
      X=cbind(1,x)
    } else {
      X=cbind(1,x[,1]*x[,2])
    }
    num.knots=NULL
    knots=NULL
    Boundary.knots=NULL

  } else {

    if (v_type=='EC'){
      temp=basis.matrix(x,num.knots,knots,Boundary.knots,degree,NorM)
      X=temp$X
      knots=temp$knots
      Boundary.knots=temp$Boundary.knots
    } else if ((v_type=='G') | (v_type=='ED')){
      n=length(x)
      x=matrix(x,n,1)
      X=cbind(1,x)
      num.knots=NULL
      knots=NULL
      Boundary.knots=NULL
    } else if (v_type=='G-EC'){
      temp=basis.matrix(x[,1],num.knots,knots,Boundary.knots,degree,NorM)
      X=temp$X
      d_temp=dim(X)[2]
      X[,2:d_temp]=X[,2:d_temp]*(x[,2]%*%matrix(1,1,d_temp-1))
      knots=temp$knots
      Boundary.knots=temp$Boundary.knots
    } else if (v_type=='G-ED'){
      n=dim(x)[1]
      X=cbind(1,x[,1]*x[,2])
      num.knots=NULL
      knots=NULL
      Boundary.knots=NULL
    }
  }
  return(list(X=X,knots=knots,Boundary.knots=Boundary.knots,degree=degree))
}