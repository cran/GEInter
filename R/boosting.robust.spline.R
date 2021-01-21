boosting.robust.spline<-function(x,y,loop_time,num.knots=NULL,knots=NULL,Boundary.knots=NULL,degree=1,Method,E_type,model_type="nonlinear",kmweight=NULL,v=0.1){
  y=as.matrix(y)
  x=as.matrix(x)
  q=dim(y)[2]+dim(x)[2]
  if (is.matrix(x)){
    n=dim(x)[1]
    p=dim(x)[2]
  } else {
    p=length(x)
    n=length(x[[1]])
  }


  if (!is.null(knots)){
    ss<-seq(from=1/(num.knots+1),to=num.knots/(num.knots+1),by=1/(num.knots+1))
    knots_temp=vector("list",p)
    for (ii in 1:q){
      knots_temp[[ii]]=stats::quantile(knots[[ii]],ss)
    }
    knots=knots_temp
    xx_dots=seq(from=0,to=1,by=0.0001)
    NorM = splines::bs(xx_dots, knots=knots[[1]],degree=degree,Boundary.knots = Boundary.knots[1,])
    NorM=colMeans(NorM)
  } else {
    xx_dots=seq(from=0,to=1,by=0.0001)
    NorM = splines::bs(xx_dots, df=num.knots+degree,degree=degree)
    NorM=colMeans(NorM)
  }



  u=y
  f=0

  result=list()
  variable=matrix(0,loop_time,1)
  BIC=matrix(-1e+20,loop_time,1)
  loglike=matrix(0,loop_time,1)
  t=1
  f_temp=f


  designmat=vector("list",p)
  v_type=E_type



  if (is.matrix(x)){

    if (is.null(knots)){

      for (j in 1:p){
        designmat[[j]]=design.matrix(x[,j],num.knots,knots,Boundary.knots,degree,v_type[j],model_type,NorM)
      }
    } else {
      for (j in 1:p){
        designmat[[j]]=design.matrix(x[,j],num.knots,knots[[j]],Boundary.knots[j,],degree,v_type[j],model_type,NorM)
      }
    }
  } else {
    if (is.null(knots)){
      for (j in 1:p){
        designmat[[j]]=design.matrix(x[[j]],num.knots,knots,Boundary.knots,degree,v_type[j],model_type,NorM)
      }
    } else {
      for (j in 1:p){
        designmat[[j]]=design.matrix(x[[j]],num.knots,knots[[j]],Boundary.knots[j,],degree,v_type[j],model_type,NorM)
      }
    }




  }


  hatmat_set=list()



  for (t in 1:loop_time){
    res=matrix(1e+20,p,1)
    for (j in 1:p){
      if (Method=="Robust"){
        temp=robust.spline(designmat[[j]],u,y-f,kmweight)
      } else if (Method=="LS") {
        temp=ls.spline(designmat[[j]],u,kmweight)
      }

      y_predict=temp$y_hat
      f1=y_predict+f_temp
      hatmat_set[[j]]=temp


      if (Method=="Robust"){
        ptemp=dim(temp$X)[2]
        if (ptemp==2){
          sdd1=sqrt(mean((temp$X[,2:ptemp]-mean(temp$X[,2:ptemp]))^2))
          sdd2=pcaPP::qn(u)
          RSS=n*(sdd2^2-(temp$estimate[2]*sdd1)^2)
        } else {
          xtemp=temp$X[,2:ptemp]
          xxtemp=xtemp-matrix(1,n,1)%*%colMeans(xtemp)
          robustcov=(t(xxtemp)%*%xxtemp)/n
          beta_temp=matrix(temp$estimates[2:ptemp],1,ptemp-1)
          RSS=n*(pcaPP::qn(u)^2-beta_temp%*%robustcov%*%t(beta_temp))
        }
      } else if (Method=="LS") {
        RSS=sum((y-f1)^2)
      }
      if (RSS<0){
        RSS=1e+100
      }
      df=length(unique(c(variable[1:(t-1)],j)))
      res[j]=log(RSS)+df*log(max(n,p))/n

    }

    id=which.min(res)
    temp=hatmat_set[[id]]
    y_predict=temp$y_hat
    coef=temp$estimates
    result[[t]]=temp
    variable[t]=id
    f=f+v*y_predict
    f_temp=f
    u=u-v*y_predict




    if (Method=="Robust") {
      absy=abs(y-f)
      MAD_v=stats::mad(y-f)

      c = 1.345*MAD_v
      Huber=absy^2
      Huber[absy>c]=2*c*(absy[absy>c]-c/2)
      RSS=sum(Huber)
    } else if (Method=="LS") {
      RSS=sum((y-f)^2)
    }

    if (RSS==0){
      RSS=1e-200
    }
    df=length(unique(variable[1:t]))
    BIC[t]=log(RSS) + df*log(n)/n
    if (t>1000) {
      if (abs((BIC[t]-BIC[t-1])/BIC[t-1])<1e-8){
        break
      }
    }
    t=t+1
  }

  BIC=BIC[1:t]
  variable=variable[1:t]
  loglike=loglike[1:t]

  id=which.min(BIC)
  output=list(spline_result=result,BIC=BIC,variable=variable,id=id,max_t=t,model_type=model_type,degree=degree,v=v,NorM=NorM)

  return(output)


}