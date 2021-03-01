MCP_Hier<-function(G,E,y,W=NULL,WW=NULL,lambda1,lambda2,gamma1,gamma2,weight=NULL){

  n<-dim(E)[1]
  p<-dim(G)[2]
  q<-dim(E)[2]

  if (is.null(W)){
    W=array(0,dim=c(n,q,p))
    for (i in 1:n){
      temp3=matrix(E[i,],q,1)%*%G[i,]
      W[i,,]=temp3
    }
  }

  if (is.null(WW)){
    WW=matrix(0,n,p*(q+1))
    WW[,seq(from=1,to=p*(q+1),by=(q+1))]=G
    for (i in 1:n){
      temp3=matrix(E[i,],q,1)%*%G[i,]
      WW[i,setdiff(seq(from=1,to=p*(q+1),by=1),seq(from=1,to=p*(q+1),by=q+1))]=matrix(temp3,p*q,1)
    }


  }


  if (!is.null(weight)){

    nonzero_id=which(weight!=0)

    y=y[weight!=0]
    E=E[weight!=0,]
    G=G[weight!=0,]

    W=W[weight!=0,,]

    weight=weight[weight!=0]

    n=length(y)

    weight=weight/sum(weight)*n


    ymean=sum(weight*y)/sum(weight)

    gmean=colSums(G*(weight%*%matrix(1,1,p)))/sum(weight)

    emean=colSums(E*(weight%*%matrix(1,1,q)))/sum(weight)

    wmean=matrix(0,q,p)
    for (k in 1:q){
      wmean[k,]=colSums(W[,k,]*(weight%*%matrix(1,1,p)))/sum(weight)
    }




    y=sqrt(weight)*(y-ymean)

    G=(sqrt(weight)%*%matrix(1,1,p))*(G-matrix(1,n,1)%*%gmean)


    E=(sqrt(weight)%*%matrix(1,1,q))*(E-matrix(1,n,1)%*%emean)


    for (k in 1:q){
      W[,k,]=(sqrt(weight)%*%matrix(1,1,p))*(W[,k,]-matrix(1,n,1)%*%wmean[k,])
    }

    WW=matrix(0,n,p*(q+1))
    WW[,seq(from=1,to=p*(q+1),by=(q+1))]=G
    for (i in 1:n){
      WW[i,setdiff(seq(from=1,to=p*(q+1),by=1),seq(from=1,to=p*(q+1),by=q+1))]=matrix(W[i,,],p*q,1)
    }

  } else {
    nonzero_id=1:n
  }



  beta0=matrix(0,p,1)
  theta0=matrix(0,q,p)
  E_temp=cbind(1,E)


  if (!is.null(weight)){

    zero_id=matrix(0,q,1)
    for (k in 1:q){
      zero_id[k]=length(unique(E[,k]))
    }
    zero_id=which(zero_id!=1)
    alpha0=matrix(0,q,1)
    tempnn=t(E[,zero_id])%*%E[,zero_id]
    ddtemp=dim(tempnn)[1]
    tempnn=tempnn+diag(1e-5,ddtemp,ddtemp)
    alpha0[zero_id]=solve(tempnn)%*%t(E[,zero_id])%*%y
    r=y-E%*%alpha0
  } else {
    zero_id=matrix(0,q+1,1)
    for (k in 1:q){
      zero_id[k+1]=length(unique(E_temp[,k+1]))
    }
    zero_id=which(zero_id!=1)
    alpha_c0=matrix(0,q+1,1)
    if (rcond(t(E_temp[,zero_id])%*%E_temp[,zero_id])>1e-16){
      temppp=t(E_temp[,zero_id])%*%E_temp[,zero_id]
      dddim=dim(temppp)[1]
      temppp=temppp+diag(1e-5,dddim,dddim)
      alpha_c0[zero_id]=solve(temppp)%*%t(E_temp[,zero_id])%*%y
    }
    alpha0=alpha_c0[2:(q+1)]
    intercept0=alpha_c0[1]
    r=y-E_temp%*%alpha_c0

  }



  # alpha0=solve(t(E)%*%E)%*%t(E)%*%y
  # r=y-E%*%alpha0

  loop_time=1
  diff_v=1
  objective=mean(r^2)/2


  beta1=beta0
  theta1=theta0




  while ((diff_v>1e-4) && (loop_time<50)){

    #para0=c(alpha0,beta0,c(theta0))

    eta=G
    for (k in 1:q){
      eta=eta+W[,k,]*(matrix(1,n,1)%*%theta0[k,])
    }

    for (j in 1:p){

      eta_j=eta[,j]

      v_j=mean(eta_j^2)

      u_j=matrix(r,1,n)%*%eta_j/n+v_j*beta0[j]

      temp_b=u_j/(v_j+1e-20)

      # if ((v_j-1/gamma1)<0){
      #   ttt[j]=1
      # }


      if (abs(temp_b)>gamma1*lambda1){
        beta1[j]=temp_b
      } else {
        c_j=u_j
        if (abs(c_j)>lambda1){
          beta1[j]=(c_j-sign(c_j)*lambda1)/(v_j-1/gamma1)
        } else {
          beta1[j]=0
        }
      }

      if (abs(beta1[j]-beta0[j])>1e-5){
        r=r-matrix(eta_j,n,1)%*%(beta1[j]-beta0[j])
        beta0=beta1
      }
    }


    active_id=which(beta0!=0)
    for (k in 1:q){

      for (j in active_id){

        eta_j=beta0[j]*W[,k,j]
        v_j=mean(eta_j^2)
        u_j=matrix(r,1,n)%*%eta_j/n+v_j*theta0[k,j] ###############

        temp_b=u_j/(v_j+1e-20)

        #print(v_j)



        if (is.finite(temp_b)){

          if ((v_j-1/gamma2)<0){
            theta1[k,j]=0
          } else {
            if (abs(temp_b)>gamma2*lambda2){
              theta1[k,j]=temp_b
            } else {
              c_j=u_j
              if (abs(c_j)>lambda2){
                theta1[k,j]=(c_j-sign(c_j)*lambda2)/(v_j-1/gamma2)
              } else {
                theta1[k,j]=0
              }
            }
          }
        } else {
          theta1[k,j]=0
        }


        if (abs(theta1[k,j]-theta0[k,j])>0){
          r=r-eta_j*(theta1[k,j]-theta0[k,j])
          theta0=theta1
        }
        #print(max(abs(r)))
      }
    }


    if (!is.null(weight)){

      temp=r+E%*%alpha0
      zero_id=matrix(0,q,1)
      for (k in 1:q){
        zero_id[k]=length(unique(E[,k]))
      }
      zero_id=which(zero_id!=1)
      alpha1=matrix(0,q,1)

      tempnn=t(E[,zero_id])%*%E[,zero_id]
      ddtemp=dim(tempnn)[1]
      tempnn=tempnn+diag(1e-5,ddtemp,ddtemp)

      alpha1[zero_id]=solve(tempnn)%*%t(E[,zero_id])%*%temp
      r=r-E%*%(alpha1-alpha0)
      alpha0=alpha1
    } else {

      temp=r+E_temp%*%alpha_c0
      zero_id=matrix(0,q+1,1)
      for (k in 1:q){
        zero_id[k+1]=length(unique(E_temp[,k+1]))
      }
      zero_id=which(zero_id!=1)
      alpha_c1=matrix(0,q+1,1)
      if (rcond(t(E_temp[,zero_id])%*%E_temp[,zero_id])>1e-16){
        tempnn=t(E_temp[,zero_id])%*%E_temp[,zero_id]
        ddtemp=dim(tempnn)[1]
        tempnn=tempnn+diag(1e-5,ddtemp,ddtemp)

        alpha_c1[zero_id]=solve(tempnn)%*%t(E_temp[,zero_id])%*%temp
      } else{
        alpha_c1=alpha_c0
      }

      #alpha_c1=solve(t(E_temp)%*%E_temp)%*%t(E_temp)%*%temp
      r=r-E_temp%*%(alpha_c1-alpha_c0)
      alpha_c0=alpha_c1
      alpha0=alpha_c1[2:(q+1)]
      intercept0=alpha_c1[1]
    }



    # temp=r+E%*%alpha0
    # alpha1=solve(t(E)%*%E)%*%t(E)%*%temp
    # r=r-E%*%(alpha1-alpha0)
    # alpha0=alpha1


    #para1=c(alpha1,beta1,c(theta1))

    temp1=pmin(abs(beta1),lambda1*gamma1)
    MCP_pen1=sum(temp1-temp1^2/(2*lambda1*gamma1))

    temp2=pmin(abs(theta1),lambda2*gamma2)
    MCP_pen2=sum(temp2-temp2^2/(2*lambda2*gamma2))


    objective1=mean(r^2)/2+lambda1*MCP_pen1+lambda2*MCP_pen2

    diff_v=abs(objective1-objective)/abs(objective)

    #diff_v=norm(para1-para0,'2')/norm(para0,'2')

    objective=objective1
    loop_time=loop_time+1
    # print(diff_v)
    #
    # print('-----------------------')

  }

  RSS=sum(r^2)

  beta1[abs(beta1)<1e-4]=0
  beta_temp=(matrix(1,q,1)%*%t(beta1))*theta1

  beta_temp[abs(beta_temp)<1e-4]=0

  beta_temp=rbind(t(beta1),beta_temp)

  df=sum(beta_temp!=0)




  b_estimate=matrix(beta_temp,p*(q+1),1)
  rr=y-matrix(WW[,b_estimate!=0],n,df)%*%b_estimate[b_estimate!=0]

  if (!is.null(weight)){
    zero_id=matrix(0,q,1)
    for (k in 1:q){
      zero_id[k]=length(unique(E[,k]))
    }
    zero_id=which(zero_id!=1)
    alpha0=matrix(0,q,1)

    tempnn=t(E[,zero_id])%*%E[,zero_id]
    ddtemp=dim(tempnn)[1]
    tempnn=tempnn+diag(1e-5,ddtemp,ddtemp)

    alpha0[zero_id]=solve(tempnn)%*%t(E[,zero_id])%*%rr
  } else {
    zero_id=matrix(0,q+1,1)
    for (k in 1:q){
      zero_id[k+1]=length(unique(E_temp[,k+1]))
    }
    zero_id=which(zero_id!=1)
    alpha_c1=matrix(0,q+1,1)
    if (rcond(t(E_temp[,zero_id])%*%E_temp[,zero_id])>1e-16){

      tempnn=t(E_temp[,zero_id])%*%E_temp[,zero_id]
      ddtemp=dim(tempnn)[1]
      tempnn=tempnn+diag(1e-5,ddtemp,ddtemp)

      alpha_c1[zero_id]=solve(tempnn)%*%t(E_temp[,zero_id])%*%rr
    }
    alpha0=alpha_c1[2:(q+1)]
    intercept0=alpha_c1[1]

  }


  if (!is.null(weight)){
    intercept0=(ymean-sum(emean*alpha0)-sum(gmean*beta_temp[1,]))
    for (k in 1:q){
      intercept0=intercept0-sum(wmean[k,]*beta_temp[k+1,])
    }
  }

  #intercept0=0

  result=list(intercept=intercept0,alpha=alpha0,beta=beta_temp,objective_set=objective,RSS=RSS,df=df)

  return(result)
}








