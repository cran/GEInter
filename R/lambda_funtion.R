lambda_funtion<-function(G,E,y,weight,nlambda1=20,nlambda2=20){   
  X=E
  Z=G
  Y=y
  
  X=X[weight>0,]
  Z=Z[weight>0,]
  Y=Y[weight>0]
  weight1=weight[weight>0]
  
  n<-dim(X)[1]
  q<-dim(X)[2]
  p<-dim(Z)[2]
  pp=p*(q+1)
  
  
  y=Y
  
  
  
  W=matrix(0,n,p*(q+1))
  W[,seq(from=1,to=p*(q+1),by=(q+1))]=Z
  for (i in 1:n){
    temp3=matrix(X[i,],q,1)%*%Z[i,]
    W[i,setdiff(seq(from=1,to=p*(q+1),by=1),seq(from=1,to=p*(q+1),by=q+1))]=matrix(temp3,p*q,1)
  }
  
  
  m_W=t(weight1)%*%W/sum(weight1)
  W = (W-matrix(1,n,1)%*%m_W)*(sqrt(weight1)%*%matrix(1,1,pp))
  
  
  
  group_inf=seq(from=1,to=p*(q+1),by=(q+1))
  
  
  R=list()
  for (j in 1:p){
    if (j<p){
      id=group_inf[j]:(group_inf[j+1]-1)
    }
    else {
      id=group_inf[j]:pp
    }
    a=W[,id]
    b=(t(a)%*%a*(1/n))
    R[[j]] = chol(b);
    W[,id]=t(solve(t(R[[j]]),t(a)))
  }
  
  alpha0=solve(t(X)%*%X,t(X)%*%y);
  r=y-X%*%alpha0;
  
  C_s=matrix(0,p,1);
  for (j in 1:p){
    if (j<p){
      id=group_inf[j]:(group_inf[j+1]-1);}
    else {
      id=group_inf[j]:pp;
    }
    S=t(W[,id])%*%r/n
    C_s[j]=norm(S,"2")/(sqrt(q+1))
  }
  
  
  lambda1_max=max(abs(C_s))/sqrt(1+q);
  e=(log(lambda1_max)-log(lambda1_max*0.01))/(nlambda1-1);
  lambda1_set=exp(seq(from=log(lambda1_max),to=log(lambda1_max*0.01),by=-e));
  
  lambda2_max=lambda1_max/(1+q);
  e=(log(lambda2_max)-log(lambda2_max*0.01))/(nlambda2-1);
  lambda2_set=exp(seq(from=log(lambda2_max),to=log(lambda2_max*0.01),by=-e));
  return=list(lambda1_set=lambda1_set,lambda2_set=lambda2_set)
}
