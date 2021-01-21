ptr1<-function(G,E,Y,lambda1,lambda2,gamma1=6,gamma2=6,max_init,h=NULL,tau=0.4,mu=2.5,family=c('continuous','survival'),W,WW){

  family = match.arg(family)
  subsamp=NULL
  y=as.matrix(Y)#y is n\times2 matrix if survival, n\times 1 otherwise
  ##survival
  if (family=='survival'){
    delta=y[,2]
    y=y[,1]
    weight=kmw(y,delta)
  } else {
    weight=NULL
  }



  n=dim(E)[1]
  q=dim(E)[2]
  p=dim(G)[2]


  pp=p*(q+1)




  max_iter=20

  GE=cbind(E,WW)
  GE=as.matrix(GE)

  if (!is.null(weight)){

    weight_new=weight/sum(weight)*sum(weight!=0)
    #weight_new=as.matrix(weight_new,n,1)
    ymean=sum(weight_new*y)/sum(weight_new)
    tempp=weight_new%*%matrix(1,1,(p+q+p*q))
    ttemp=GE*tempp
    gmean=base::colSums(ttemp)/sum(weight_new)
    y_new=sqrt(weight_new)*(y-ymean)
    GE_new=(sqrt(weight_new)%*%matrix(1,1,(p+q+p*q)))*(GE-matrix(1,n,1)%*%gmean)
    nonzero_id=which(weight!=0)
  } else {
    y_new=y
    GE_new=GE
    nonzero_id=1:n
  }

  if (is.null(subsamp)){
    obj_init=matrix(1e+20,max_init,1)
    subsamp=matrix(0,max_init,n)

    for (ii in 1:max_init){
      #nsub=sample(n)[1:(q+10)]
      if (!is.null(weight)){
        nsub=sample(length(nonzero_id))[1:(q+10)]
        nsub=nonzero_id[nsub]
      } else {
        nsub=sample(n)[1:(q+10)]
      }
      tt=1
      diff_c=1
      objective=100
      while ((tt<=max_iter) && (diff_c>1e-5)){
        nsub1=nsub

        # if (tt==1){
        #   result_temp<-ncvreg(GE_new[nsub,],y_new[nsub],family='gaussian',penalty='MCP',gamma1)
        #   df_set=base::colSums(result_temp$beta!=0)
        #   BIC=log(result_temp$loss)+df_set*log(n)/n
        #
        #   idx=which.min(BIC)
        #   temp=result_temp$beta[,idx]
        #   if (!is.null(weight)){
        #     temp[1]=(ymean-sum(gmean*temp[2:(p+q+p*q+1)]))
        #   }
        #   id=temp!=0
        #   temp2=cbind(1,GE)
        #   res_temp=y-temp2[,id]%*%temp[id]
        #   objective1=200
        # } else {

        result_temp=MCP_Hier(G[nsub,],E[nsub,],y[nsub],W[nsub,,],WW[nsub,],lambda1,lambda2,gamma1,gamma2,weight[nsub])
        temp=matrix(result_temp$beta,p*(q+1),1)
        id=temp!=0
        if (!is.null(weight)){
          res_temp=sqrt(weight)*(y-result_temp$intercept-E%*%result_temp$alpha-matrix(WW[,id],n,sum(id))%*%temp[id])
          #res_temp=(y-result_temp$intercept-E%*%result_temp$alpha-matrix(WW[,id],n,sum(id))%*%temp[id])
        } else {
          res_temp=y-result_temp$intercept-E%*%result_temp$alpha-matrix(WW[,id],n,sum(id))%*%temp[id]
        }
        objective1=result_temp$objective_set
        # }

        if (is.null(h)){
          if (is.null(weight)){
            nsub=abs(res_temp-stats::median(res_temp))<mu*stats::mad(res_temp)
          } else {
            temp=abs(res_temp[nonzero_id]-stats::median(res_temp[nonzero_id]))<mu*stats::mad(res_temp[nonzero_id])
            temp2=nonzero_id[temp]
            nsub=matrix(0,n,1)
            nsub[temp2]=1
            nsub=nsub==1
          }
        } else {
          temp=order(abs(res_temp))[1:h]
          nsub=matrix(0,n,1)
          nsub[temp]=1
          nsub=nsub==1
        }

        diff_c=abs(objective1-objective)
        objective=objective1
        tt=tt+1
      }
      obj_init[ii]=objective1
      subsamp[ii,]=nsub1
      #print(paste0('The ',ii,'th init, sample number=',sum(nsub),' objective=',obj_init[ii]))
    }
  }


  final_samp=which(base::colSums(subsamp)>=max_init*tau)

  result_temp=MCP_Hier(G[final_samp,],E[final_samp,],y[final_samp],W[final_samp,,],WW[final_samp,],lambda1,lambda2,gamma1,gamma2,weight[final_samp])
  beta_estimate=result_temp$beta
  alpha_estimate=result_temp$alpha
  b_estimate=matrix(beta_estimate,p*(q+1),1)

  df=sum(abs(b_estimate)>0)
  BIC=log(result_temp$RSS)+log(n)*sum(abs(b_estimate)>0)/n
  beta=result_temp$beta
  alpha=result_temp$alpha
  intercept=result_temp$intercept
  result=list(intercept=intercept,alpha=alpha,beta=beta,df=df,BIC=BIC,select_sample=final_samp,family=family)
  return(result)
}