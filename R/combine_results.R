combine_results<-function(para_value,tau,p,q){



  im_time=length(para_value)

  for (i in 1:im_time){
    if (i==1){
      unique_variable=para_value[[i]]$unique_variable
      unique_knots=para_value[[i]]$unique_knots
      unique_Boundary.knots=para_value[[i]]$unique_Boundary.knots
      intercept=para_value[[i]]$intercept
      unique_vtype=para_value[[i]]$unique_vtype
    } else {
      unique_variable=rbind(unique_variable,para_value[[i]]$unique_variable)
      unique_knots=c(unique_knots,para_value[[i]]$unique_knots)
      unique_Boundary.knots=c(unique_Boundary.knots,para_value[[i]]$unique_Boundary.knots)
      unique_vtype=c(unique_vtype,para_value[[i]]$unique_vtype)
      intercept=intercept+para_value[[i]]$intercept
    }
  }

  intercept=intercept/im_time

  dp=duplicated(unique_variable,MARGIN=1)
  unique_variable=unique_variable[!dp,]
  unique_knots=unique_knots[!dp]
  unique_Boundary.knots=unique_Boundary.knots[!dp]
  unique_vtype=unique_vtype[!dp]
  unique_coef_set=list()

  if (!is.matrix(unique_variable)){
    unique_variable=matrix(unique_variable,1,2)
  }



  for (i in 1:dim(unique_variable)[1]){
    unique_coef_set[[i]]=0
  }
  unique_coef_id=matrix(0,sum(!dp),im_time)
  for (ii in 1:dim(unique_variable)[1]){
    for (i in 1:im_time){
      temp=which.column(unique_variable[ii,],para_value[[i]]$unique_variable)
      if (length(temp)!=0){
        unique_coef_set[[ii]]=unique_coef_set[[ii]]+para_value[[i]]$unique_coef[[temp]]/im_time
        unique_coef_id[ii,i]=1
      }
    }
  }

  id=which(rowMeans(unique_coef_id)>tau)

  unique_variable=unique_variable[id,]
  unique_knots=unique_knots[id]
  unique_Boundary.knots=unique_Boundary.knots[id]
  unique_vtype=unique_vtype[id]
  unique_coef=unique_coef_set[id]

  if (!is.matrix(unique_variable)){
    unique_variable=matrix(unique_variable,1,2)
  }


  variable_pair=unique_variable
  if (!is.matrix(variable_pair)){
    variable_pair=matrix(variable_pair,1,2)
  }

  alpha0=matrix(0,q,1)
  G_main=matrix(0,p,1)
  GE_interaction=matrix(0,q,p)
  temp=variable_pair[which(variable_pair[,2]==0),1]
  alpha0[temp]=1
  temp=variable_pair[which(variable_pair[,1]==0),2]
  G_main[temp]=1

  temp=variable_pair[((variable_pair[,1]!=0) & (variable_pair[,2]!=0)),]

  GE_interaction[temp]=1

  beta0=matrix(0,p+p*q,1)

  for (j in 1:p) {
    beta0[(j-1)*(q+1)+1]=G_main[j]
    beta0[((j-1)*(q+1)+2):(j*(q+1))]=GE_interaction[,j]
  }

  unique_set=list(intercept=intercept,unique_variable=unique_variable,unique_coef=unique_coef,unique_knots=unique_knots,unique_Boundary.knots=unique_Boundary.knots,unique_vtype=unique_vtype)


  return(list(para_id=c(alpha0,beta0),alpha0=alpha0,beta0=beta0,unique_set=unique_set))

}
