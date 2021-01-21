boosting_impute<-function(E,im_time,E_type,loop_time,num.knots=NULL,Boundary.knots=NULL,degree=1,knots=NULL){



  model_type="nonlinear"
  Method="LS"
  E_cc=E
  n=dim(E_cc)[1]
  q=dim(E_cc)[2]

  # if(is.null(E_type)){
  # E_type=vector()
  # for (j in 1:q){
  #   if (length(unique(E[,j]))<n/2){
  #     E_type[j]='ED'
  #   }else E_type[j]='EC'
  # }
  # }



  miss_id=which(colSums(is.na(E_cc))!=0)
  if(length(miss_id)>1){
  if(max(miss_id)>length(miss_id)){
      stop("You need to put the missing E factors in the first few columns!")
  }
  }

  p1=1
  p2=q-1



  E_impute=list()

  fit_set=list()

  if (length(miss_id)==1){

    miss=matrix(0,n,1)
    miss[rowSums(is.na(E_cc))==0]=1

    nonmiss_id=which(colSums(is.na(E_cc))==0)
    n_m=sum(miss==0)
    n_c=sum(miss==1)

    E_miss_obs=matrix(E_cc[miss==1,miss_id],n_c,p1) # observed samples for the miss variable  y_train for boosting spline
    E_miss_miss=matrix(E_cc[miss==0,miss_id],n_m,p1) # missing samples for the miss variable  y_test for boosting spline

    E_obs_obs=matrix(E_cc[miss==1,nonmiss_id],n_c,p2) # observed samples for the nonmiss variable x_train for boosting spline
    E_obs_miss=matrix(E_cc[miss==0,nonmiss_id],n_m,p2) # missing samples for the nonmiss variable  x_test for boosting spline

    fit<-boosting.robust.spline(E_obs_obs,E_miss_obs,loop_time,num.knots,knots,Boundary.knots,degree,Method,E_type[nonmiss_id],model_type)

    fit_set[[1]]=fit

    im_mean=predict_boosting2(fit,E_obs_miss,E_type[nonmiss_id])


    nnn=n_c-ceiling(n_c/2)

    err=matrix(0,nnn,1)

    repeat_time=10


    for (rr in 1:repeat_time){
      D1=sample(n_c)[1:ceiling(n_c/2)]
      D2=setdiff(1:n_c,D1)
      fit1<-boosting.robust.spline(E_obs_obs[D1,],E_miss_obs[D1,],loop_time,num.knots,knots,Boundary.knots,degree,Method,E_type[nonmiss_id],model_type)
      predict_cv=predict_boosting2(fit1,E_obs_obs[D2,],E_type[nonmiss_id])
      err[((rr-1)*nnn+1):(rr*nnn)]=E_miss_obs[D2,]-predict_cv
    }

    im_value=matrix(0,n_m,im_time)



    for (ii in 1:im_time){
      idd=1:n_m
      temp=matrix(0,1,n_m)
      nniter=0
      while ((length(idd)>0)&(nniter<200)){
        nniter=nniter+1
        noise_random=emprand(err,length(idd))
        temp[idd]=im_mean[idd]+noise_random
        idd=which((temp>1) | (temp<0))
      }
      im_value[,ii]=temp
    }

    for (i in 1:im_time){
      temp=E_cc
      temp[miss==0,miss_id]=im_value[,i]
      E_impute[[i]]=temp
    }
  } else {


    for (i in 1:im_time){
      E_impute[[i]]=E_cc
    }

    miss_id_s=miss_id

    miss_each=matrix(0,n,q)
    miss_each[!is.na(E_cc)]=1

    miss_rate_each=colSums(miss_each[,miss_id])
    miss_rank=sort(miss_rate_each,index.return=T,decreasing=T)$ix



    for (mm in 1:im_time){

      for (ii in 1:length(miss_id_s)) {

        ii=miss_rank[ii]

        E_temp=E_impute[[mm]]

        miss=matrix(0,n,1)
        miss[rowSums(is.na(E_temp))==0]=1

        nonmiss_id=which(colSums(is.na(E_temp))==0)
        nn_c=sum(miss==1)


        nn_m=sum(miss_each[,ii]==0)

        p2=length(nonmiss_id)

        E_obs_obs=matrix(E_temp[miss==1,nonmiss_id],nn_c,p2) # observed samples for the nonmiss variable x_train for boosting spline

        E_obs_miss=matrix(E_temp[miss_each[,ii]==0,nonmiss_id],nn_m,p2) # missing samples for the nonmiss variable  x_test for boosting spline
        E_miss_miss=matrix(E_temp[miss_each[,ii]==0,ii],nn_m,1) # missing samples for the miss variable  y_test for boosting spline
        E_miss_obs=matrix(E_temp[miss==1,ii],nn_c,p1) # observed samples for the miss variable  y_train for boosting spline


        fit<-boosting.robust.spline(E_obs_obs,E_miss_obs,loop_time,num.knots,knots,Boundary.knots,degree,Method,E_type[nonmiss_id],model_type)

        fit_set[[ii]]=fit
        im_mean=predict_boosting2(fit,E_obs_miss,E_type[nonmiss_id])


        nnn=nn_c-ceiling(nn_c/2)

        err=matrix(0,nnn,1)

        repeat_time=10


        for (rr in 1:repeat_time){
          D1=sample(nn_c)[1:ceiling(nn_c/2)]
          D2=setdiff(1:nn_c,D1)
          fit1<-boosting.robust.spline(E_obs_obs[D1,],E_miss_obs[D1,],loop_time,num.knots,knots,Boundary.knots,degree,Method,E_type[nonmiss_id],model_type)
          predict_cv=predict_boosting2(fit1,E_obs_obs[D2,],E_type[nonmiss_id])
          err[((rr-1)*nnn+1):(rr*nnn)]=E_miss_obs[D2,]-predict_cv
        }


        idd=1:nn_m
        temp=matrix(0,1,nn_m)
        while (length(idd)>0){
          noise_random=emprand(err,length(idd))
          temp[idd]=im_mean[idd]+noise_random
          idd=which((temp>1) | (temp<0))
        }
        im_value=temp



        if (E_type[ii]=='ED'){
          # im_value[im_value>=0.5]=1
          # im_value[im_value<0.5]=0
          im_value=round(im_value)
        }

        temp=E_impute[[mm]]
        temp[miss_each[,ii]==0,ii]=im_value
        E_impute[[mm]]=temp
      }
    }
  }
  return(list(E_impute=E_impute))
}
