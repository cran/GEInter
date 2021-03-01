#' Accommodating missingness in environmental measurements in gene-environment interaction analysis
#'
#' We consider the scenario with missingness in environmental (E) measurements. Our approach
#' consists of two steps. We first develop a nonparametric kernel-based data augmentation
#' approach to accommodate missingness. Then, we adopt a penalization approach \code{BLMCP}
#' for regularized estimation and selection of important interactions and main genetic (G) effects,
#' where the "main effects-interactions" hierarchical structure is respected.
#' As E variables are usually preselected and have a low dimension, selection is not conducted on E
#' variables. With a well-designed weighting scheme, a nice "byproduct" is that the proposed
#' approach enjoys a certain robustness property.
#'
#' @param G Input matrix of \code{p} genetic measurements consisting of \code{n} rows. Each row is an observation vector.
#' @param E Input matrix of \code{q} environmental risk factors. Each row is an observation vector.
#' @param Y Response variable. A quantitative vector for \code{family="continuous"}.
#' For \code{family="survival"}, \code{Y} should be a two-column matrix with the first column
#' being the log(survival time) and the second column being the censoring indicator.
#' The indicator is a binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @param h The bandwidths of the kernel functions with the first and second elements corresponding
#' to the discrete and continuous E factors.
#' @param family Response type of \code{Y} (see above).
#' @param E_type A vector indicating the type of each E factor, with "ED" representing discrete E factor, and "EC" representing continuous E factor.
#' @return
#' \item{E_w}{The augmented data corresponding to \code{E}.}
#' \item{G_w}{The augmented data corresponding to \code{G}.}
#' \item{y_w}{The augmented data corresponding to response \code{y}.}
#' \item{weight}{The weights of the augmented observation data for accommodating missingness and also
#' right censoring if \code{family="survival"}.}
#' @references Mengyun Wu, Yangguang Zang, Sanguo Zhang, Jian Huang, and Shuangge Ma.
#' \emph{Accommodating missingness in environmental measurements in gene-environment interaction
#' analysis. Genetic Epidemiology, 41(6):523-554, 2017.}\cr Jin Liu, Jian Huang, Yawei Zhang, Qing
#' Lan, Nathaniel Rothman, Tongzhang Zheng, and Shuangge Ma.
#' \emph{Identification of gene-environment interactions in cancer studies using penalization.
#' Genomics, 102(4):189-194, 2013.}
#' @export Augmented.data
#'
#' @examples
#' set.seed(100)
#' sigmaG=AR(0.3,50)
#' G=MASS::mvrnorm(100,rep(0,50),sigmaG)
#' E=matrix(rnorm(100*5),100,5)
#' E[,2]=E[,2]>0
#' E[,3]=E[,3]>0
#' alpha=runif(5,2,3)
#' beta=matrix(0,5+1,50)
#' beta[1,1:7]=runif(7,2,3)
#' beta[2:4,1]=runif(3,2,3)
#' beta[2:3,2]=runif(2,2,3)
#' beta[5,3]=runif(1,2,3)
#'
#' # continuous with Normal error N(0,4)
#' y1=simulated_data(G=G,E=E,alpha=alpha,beta=beta,error=rnorm(100,0,4),family="continuous")
#'
#' # survival with Normal error N(0,1)
#' y2=simulated_data(G,E,alpha,beta,rnorm(100,0,1),family="survival",0.7,0.9)
#'
#' # generate E measurements with missingness
#' miss_label1=c(2,6,8,15)
#' miss_label2=c(4,6,8,16)
#' E1=E2=E;E1[miss_label1,1]=NA;E2[miss_label2,1]=NA
#'
#' # continuous
#' data_new1<-Augmented.data(G,E1,y1,h=c(0.5,1), family="continuous",
#' E_type=c("EC","ED","ED","EC","EC"))
#' fit1<-BLMCP(data_new1$G_w, data_new1$E_w, data_new1$y_w, data_new1$weight,
#' lambda1=0.025,lambda2=0.06,gamma1=3,gamma2=3,max_iter=200)
#' coef1=coef(fit1)
#' y1_hat=predict(fit1,E[c(1,2),],G[c(1,2),])
#' plot(fit1)
#'
#' ## survival
#' data_new2<-Augmented.data(G,E2,y2, h=c(0.5,1), family="survival",
#' E_type=c("EC","ED","ED","EC","EC"))
#' fit2<-BLMCP(data_new2$G_w, data_new2$E_w, data_new2$y_w, data_new2$weight,
#' lambda1=0.04,lambda2=0.05,gamma1=3,gamma2=3,max_iter=200)
#' coef2=coef(fit2)
#' y2_hat=predict(fit2,E[c(1,2),],G[c(1,2),])
#' plot(fit2)
Augmented.data<-function(G,E,Y,h,family=c("continuous","survival"),E_type){
  n<-dim(E)[1]
  q<-dim(E)[2]
  p<-dim(G)[2]
  y=Y


  E_id=matrix(1,q,1)
  # for (j in 1:q){
  #   if (length(unique(E[,j]))<n/2){
  #     E_id[j,]=0
  #   }
  # }
  for(j in 1:q){
    if (E_type[j]=="ED"){
        E_id[j,]=0
       }
  }

  if(family=="survival"){
    delta=y[,2]
    y=y[,1]
  } else delta=matrix(1,n,1)
  miss=!is.na(rowSums(E))
  nonmiss_id=which(is.na(colSums(E))==0)
  miss_id=which(is.na(colSums(E))==1)

  miss_each=!is.na(E)

  weight_c<-kmw(y,delta)

  y_estimate<-KMExpection(y,delta)

  xx<-cbind(E[,nonmiss_id],y_estimate)
  E_id=c(E_id[setdiff(1:q,miss_id)],1)



  W=matrix(0,n,p*(q+1))
  W[,seq(from=1,to=p*(q+1),by=(q+1))]=G
  for (i in 1:n){
    temp3=matrix(E[i,],q,1)%*%G[i,]
    W[i,setdiff(seq(from=1,to=p*(q+1),by=1),seq(from=1,to=p*(q+1),by=q+1))]=matrix(temp3,p*q,1)
  }

  X_cc=E
  G_cc=G
  W_cc=W
  y_cc=y
  y_cc=as.matrix(y_cc)
  X_obs=X_cc[miss==1,];# E variables that are completely observed
  G_obs=G_cc[miss==1,];# G variables that are completely observed
  W_obs=W_cc[miss==1,];# G*E variables that are completely observed
  U_obs=cbind(X_obs,W_obs);
  y_obs=y_cc[miss==1]; # responses that are completely observed
  weight_obs_c=weight_c[miss==1]
  delta_obs=delta[miss==1]

  X_miss=X_cc[miss==0,]; #E variables that are partially missing
  y_miss=y_cc[miss==0,]; # responses that are partially missing
  G_miss=G_cc[miss==0,]; # G variables that are completely observed


  miss_each=miss_each[miss==0,]

  xx_obs=xx[miss==1,]; # completely observed variables for completely observed samples
  xx_miss=xx[miss!=1,];# completely observed variables for partially missing samples

  delta_miss=delta[miss==0]
  weight_miss_c=weight_c[miss==0]

  nc=sum(miss==1)

  #######E_id represent ED(0) or EC(1)
  ##nonmiss_id reprents index which is no missing
  temp=rep(0,length(nonmiss_id))
  jj=1
  for(i in nonmiss_id){
    if(E_id[i]==1){
      temp[jj]=h[2]
    }else temp[jj]=h[1]
    jj=jj+1
  }
  aa_c=h[2]
  h=rep(0,(length(nonmiss_id)+1))
  h[1:length(nonmiss_id)]=temp
  h[(length(nonmiss_id)+1)]=aa_c


  w_new=matrix(0,nc,n-nc)
  for (j in 1:nc){
    for (i in 1:(n-nc)){
      temp1=kernel_function(matrix(xx_obs[j,]-xx_miss[i,],1,dim(xx)[2]),E_id,h);
      temp2=kernel_function(xx-matrix(1,n,1)%*%xx_miss[i,],E_id,h);
      w_new[j,i]=temp1/sum(temp2);
    }
  }

  pi=matrix(0,n-nc,1);

  for (i in 1:(n-nc)){
    pi[i]=sum(w_new[,i])
  }

  ww_new=matrix(0,nc*(n-nc),1);
  X_new=matrix(0,nc*(n-nc),q);
  G_new=matrix(0,nc*(n-nc),p);

  y_new=matrix(0,nc*(n-nc),1);

  for (j in 1:nc){
    for (i in 1:(n-nc)){
      temp=matrix(0,1,q)
      temp[miss_each[i,]==0]=X_obs[j,miss_each[i,]==0]
      temp[miss_each[i,]==1]= X_miss[i,miss_each[i,]==1]
      X_new[i+(j-1)*(n-nc),]=temp

      G_new[i+(j-1)*(n-nc),]=G_miss[i,]
      y_new[i+(j-1)*(n-nc)]=y_miss[i]
      ww_new[i+(j-1)*(n-nc)]=w_new[j,i]/sum(pi)*(n-nc)
    }
  }


  X_w=rbind(X_obs, X_new)
  G_w=rbind(G_obs, G_new)
  y_w=rbind(matrix(y_obs,nc,1), y_new)
  ww_new2=rbind(matrix(1,nc,1),ww_new)


  temp3 = matrix(weight_miss_c%*%matrix(1,1,nc),(n-nc)*nc,1);
  kmweight = rbind(matrix(weight_obs_c,nc,1), temp3)
  total_weight=ww_new2*kmweight

  if(!(is.null(colnames(E)))) colnames(X_w)=colnames(E)
  if(!(is.null(colnames(G)))) colnames(G_w)=colnames(G)
  result=list(E_w=X_w,G_w=G_w,y_w=y_w,weight=total_weight)
  return(result)

}
