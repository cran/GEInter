#' P-values of the "QPCorr.matrix" obtained using a permutation approach
#'
#' P-values of the \code{"QPCorr.matrix "} obtained using a permutation approach, the
#' interactions with smaller p-values are regarded as more important.
#' @param G Input matrix of \code{p} genetic (G) measurements consisting of \code{n} rows. Each row is an observation vector.
#' @param E Input matrix of \code{q} environmental (E) risk factors, each row is an observation
#' vector.
#' @param Y Response variable. A quantitative vector for \code{family="continuous"}. For
#' \code{family="survival"}, \code{Y} should be a two-column matrix with the first column being
#' the log(survival time) and the second column being the censoring indicator. The indicator is a
#' binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @param tau Quantile.
#' @param w Weight for accommodating censoring if \code{family="survival"}. Default is NULL and a
#' Kaplan-Meier estimator-based weight is used.
#' @param permutation_t Permutation times.
#' @param family Response type of \code{Y} (see above).
#'
#' @return Matrix of p-value, with the element in the \code{i}th row and the \code{j} column
#' represents the p-value of the (censored) quantile partial correlation corresponding to the
#' \code{i}th E and the \code{j}th G.
#' @seealso \code{QPCorr.matrix } method.
#' @references Yaqing Xu, Mengyun Wu, Qingzhao Zhang, and Shuangge Ma.
#' \emph{Robust identification of gene-environment interactions for prognosis using a quantile
#' partial correlation approach. Genomics, 111(5):1115-1123, 2019.}
#' @examples
#' n=50
#' alpha=matrix(0,5,1)
#' alpha[1:2]=1
#' beta=matrix(0,6,20)
#' beta[1,1:4]=1
#' beta[2:3,1:4]=2
#' sigmaG<-AR(rho=0.3,20)
#' sigmaE<-AR(rho=0.3,5)
#' G<-MASS::mvrnorm(n,rep(0,20),sigmaG)
#' E<-MASS::mvrnorm(n,rep(0,5),sigmaE)
#' e1<-rnorm(n*.05,50,1);e2<-rnorm(n*.05,-50,1);e3<-rnorm((n-length(e1)-length(e2)))
#' e<-c(e1,e2,e3)
#'
#' \donttest{
#' # continuous
#' y1=simulated_data(G=G,E=E,alpha=alpha,beta=beta,error=e,family="continuous")
#' cpqcorr_pvalue1<-QPCorr.pval(G,E,y1,tau=0.5,permutation_t=500,family="continuous")
#'
#' # survival
#' y2=simulated_data(G,E,alpha,beta,rnorm(n,0,1),family="survival",0.7,0.9)
#' cpqcorr_pvalue2<-QPCorr.pval(G,E,y2,tau=0.5,permutation_t=500,family="survival")
#' }
#' @export
#' @export QPCorr.pval
#'
QPCorr.pval<-function(G,E,Y,tau,w=NULL,permutation_t=1000,family=c("continuous","survival")) {
  family = match.arg(family)
  n=dim(Y)[1]
  if(family=="survival"){
    delta=Y[,2]
    y=Y[,1]
  } else {delta=matrix(1,n,1);y=Y}

  p=dim(G)[2]
  q=dim(E)[2]
  x=E;z=G
  if (is.null(w)){
    weight<-reweight(y,delta,tau)
  }

  all<-t(QPCorr.matrix (z,x,Y,tau,weight,family))

  ##### to reduce time of permutation
  #all_corr<-lapply(1:permutation_t,function(x) matrix(rep(NA,p*q),nrow=p))

  response<-cbind(y,delta)

  pvalf<-matrix(1e+10,p,q)

  p_freq<-matrix(0,p,q)

  for(i in 1:100){
    porder<-sample(1:n,n,replace=F)
    new<-response[porder,]
    if(family=="continuous") new=new[,1]
    temp<-t(QPCorr.matrix (z,x,new,tau,weight[porder],family))
    p_freq<-p_freq+(abs(all)<=abs(temp))

    #all_corr[[i]]<-cpqcorr.matrix(new[,1],new[,2],x,z,tau,weight[porder])
    #print(i)
  }
  ### compute p-val and then delete half

  temp<-p_freq/100

  pvalf[temp>0.5]=temp[temp>0.5]

  keepp1<-matrix(as.numeric(temp<=0.5),nrow=p)

  weights<-reweight(response[,1],response[,2],tau)

  for(i in 1:200){
    porder<-sample(1:n,n,replace=F)
    new<-response[porder,]
    weig<-weights[porder]
    temp=matrix(1e+10,p,q)
    for(k in 1:p){
      for(m in 1:q){
        if(keepp1[k,m]==1){
          temp[k,m]<-cqpcorr(new[,1],new[,2],x[,m],z[,k],tau,weig)
        }
      }
    }
    p_freq<-p_freq+(abs(all)<=abs(temp))
    #print(c("step2",i))
  }


  temp<-p_freq/300
  pvalf[(keepp1==1) & (temp>0.15)]=temp[(keepp1==1) & (temp>0.15)]
  keepp2<-matrix(as.numeric(((keepp1==1) & (temp<=0.15))),nrow=p)


  keepp4<-keepp2

  #print('step3')

  for(i in 1:(permutation_t-300)){
    if(i%%1000==0){
      print(c("step3",i))
    }
    porder<-sample(1:n,n,replace=F)
    new<-response[porder,]
    weig<-weights[porder]
    temp=matrix(1e+10,p,q)
    for(k in 1:p){
      for(m in 1:q){
        if(keepp4[k,m]==1){
          temp[k,m]<-cqpcorr(new[,1],new[,2],x[,m],z[,k],tau,weig)
        }
      }
    }
    p_freq<-p_freq+(abs(all)<=abs(temp))
  }

  temp=p_freq/permutation_t

  pvalf[keepp4==1]=temp[keepp4==1]

  return(pvalf=t(pvalf))
}
