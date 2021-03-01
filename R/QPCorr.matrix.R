#' Robust identification of gene-environment interactions using a
#' quantile partial correlation approach
#'
#'
#' A robust gene-environment interaction identification approach using the quantile partial
#' correlation technique. This approach is a marginal analysis approach built on the quantile
#' regression technique, which can accommodate long-tailed or contaminated outcomes. For response
#' with right censoring, Kaplan-Meier (KM) estimator-based weights are adopted to easily
#' accommodate censoring. In addition, it adopts partial correlation to identify important
#' interactions while properly controlling for the main genetic (G) and environmental (E) effects.
#' @param G Input matrix of \code{p} G measurements consisting of \code{n} rows. Each row is an
#' observation vector.
#' @param E Input matrix of \code{q} E risk factors. Each row is an observation vector.
#' @param Y Response variable. A quantitative vector for \code{family="continuous"}. For
#' \code{family="survival"}, \code{Y} should be a two-column matrix with the first column being
#' the log(survival time) and the second column being the censoring indicator. The indicator is a
#' binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @param tau Quantile.
#' @param w Weight for accommodating censoring if \code{family}="survival". Default is NULL and a
#' Kaplan-Meier estimator-based weight is used.
#' @param family Response type of \code{Y} (see above).
#'
#' @return Matrix of (censored) quantile partial correlations for interactions.
#' @seealso \code{QPCorr.pval} method.
#' @references Yaqing Xu, Mengyun Wu, Qingzhao Zhang, and Shuangge Ma.
#' \emph{Robust identification of gene-environment interactions for prognosis using a quantile
#' partial correlation approach. Genomics, 111(5):1115-1123, 2019.}
#' @export
#' @export QPCorr.matrix
#'
#' @examples
#' alpha=matrix(0,5,1)
#' alpha[1:2]=1
#' beta=matrix(0,6,100)
#' beta[1,1:5]=1
#' beta[2:3,1:5]=2
#' beta[4:6,6:7]=2
#' sigmaG<-AR(rho=0.3,100)
#' sigmaE<-AR(rho=0.3,5)
#' G<-MASS::mvrnorm(200,rep(0,100),sigmaG)
#' E<-MASS::mvrnorm(200,rep(0,5),sigmaE)
#' e1<-rnorm(200*.05,50,1);e2<-rnorm(200*.05,-50,1);e3<-rnorm(200*.9)
#' e<-c(e1,e2,e3)
#'
#' # continuous
#' y1=simulated_data(G=G,E=E,alpha=alpha,beta=beta,error=e,family="continuous")
#' cpqcorr_stat1<-QPCorr.matrix(G,E,y1,tau=0.5,w=NULL,family="continuous")
#'
#' # survival
#' y2=simulated_data(G,E,alpha,beta,rnorm(200,0,1),family="survival",0.7,0.9)
#' cpqcorr_stat<-QPCorr.matrix(G,E,y2,tau=0.5,w=NULL,family="survival")
#'
QPCorr.matrix<-function(G,E,Y,tau,w=NULL,family=c("continuous","survival")){
  Y=as.matrix(Y)
  n=dim(Y)[1]
  if(family=="survival"){
    delta=Y[,2]
    y=Y[,1]
  } else {delta=matrix(1,n,1);y=Y}

  x=E
  z=G
  q=dim(x)[2]
  p=dim(z)[2]
  all<-matrix(nrow=p,ncol=q)
  ## x_m z_k
  if (is.null(w)){
    w<-reweight(y,delta,tau)
  }
  for(m in 1:q){
    for(k in 1:p){
      all[k,m]<-cqpcorr(y,delta,x[,m],z[,k],tau,w)
    }
  }
  all<-t(all)

  if(!(is.null(colnames(G)))) colnames(all)=colnames(G)
  if(!(is.null(colnames(E)))) rownames(all)=colnames(E)

  return(all)
}
