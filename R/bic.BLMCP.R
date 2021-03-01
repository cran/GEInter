#' BIC for BLMCP
#'
#' Selects a point along the regularization path of a fitted \code{BLMCP} object according to
#' the BIC.
#' @param G Input matrix of \code{p} genetic (G) measurements consisting of \code{n} rows.
#' Each row is an observation vector.
#' @param E Input matrix of \code{q} environmental (E) risk factors. Each row is an observation
#' vector.
#' @param Y Response variable. A quantitative vector for continuous response. For survival response, \code{Y} should be a two-column matrix with the first column being
#' the log(survival time) and the second column being the censoring indicator. The indicator is a
#' binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @param weight Observation weights.
#' @param lambda1_set A user supplied lambda sequence for group minimax concave penalty (MCP),
#' where each main G effect and its corresponding interactions are regarded as a group.
#' @param lambda2_set A user supplied lambda sequence for MCP accommodating interaction selection.
#' @param nlambda1 The number of lambda1 values.
#' @param nlambda2 The number of lambda2 values.
#' @param gamma1 The regularization parameter of the group MCP penalty.
#' @param gamma2 The regularization parameter of the MCP penalty.
#' @param max_iter Maximum number of iterations.
#' @return An object with S3 class \code{"bic.BLMCP"} is returned, which is a list with the ingredients of the BIC fit.
#' \item{call}{The call that produced this object.}
#' \item{alpha}{The matrix of the coefficients for main E effects, each column corresponds to one
#' combination of (lambda1,lambda2).}
#' \item{beta}{The coefficients for main G effects and G-E interactions, each column corresponds to
#' one combination of (lambda1,lambda2). For each column, the first element is the first G effect and
#' the second to (\code{q+1}) elements are the interactions for the first G factor, and so on.}
#' \item{df}{The number of nonzeros for each value of (lambda1,lambda2).}
#' \item{BIC}{Bayesian Information Criterion for each value of (lambda1,lambda2).}
#' \item{alpha_estimate}{Final alpha estimate using Bayesian Information Criterion.}
#' \item{beta_estimate}{Final beta estimate using Bayesian Information Criterion.}
#' \item{lambda_combine}{The matrix of (lambda1, lambda2), with the first column being the values of
#' lambda1, the second being the values of lambda2.}
#' @seealso \code{predict}, \code{coef} and \code{plot} methods,
#' and the \code{BLMCP} function.
#' @references Mengyun Wu, Yangguang Zang, Sanguo Zhang, Jian Huang, and Shuangge Ma.
#' \emph{Accommodating missingness in environmental measurements in gene-environment interaction
#' analysis. Genetic Epidemiology, 41(6):523-554, 2017.}\cr Jin Liu, Jian Huang, Yawei Zhang, Qing
#' Lan, Nathaniel Rothman, Tongzhang Zheng, and Shuangge Ma.
#' \emph{Identification of gene-environment interactions in cancer studies using penalization.
#' Genomics, 102(4):189-194, 2013.}
#' @export bic.BLMCP
#' @examples
#' set.seed(100)
#' sigmaG=AR(0.3,50)
#' G=MASS::mvrnorm(150,rep(0,50),sigmaG)
#' E=matrix(rnorm(150*5),150,5)
#' E[,2]=E[,2]>0;E[,3]=E[,3]>0
#' alpha=runif(5,2,3)
#' beta=matrix(0,5+1,100);beta[1,1:8]=runif(8,2,3)
#' beta[2:4,1]=runif(3,2,3)
#' beta[2:3,2]=runif(2,2,3)
#' beta[5,3]=runif(1,2,3)
#'
#' # continuous with Normal error
#' y1=simulated_data(G=G,E=E,alpha=alpha,beta=beta,error=rnorm(150),family="continuous")
#'
#' # survival with Normal error
#' y2=simulated_data(G,E,alpha,beta,rnorm(150,0,1),family="survival",0.8,1)
#'
#' # continuous
#' fit1<-bic.BLMCP(G,E,y1,weight=NULL,lambda1_set=NULL,lambda2_set=NULL,
#' nlambda1=10,nlambda2=10,gamma1=6,gamma2=6,max_iter=200)
#' coef1=coef(fit1)
#' y1_hat=predict(fit1,E,G)
#' plot(fit1)
#'
#' \donttest{
#' ## survival
#' fit2<-bic.BLMCP(G,E,y2,weight=NULL,lambda1_set=NULL,lambda2_set=NULL,
#' nlambda1=20,nlambda2=20,gamma1=6,gamma2=6,max_iter=200)
#' coef2=coef(fit2)
#' y2_hat=predict(fit2,E,G)
#' plot(fit2)
#' }
bic.BLMCP<-function(G,E,Y,weight=NULL,lambda1_set=NULL,lambda2_set=NULL,nlambda1=20,nlambda2=20,gamma1=6,gamma2=6,max_iter=200){

  # get call and family
  thisCall = match.call()
  #response = match.arg(response)
  n=dim(G)[1]
  q=dim(E)[2]
  p=dim(G)[2]
  y=Y
  y=as.matrix(y)
  if (dim(y)[2]==2){
    delta=y[,2]
    y=y[,1]

  } else {
    delta=matrix(1,n,1)
  }
  if(is.null(weight)) {
    weight=kmw(y,delta)
  }

  if (is.null(lambda1_set)){
    lambda_set<-lambda_funtion(G,E,y,weight,nlambda1,nlambda2)
    lambda1_set=lambda_set$lambda1_set
    lambda2_set=lambda_set$lambda2_set
  }



  lambda1_n=length(lambda1_set);
  lambda2_n=length(lambda2_set)
  alpha_f=matrix(0,q,(lambda1_n*lambda2_n));
  pp=p*(q+1)
  beta_f=matrix(0,pp,lambda1_n*lambda2_n);
  lambda_combine=matrix(0,lambda1_n*lambda2_n,2);
  df=matrix(0,lambda1_n*lambda2_n,1);
  BIC=1e+100*matrix(1,lambda1_n*lambda2_n,1);
  for (kk1 in 1:length(lambda1_set)){
    lambda1=lambda1_set[kk1]
    for (kk2 in 1:length(lambda2_set)){
      lambda2=lambda2_set[kk2]

      fit<-BLMCP(G,E,y,weight,lambda1,lambda2,gamma1,gamma2,max_iter)

      alpha_f[,(kk1-1)*nlambda2+kk2]=fit$alpha;
      beta_f[,(kk1-1)*nlambda2+kk2]=matrix(fit$beta,pp,1)
      lambda_combine[(kk1-1)*nlambda2+kk2,1]=lambda1;
      lambda_combine[(kk1-1)*nlambda2+kk2,2]=lambda2;
      df[(kk1-1)*nlambda2+kk2]=fit$df
      BIC[(kk1-1)*nlambda2+kk2]=fit$BIC

      if (fit$aa==1){
        break
      }
    }
  }
  id=which.min(BIC)
  beta_estimate=matrix(beta_f[,id],(1+q),p)
  alpha_estimate=matrix(alpha_f[,id],q,1)

  ##change
  if(!(is.null(colnames(G)))){
    colnames(beta_estimate)=colnames(G)
  }else{
    cnames=paste("G",1:p,sep="")
    colnames(beta_estimate)=cnames
  }

  if(!(is.null(colnames(E)))){
    rownames(alpha_estimate)=colnames(E)
    cnames=c("G",colnames(E))
    rownames(beta_estimate)=cnames
  }else{
    cnames=paste("E",1:q,sep="")
    rownames(alpha_estimate)=cnames
    cnames=c("G",cnames)
    rownames(beta_estimate)=cnames
  }
  #########


  colnames(lambda_combine)=c("lambda1","lambda2")
  result=list(call=thisCall,alpha=alpha_f,beta=beta_f,df=df,BIC=BIC,alpha_estimate=alpha_estimate,beta_estimate=beta_estimate,lambda_combine=lambda_combine)

  class(result) = "bic.BLMCP"
  return(result)
}
