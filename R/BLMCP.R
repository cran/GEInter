#' Accommodating missingness in environmental measurements in gene-environment interaction
#' analysis: penalized estimation and selection
#'
#' The joint gene-environment (G-E) interaction analysis approach developed in Liu et al, 2013.
#' To accommodate "main effects, interactions" hierarchy, two types of penalty, group minimax
#' concave penalty (MCP) and MCP are adopted. Specifically, for each G factor, its main effect
#' and corresponding G-E interactions are regarded as a group, where the group MCP is imposed to
#' identify whether this G factor has any effect at all. In addition, the MCP is imposed on the
#' interaction terms to further identify important interactions.
#'
#' @param G Input matrix of \code{p} G measurements consisting of \code{n} rows. Each row is an
#' observation vector.
#' @param E Input matrix of \code{q} environmental risk factors. Each row is an observation
#' vector.
#' @param Y Response variable. A quantitative vector for \code{family="continuous"}. For
#' \code{family="survival"}, \code{Y} should be a two-column matrix with the first column being
#' the log(survival time) and the second column being the censoring indicator. The indicator is a
#' binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @param weight Observation weights.
#' @param lambda1 A user supplied lambda for group MCP, where each main G effect and its
#' corresponding interactions are regarded as a group.
#' @param lambda2 A user supplied lambda for MCP accommodating interaction selection.
#' @param gamma1 The regularization parameter of the group MCP penalty.
#' @param gamma2 The regularization parameter of the MCP penalty.
#' @param max_iter Maximum number of iterations.
#' @return An object with S3 class \code{"BLMCP"} is returned, which is a list with the following components.
#' \item{call}{The call that produced this object.}
#' \item{alpha}{Matrix of the coefficients for main E effects.}
#' \item{beta}{The matrix of the regression coefficients for all main G effects (the first row)
#' and interactions.}
#' \item{df}{The number of nonzeros.}
#' \item{BIC}{Bayesian Information Criterion.}
#' \item{aa}{The indicator representing whether the algorithm reaches convergence.}
#' @seealso \code{predict}, and \code{coef}, and \code{plot}, and \code{bic.BLMCP} and
#' \code{Augmentated.data} methods.
#' @references Mengyun Wu, Yangguang Zang, Sanguo Zhang, Jian Huang, and Shuangge Ma.
#' \emph{Accommodating missingness in environmental measurements in gene-environment interaction
#' analysis. Genetic Epidemiology, 41(6):523-554, 2017.}\cr  Jin Liu, Jian Huang, Yawei Zhang,
#' Qing Lan, Nathaniel Rothman, Tongzhang Zheng, and Shuangge Ma.
#' \emph{Identification of gene-environment interactions in cancer studies using penalization.
#' Genomics, 102(4):189-194, 2013.}
#' @export
#' @export BLMCP
#'
#' @examples
#' set.seed(100)
#' sigmaG=AR(0.3,100)
#' G=MASS::mvrnorm(250,rep(0,100),sigmaG)
#' E=matrix(rnorm(250*5),250,5)
#' E[,2]=E[,2]>0;E[,3]=E[,3]>0
#' alpha=runif(5,2,3)
#' beta=matrix(0,5+1,100);beta[1,1:8]=runif(8,2,3)
#' beta[2:4,1]=runif(3,2,3);beta[2:3,2]=runif(2,2,3);beta[5,3]=runif(1,2,3)
#'
#' # continuous with Normal error
#' y1=simulated_data(G,E,alpha,beta,error=rnorm(250),family="continuous")
#' fit1<-BLMCP(G,E,y1,weight=NULL,lambda1=0.05,lambda2=0.06,gamma1=3,gamma2=3,max_iter=200)
#' coef1=coef(fit1)
#' y1_hat=predict(fit1,E,G)
#' plot(fit1)
#'
#' # survival with Normal error
#' y2=simulated_data(G,E,alpha,beta,rnorm(250,0,1),family="survival",0.7,0.9)
#' fit2<-BLMCP(G,E,y2,weight=NULL,lambda1=0.05,lambda2=0.06,gamma1=3,gamma2=3,max_iter=200)
#' coef2=coef(fit2)
#' y2_hat=predict(fit2,E,G)
#' plot(fit2)
BLMCP<-function(G,E,Y,weight=NULL,lambda1,lambda2,gamma1=6,gamma2=6,max_iter=200){
  # get call and family
  thisCall = match.call()
  y=as.matrix(Y)
  #response = match.arg(response)
  n=dim(G)[1]
  q=dim(E)[2]
  p=dim(G)[2]
  if (dim(y)[2]==2){
    delta=y[,2]
    y=y[,1]

  } else {
    delta=matrix(1,n,1)
  }
  if(is.null(weight)) {
    weight=kmw(y,delta)
  }

  X=E
  Z=G
  Y=y

  X=X[weight>0,]
  Z=Z[weight>0,]
  Y=Y[weight>0]
  weight=weight[weight>0]

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

  weigth=matrix(weight,n,1)
  m_X=t(weight)%*%X/sum(weight)
  m_W=t(weight)%*%W/sum(weight)
  m_y=t(weight)%*%y/sum(weight)
  X = (X-matrix(1,n,1)%*%m_X)*(sqrt(weight)%*%matrix(1,1,q))
  W = (W-matrix(1,n,1)%*%m_W)*(sqrt(weight)%*%matrix(1,1,pp))
  y = (y-matrix(1,n,1)%*%m_y)*sqrt(weight)


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


  beta0=matrix(0,pp,1);
  alpha0=solve(t(X)%*%X,t(X)%*%y);
  r=y-X%*%alpha0;
  aa=0
  beta1=beta0

  i=1;
  diff=1;
  while (i<=max_iter && diff>=1e-5){
    for (j in 1:p){
      if (j<p){
        id=group_inf[j]:(group_inf[j+1]-1)
      }
      else {
        id=group_inf[j]:pp
      }
      W_j=W[,id]
      r_j=r+W_j%*%beta0[id]
      b_j=beta0[id]
      u_j=1/n*t(W_j)%*%r_j
      g=1+sqrt(q+1)*lambda1*soft_threshold(1/(sqrt(sum(b_j^2))+1e-100),sqrt(sum(b_j^2))/(gamma1*sqrt(q+1)*lambda1));
      v_j=u_j;
      for (k in 2:(q+1)){
        if (u_j[k]<=gamma2*lambda2*g){
          v_j[k]=soft_threshold(u_j[k],lambda2/abs(u_j[k]))/(1-1/(gamma2*g));
        }
      }
      if (sqrt(sum(v_j^2))>gamma1*sqrt(q+1)*lambda1){
        beta1[id]=v_j;}
      else {
        beta1[id]=gamma1/(gamma1-1)*soft_threshold(v_j,sqrt(q+1)*lambda1/sqrt(sum(v_j^2)))
      }
      r=r-W_j%*%(beta1[id]-beta0[id]);
    }
    temp=r+X%*%alpha0;
    alpha1=solve(t(X)%*%X,t(X)%*%temp);
    r=r-X%*%(alpha1-alpha0);
    diff=mean(abs(beta1)-abs(beta1));


    beta0=beta1;
    alpha0=alpha1;


    if (length(which(abs(beta1)>0))>200){  #the algorithm may be not convergence if the values of lambda1 and lambda2 are not appropriate.
      aa=1;
      break
    }

    i=i+1;
  }


  beta_print=beta1;

  for (j in 1:p){
    if (j<p){
      id=group_inf[j]:(group_inf[j+1]-1)
    }
    else {
      id=group_inf[j]:pp
    }
    beta_print[id]=solve(R[[j]],beta1[id])
  }





  df=sum(abs(beta_print)>0)
  BIC=log(sum(r^2))+log(n)*sum(abs(beta_print)>0)/n
  bb=matrix(beta_print,q+1,p)
  result=list(alpha=alpha1,beta=bb,df=df,BIC=BIC,aa=aa)
  class(result) = "BLMCP"
  return(result)
}
