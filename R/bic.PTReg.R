#' BIC for PTReg
#'
#' Selects a point along the regularization path of a fitted \code{PTReg} object according to
#' the BIC.
#'
#' @param G Input matrix of \code{p} genetic (G) measurements consisting of \code{n} rows.
#' Each row is an observation vector.
#' @param E Input matrix of \code{q} environmental (E) risk factors. Each row is an
#' observation vector.
#' @param Y Response variable. A quantitative vector for \code{family="continuous"}.
#' For \code{family="survival"}, \code{Y} should be a two-column matrix with the first column being
#' the log(survival time) and the second column being the censoring indicator. The indicator is a
#' binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @param lambda1_set A user supplied lambda sequence for minimax concave penalty (MCP)
#' accommodating main G effect selection.
#' @param lambda2_set A user supplied lambda sequence for MCP accommodating interaction selection.
#' @param gamma1 The regularization parameter of the MCP penalty corresponding to G effects.
#' @param gamma2 The regularization parameter of the MCP penalty corresponding to G-E interactions.
#' @param max_init The number of initializations.
#' @param h The number of the trimmed samples if the parameter \code{mu} is not given.
#' @param tau  The threshold value used in stability selection.
#' @param mu The parameter for screening outliers with extreme absolute residuals if the number of
#' the trimmed samples \code{h} is not given.
#' @param family Response type of \code{Y} (see above).
#' @return An object with S3 class \code{"bic.PTReg"} is returned, which is a list with the ingredients of the BIC fit.
#' \item{call}{The call that produced this object.}
#' \item{alpha}{The matrix of the coefficients for main E effects, each column corresponds to one
#' combination of (lambda1,lambda2).}
#' \item{beta}{The coefficients for main G effects and G-E interactions, each column corresponds to
#' one combination of (lambda1,lambda2). For each column, the first element is the first G effect and
#' the second to (\code{q+1}) elements are the interactions for the first G factor, and so on.}
#' \item{intercept}{Matrix of the intercept estimate, each column corresponds to one combination of
#' (lambda1,lambda2).}
#' \item{df}{The number of nonzeros for each value of (lambda1,lambda2).}
#' \item{BIC}{Bayesian Information Criterion for each value of (lambda1,lambda2).}
#' \item{family}{The same as input \code{family}.}
#' \item{intercept_estimate}{Final intercept estimate using Bayesian Information Criterion.}
#' \item{alpha_estimate}{Final alpha estimate using Bayesian Information Criterion.}
#' \item{beta_estimate}{Final beta estimate using Bayesian Information Criterion.}
#' \item{lambda_combine}{Matrix of (lambda1, lambda2), with the first column being the values of
#' lambda1, the second being the values of lambda2.}
#' @references Yaqing Xu, Mengyun Wu, Shuangge Ma, and Syed Ejaz Ahmed.
#' \emph{Robust gene-environment interaction analysis using penalized trimmed regression. Journal of
#' Statistical Computation and Simulation, 88(18):3502-3528, 2018.}
#'
#' @export bic.PTReg
#' @examples
#' sigmaG<-AR(rho=0.3,p=30)
#' sigmaE<-AR(rho=0.3,p=3)
#' set.seed(300)
#' G=MASS::mvrnorm(150,rep(0,30),sigmaG)
#' EC=MASS::mvrnorm(150,rep(0,2),sigmaE[1:2,1:2])
#' ED = matrix(rbinom((150),1,0.6),150,1)
#' E=cbind(EC,ED)
#' alpha=runif(3,0.8,1.5)
#' beta=matrix(0,4,30)
#' beta[1,1:4]=runif(4,1,1.5)
#' beta[2,c(1,2)]=runif(2,1,1.5)
#' lambda1_set=lambda2_set=c(0.2,0.25,0.3,0.35,0.4,0.5)
#'
#' \donttest{
#' #continuous response with outliers/contaminations in response variable
#' y1=simulated_data(G,E,alpha,beta,error=c(rnorm(140),rcauchy(10,0,5)),family="continuous")
#' fit1<-bic.PTReg(G,E,y1,lambda1_set,lambda2_set,gamma1=6,gamma2=6,
#' max_init=50,tau=0.6,mu=2.5,family="continuous")
#' coefficients1=coefficients(fit1)
#' y_predict=predict(fit1,E,G)
#' plot(fit1)
#'
#' # survival with Normal error
#' y2=simulated_data(G,E,alpha,beta,rnorm(150,0,1),family="survival",0.7,0.9)
#' fit2<-bic.PTReg(G,E,y2,lambda1_set,lambda2_set,gamma1=6,gamma2=6,
#' max_init=50,tau=0.6,mu=2.5,family="survival")
#' coefficients2=coefficients(fit2)
#' y_predict=predict(fit2,E,G)
#' plot(fit2)
#' }
bic.PTReg<-function(G,E,Y,lambda1_set,lambda2_set,gamma1,gamma2,max_init,h=NULL,tau=0.4,mu=2.5,family=c("continuous","survival")){
  # get call and family
  thisCall = match.call()
  family = match.arg(family)
  y=as.matrix(Y)

  if (family=="survival"){
    y_s=y[,1]
    delta=y[,2]
    y=y[,1]
    weight=kmw(y,delta)
  } else {
    weight=NULL
  }



  nlam1=length(lambda1_set)

  nlam2=length(lambda2_set)

  n<-dim(E)[1]
  p<-dim(G)[2]
  q<-dim(E)[2]


  W=array(0,dim=c(n,q,p))
  for (i in 1:n){
    temp3=matrix(E[i,],q,1)%*%G[i,]
    W[i,,]=temp3
  }



  WW=matrix(0,n,p*(q+1))
  WW[,seq(from=1,to=p*(q+1),by=(q+1))]=G
  for (i in 1:n){
    temp3=matrix(E[i,],q,1)%*%G[i,]
    WW[i,setdiff(seq(from=1,to=p*(q+1),by=1),seq(from=1,to=p*(q+1),by=q+1))]=matrix(temp3,p*q,1)
  }



  fit_LTS<-ptr1(G,E,Y,lambda1_set[floor(nlam1/2)],lambda2_set[floor(nlam2/2)],gamma1,gamma2,max_init,h,tau,mu=mu,family,W,WW)

  if (!is.null(weight)){
    nn=length(fit_LTS$select_sample)+sum(weight!=0)
  } else {
    nn=length(fit_LTS$select_sample)
  }

  select_sample=fit_LTS$select_sample

  BIC=matrix(0,length(lambda1_set)*length(lambda2_set),1)

  t=1

  fit=list()

  intercept_f=matrix(0,1,length(BIC));
  alpha_f=matrix(0,q,length(BIC));
  pp=p*(q+1)
  beta_f=matrix(0,pp,length(BIC));
  lambda_combine=matrix(0,length(BIC),2);
  colnames(lambda_combine)=c("lambda1","lambda2")
  df=matrix(0,length(BIC),1);



  for (l1 in 1:length(lambda1_set)){
    #print(paste0('l1-',l1))
    for (l2 in 1:length(lambda2_set)){
      #print(paste0('l2-',l2))
      lambda1=lambda1_set[l1]
      lambda2=lambda2_set[l2]


      fit_temp=MCP_Hier(G[select_sample,],E[select_sample,],y[select_sample],W[select_sample,,],WW[select_sample,],lambda1=lambda1,lambda2=lambda2,gamma1=gamma1,gamma2=gamma2,weight=weight[select_sample])


      BIC[t]=log(fit_temp$RSS)+fit_temp$df*log(nn)/nn

      fit[[t]]=fit_temp

      alpha_f[,t]=fit_temp$alpha;
      beta_f[,t]=matrix(fit_temp$beta,pp,1)
      intercept_f[,t]=fit_temp$intercept
      lambda_combine[(l1-1)*length(lambda2_set)+l2,1]=lambda1;
      lambda_combine[(l1-1)*length(lambda2_set)+l2,2]=lambda2;
      df[t]=fit_temp$df

      t=t+1
    }
  }


  id=which.min(BIC)
  beta_estimate=fit[[id]]$beta
  alpha_estimate=matrix(fit[[id]]$alpha,q,1)
  intercept_estimate=fit[[id]]$intercept


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
  result=list(call=thisCall,alpha=alpha_f,beta=beta_f,intercept=intercept_f,df=df,BIC=BIC,family=family,intercept_estimate=intercept_estimate,alpha_estimate=alpha_estimate,beta_estimate=beta_estimate,lambda_combine=lambda_combine)
  class(result) = "bic.PTReg"
  return(result)


}
