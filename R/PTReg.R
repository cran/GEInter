#' Robust gene-environment interaction analysis using penalized trimmed regression
#'
#' Gene-environment interaction analysis using penalized trimmed regression, which is robust to
#' outliers in both predictor and response spaces. The objective function is based on trimming
#' technique, where the samples with extreme absolute residuals are trimmed. A decomposition
#' framework is adopted for accommodating "main effects-interactions" hierarchy, and minimax
#' concave penalty (MCP) is adopted for regularized estimation and interaction (and main genetic
#' effect) selection.
#'
#' @param G Input matrix of \code{p} genetic (G) measurements consisting of \code{n} rows. Each
#' row is an observation vector.
#' @param E Input matrix of \code{q} environmental (E) risk factors. Each row is an observation
#' vector.
#' @param Y Response variable. A quantitative vector for \code{family="continuous"}. For
#' \code{family="survival"}, \code{Y} should be a two-column matrix with the first column being
#' the log(survival time) and the second column being the censoring indicator. The indicator is a
#' binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @param lambda1 A user supplied lambda for MCP accommodating main G effect selection.
#' @param lambda2 A user supplied lambda for MCP accommodating G-E interaction selecton.
#' @param gamma1 The regularization parameter of the MCP penalty corresponding to G effects.
#' @param gamma2 The regularization parameter of the MCP penalty corresponding to G-E
#' interactions.
#' @param max_init The number of initializations.
#' @param h The number of the trimmed samples if the parameter \code{mu} is not given.
#' @param tau The threshold value used in stability selection.
#' @param mu The parameter for screening outliers with extreme absolute residuals if the number
#' of the trimmed samples \code{h} is not given.
#' @param family Response type of \code{Y} (see above).

#' @return An object with S3 class \code{"PTReg"} is returned, which is a list with the following components.
#' \item{call}{The call that produced this object.}
#' \item{intercept}{The intercept estimate.}
#' \item{alpha}{The matrix of the coefficients for main E effects.}
#' \item{beta}{The matrix of the regression coefficients for all main G effects (the first row)
#' and interactions.}
#' \item{df}{The number of nonzeros.}
#' \item{BIC}{Bayesian Information Criterion.}
#' \item{select_sample}{Selected samples where samples with extreme absolute residuals are
#' trimmed.}
#' \item{family}{The same as input \code{family}.}
#' @seealso \code{coef}, \code{predict}, and \code{plot} methods, and \code{bic.PTReg} method.
#' @references Yaqing Xu, Mengyun Wu, Shuangge Ma, and Syed Ejaz Ahmed.
#' \emph{Robust gene-environment interaction analysis using penalized trimmed regression. Journal
#' of Statistical Computation and Simulation, 88(18):3502-3528, 2018.}
#' @importFrom stats median
#' @importFrom stats mad
#' @export
#' @export PTReg
#'
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
#'
#' \donttest{
#' #continuous response
#' y1=simulated_data(G=G,E=E,alpha=alpha,beta=beta,error=c(rnorm(130),
#' rcauchy(20,0,5)),family="continuous")
#' fit1<-PTReg(G=G,E=E,y1,lambda1=0.3,lambda2=0.3,gamma1=6,gamma2=6,
#' max_init=50,h=NULL,tau=0.6,mu=2.5,family="continuous")
#' coef1=coef(fit1)
#' y_hat1=predict(fit1,E,G)
#' plot(fit1)
#'
#' # survival response
#' y2=simulated_data(G,E,alpha,beta,rnorm(150,0,1),
#' family="survival",0.7,0.9)
#' fit2<-PTReg(G=G,E=E,y2,lambda1=0.3,lambda2=0.3,gamma1=6,gamma2=6,
#' max_init=50,h=NULL,tau=0.6,mu=2.5,family="survival")
#' coef2=coef(fit2)
#' y_hat2=predict(fit2,E,G)
#' plot(fit2)
#' }
#'



PTReg<-function(G,E,Y,lambda1,lambda2,gamma1=6,gamma2=6,max_init,h=NULL,tau=0.4,mu=2.5,family=c("continuous","survival")){
  # get call and family
  thisCall = match.call()
  family = match.arg(family)
  subsamp=NULL
  W=NULL;WW=NULL;
  y=as.matrix(Y)
  ##survival
  if (family=="survival"){
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




  if (is.null(W)){
    W=array(0,dim=c(n,q,p))
    for (i in 1:n){
      temp3=matrix(E[i,],q,1)%*%G[i,]
      W[i,,]=temp3
    }
  }

  if (is.null(WW)){
    WW=matrix(0,n,(p*(q+1)))
    WW[,seq(from=1,to=p*(q+1),by=(q+1))]=G
    for (i in 1:n){
      temp3=matrix(E[i,],q,1)%*%G[i,]
      WW[i,setdiff(seq(from=1,to=p*(q+1),by=1),seq(from=1,to=p*(q+1),by=q+1))]=matrix(temp3,p*q,1)
    }
  }

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



        result_temp=MCP_Hier(G[nsub,],E[nsub,],y[nsub],W[nsub,,],WW[nsub,],lambda1,lambda2,gamma1,gamma2,weight[nsub])
        temp=matrix(result_temp$beta,p*(q+1),1)
        id=temp!=0
        if (!is.null(weight)){
          res_temp=sqrt(weight)*(y-result_temp$intercept-E%*%result_temp$alpha-matrix(WW[,id],n,sum(id))%*%temp[id])

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


  alpha=matrix(alpha,ncol=1)
  ##change
  if(!(is.null(colnames(G)))){
    colnames(beta)=colnames(G)
  }else{
    cnames=paste("G",1:p,sep="")
    colnames(beta)=cnames
  }

  if(!(is.null(colnames(E)))){
    rownames(alpha)=colnames(E)
    cnames=c("G",colnames(E))
    rownames(beta)=cnames
  }else{
    cnames=paste("E",1:q,sep="")
    rownames(alpha)=cnames
    cnames=c("G",cnames)
    rownames(beta)=cnames
  }


  result=list(call=thisCall,intercept=intercept,alpha=alpha,beta=beta,df=df,BIC=BIC,select_sample=final_samp,family=family)
  class(result) = "PTReg"
  return(result)
}



