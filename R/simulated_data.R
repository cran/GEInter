#' Simulated data for generating response
#'
#' Generate simulated response.
#' @param G Input matrix of \code{p} genetic (G) measurements consisting of \code{n} rows. Each
#' row is an observation vector.
#' @param E Input matrix of \code{q} environmental (E) risk factors. Each row is an observation
#' vector.
#' @param alpha  Matrix of the true coefficients for main E effects.
#' @param beta Matrix of the true regression coefficients for all main G effects (the first row) and interactions.
#' @param error Error terms.
#' @param family Type of the response variable. If \code{family="continuous"}, a quantitative
#' vector is generated. If \code{family="survival"}, a two-column matrix with the first column
#' being the log(survival time) and the second column being the censoring indicator is
#' generated.The indicator is a binary variable, with "1" indicating dead, and "0" indicating
#' right censored.
#' @param a1 If \code{family}="survival", we generate the censoring time from a uniform
#' distribution where \code{a1} is the left endpoint.
#' @param a2 If \code{family}="survival", we generate the censoring time from a uniform
#' distribution where \code{a2} is the right endpoint.
#' @return Response variable. A quantitative vector for \code{family="continuous"}. For
#' \code{family="survival"}, it would be a two-column matrix with the first column being the
#' log(survival time) and the second column being the censoring indicator. The indicator
#' is a binary variable, with "1" indicating dead, and "0" indicating right censored.
#' @import stats
#' @export simulated_data
simulated_data<-function(G,E,alpha,beta,error,family= c("continuous","survival"),a1=NULL,a2=NULL){
  n=dim(G)[1]
  p=dim(G)[2]
  q=dim(E)[2]
  b=matrix(beta,(q+1)*p,1);
  W=matrix(0,n,p*(q+1));
  W[,seq(from=1,to=p*(q+1),by=(q+1))]=G
  for (i in 1:n){
    temp3=matrix(E[i,],q,1)%*%G[i,]
    W[i,setdiff(seq(from=1,to=p*(q+1),by=1),seq(from=1,to=p*(q+1),by=q+1))]=matrix(temp3,p*q,1)
  }
  WW=E%*%matrix(alpha,q,1)+W%*%b
  TT=WW+error;
  if (family=="continuous"){
    y=TT
    }else {
    C=stats::runif(n,stats::quantile(TT,a1),stats::quantile(TT,a2))
    y=matrix(0,n,2)
    y[,1]=pmin(TT,C);
    y[,2]=TT<=C+0
  }
  return(y=y)
}
