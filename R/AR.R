#' The covariance matrix with autoregressive (AR) structure among variables
#'
#' The covariance matrix with AR structure among variables, where the marginal variances are 1 and the \code{j}th and \code{k}th variables have correlation coefficient \code{rho^abs(j-k)}.
#'
#' @param rho The correlation coefficient indicating the AR relationship between the variables.
#' @param p Dimension of variables.
#'
#' @return A covariance matrix.
#' @export AR
AR<-function(rho,p){
  sigma<-matrix(0,nrow=p,ncol=p)
for(i in 1:p){
  for(j in 1:i){
    sigma[i,j]<-sigma[j,i]<-rho^(i-j)
  }
}
return(sigma)
}