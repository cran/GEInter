#' Make predictions from a "BLMCP" object
#'
#'
#' This function makes predictions from a  BLMCP
#' model, using the stored \code{"BLMCP"} object.
#'
#'
#' @param object Fitted \code{"BLMCP"} object.
#' @param newE Matrix of new values for \code{E} at which predictions are to be
#' made.
#' @param newG Matrix of new values for \code{G} at which predictions are to be
#' made.
#' @param \dots Not used. Other arguments to predict.
#'
#' @return The object returned depends on the \dots{} argument which is passed
#' on to the \code{predict} method for \code{BLMCP} objects.
#' @seealso \code{BLMCP}, \code{coef}, and \code{plot} methods, and \code{bic.BLMCP} method.
#' @references Mengyun Wu, Yangguang Zang, Sanguo Zhang, Jian Huang, and Shuangge Ma.
#' \emph{Accommodating missingness in environmental measurements in gene-environment interaction
#' analysis. Genetic Epidemiology, 41(6):523-554, 2017.}\cr  Jin Liu, Jian Huang, Yawei Zhang,
#' Qing Lan, Nathaniel Rothman, Tongzhang Zheng, and Shuangge Ma.
#' \emph{Identification of gene-environment interactions in cancer studies using penalization. Genomics, 102(4):189-194, 2013.}
#'
#' @method predict BLMCP
#' @export
#' @export predict.BLMCP
#'
predict.BLMCP<-function(object,newE,newG, ...){
  if((missing(newE))|(missing(newG))){
    stop("You need to supply a value for 'newE' and 'newG'")
  }
  newE=as.matrix(newE)
  newG=as.matrix(newG)
  q=dim(newE)[2]
  p=dim(newG)[2]
  n=dim(newE)[1]
  alpha=coef.BLMCP(object)$alpha
  beta=coef.BLMCP(object)$beta
  bb=matrix(beta,q+1,p)
  aa=matrix(alpha,q,1)

  W=matrix(0,n,p*(q+1))
  W[,seq(from=1,to=p*(q+1),by=(q+1))]=newG
  for (i in 1:n){  temp3=matrix(newE[i,],q,1)%*%newG[i,]
  W[i,setdiff(seq(from=1,to=p*(q+1),by=1),seq(from=1,to=p*(q+1),by=q+1))]=matrix(temp3,p*q,1)}

  b=matrix(bb,(q+1)*p,1)
  y.predict=newE%*%matrix(aa,q,1)+W%*%b
  result=y.predict

  return(result)
}
