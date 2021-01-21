#' Extract coefficients from a "bic.BLMCP" object
#'
#' This function extracts the coefficients of main effects and interactions
#' from a BIC BLMCP model, using the stored \code{"bic.BLMCP"} object.
#' @method coef bic.BLMCP
#' @param object Fitted \code{"bic.BLMCP"} model object.
#' @param \dots Not used. Other arguments to get coefficients.
#'
#' @return The object returned depends on the \dots{} argument which is passed on to the \code{coef}
#' method for \code{bic.BLMCP} objects.
#' \item{alpha}{Matrix of the coefficients for main environmental effects.}
#' \item{beta}{The matrix of the regression coefficients for all main genetic effects (the first row)
#' and interactions.}
#' @seealso \code{bic.BLMCP}, and \code{predict}, and \code{plot} methods, and the
#' \code{BLMCP} function.
#' @references Mengyun Wu, Yangguang Zang, Sanguo Zhang, Jian Huang, and Shuangge Ma.
#' \emph{Accommodating missingness in environmental measurements in gene-environment interaction
#' analysis. Genetic Epidemiology, 41(6):523-554, 2017.}\cr  Jin Liu, Jian Huang, Yawei Zhang, Qing
#' Lan, Nathaniel Rothman, Tongzhang Zheng, and Shuangge Ma.
#' \emph{Identification of gene-environment interactions in cancer studies using penalization.
#' Genomics, 102(4):189-194, 2013.}
#' @export
#' @export  coef.bic.BLMCP
coef.bic.BLMCP<-function(object, ...){
  q=length(object$alpha_estimate)
  p=length(object$beta_estimate)/(q+1)
  alpha=matrix(object$alpha_estimate,q,1)
  beta=(matrix(object$beta_estimate,q+1,p))
  cnames1=paste("E",1:q,sep="")
  rownames(alpha)=cnames1
  cnames=paste("G",1:p,sep="")
  rownames(beta)=c("-",cnames1)
  colnames(beta)=cnames
 return(list(alpha=alpha,beta=beta))

}
