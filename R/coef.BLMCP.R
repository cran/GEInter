#' Extract coefficients from a "BLMCP" object
#'
#' This function extracts the coefficients of main effects and interactions from a BLMCP model,
#' using the stored \code{"BLMCP"} object.
#' @method coef BLMCP
#' @param object Fitted \code{"BLMCP"} model object.
#' @param \dots Not used. Other arguments to get coefficients.
#'
#' @return The object returned depends on the \dots{} argument which is passed on to the \code{coef}
#' method for \code{BLMCP} objects.
#' \item{alpha}{The matrix of the coefficients for main environmental effects.}
#' \item{beta}{The matrix of the regression coefficients for all main genetic effects (the first row)
#' and interactions.}
#' @seealso \code{BLMCP}, and \code{predict}, \code{plot} methods, and
#' \code{bic.BLMCP}.
#' @references Mengyun Wu, Yangguang Zang, Sanguo Zhang, Jian Huang, and Shuangge Ma.
#' \emph{Accommodating missingness in environmental measurements in gene-environment interaction analysis. Genetic Epidemiology, 41(6):523-554, 2017.}\cr Jin Liu, Jian Huang, Yawei Zhang, Qing Lan, Nathaniel Rothman, Tongzhang Zheng, and Shuangge Ma.
#' \emph{Identification of gene-environment interactions in cancer studies using penalization.
#' Genomics, 102(4):189-194, 2013.}
#' @export
#' @export coef.BLMCP

coef.BLMCP<-function(object, ...){
  # q=length(object$alpha)
  # p=length(object$beta)/(q+1)
  bb=object$beta
  aa=object$alpha

  return(list(alpha=aa,beta=bb))

}
