#' Extract coefficients from a "PTReg" object
#'
#' This function extracts main effect and interaction coefficients from a PTReg model, using the
#' stored \code{"PTReg"} object.
#' @method coef PTReg
#' @param object Fitted \code{"PTReg"} model object.
#' @param \dots Not used. Other arguments to get coefficients.
#'
#' @return The object returned depends on the \dots{} argument which is passed on to the \code{coef}
#' method for \code{PTReg} objects.
#' \item{intercept}{The intercept estimate.}
#' \item{alpha}{Matrix of the coefficients for main environmental effects.}
#' \item{beta}{The matrix of the regression coefficients for all main genetic effects (the first row) and interactions.}
#' @seealso \code{PTReg}, and \code{predict} methods, and
#' \code{bic.PTReg}.
#' @references Yaqing Xu, Mengyun Wu, Shuangge Ma, and Syed Ejaz Ahmed.
#' \emph{Robust gene-environment interaction analysis using penalized trimmed regression. Journal of
#' Statistical Computation and Simulation, 88(18):3502-3528, 2018.}
#' @export
#' @export coef.PTReg
coef.PTReg<-function(object,...){
  beta_estimate=object$beta
  p=dim(beta_estimate)[2]
  q=dim(beta_estimate)[1]-1
  alpha_estimate=object$alpha
  intercept_estimate=object$intercept
  #b_estimate=matrix(beta_estimate,p*(q+1),1)
  return(list(intercept=intercept_estimate,alpha=alpha_estimate,beta=beta_estimate))
  }
