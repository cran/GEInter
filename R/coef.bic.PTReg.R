#' Extract coefficients from a "bic.PTReg" object
#'
#' This function extracts the coefficients of main effects and interactions from a BIC PTReg model,
#' using the stored \code{"bic.PTReg"} object.
#'
#' @method coef bic.PTReg
#' @param object Fitted "bic.PTReg" model object.
#' @param \dots Not used. Other arguments to get coefficients.
#'
#' @return The object returned depends on the \dots{} argument which is passed on to the \code{coef}
#' method for \code{bic.PTReg} objects.
#' \item{intercept}{The intercept estimate.}
#' \item{alpha}{Matrix of the coefficients for main environmental effects.}
#' \item{beta}{The matrix of the regression coefficients for all main genetic effects (the first row) and interactions.}
#' @seealso \code{bic.PTReg}, and \code{predict}, and \code{plot} methods, and
#' \code{PTReg}.
#' @references Yaqing Xu, Mengyun Wu, Shuangge Ma, and Syed Ejaz Ahmed.
#' \emph{Robust gene-environment interaction analysis using penalized trimmed regression. Journal of
#' Statistical Computation and Simulation, 88(18):3502-3528, 2018.}
#'
#' @export
#' @export coef.bic.PTReg
coef.bic.PTReg<-function(object,...){
  intercept_estimate=object$intercept_estimate
  alpha_estimate=object$alpha_estimate
  beta_estimate=object$beta_estimate
  return(list(intercept=intercept_estimate,alpha=alpha_estimate,beta=beta_estimate))
}
