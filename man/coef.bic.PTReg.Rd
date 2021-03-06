% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.bic.PTReg.R
\name{coef.bic.PTReg}
\alias{coef.bic.PTReg}
\title{Extract coefficients from a "bic.PTReg" object}
\usage{
\method{coef}{bic.PTReg}(object, ...)
}
\arguments{
\item{object}{Fitted "bic.PTReg" model object.}

\item{\dots}{Not used. Other arguments to get coefficients.}
}
\value{
The object returned depends on the \dots{} argument which is passed on to the \code{coef}
method for \code{bic.PTReg} objects.
\item{intercept}{The intercept estimate.}
\item{alpha}{The matrix of the coefficients for main environmental effects.}
\item{beta}{The matrix of the regression coefficients for all main genetic effects (the first row) and interactions.}
}
\description{
This function extracts the coefficients of main effects and interactions from a BIC PTReg model,
using the stored \code{"bic.PTReg"} object.
}
\references{
Yaqing Xu, Mengyun Wu, Shuangge Ma, and Syed Ejaz Ahmed.
\emph{Robust gene-environment interaction analysis using penalized trimmed regression. Journal of
Statistical Computation and Simulation, 88(18):3502-3528, 2018.}
}
\seealso{
\code{bic.PTReg}, and \code{predict}, and \code{plot} methods, and
\code{PTReg}.
}
