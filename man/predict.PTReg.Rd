% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.PTReg.R
\name{predict.PTReg}
\alias{predict.PTReg}
\title{Make predictions from a "PTReg" object}
\usage{
\method{predict}{PTReg}(object, newE, newG, ...)
}
\arguments{
\item{object}{Fitted \code{"PTReg"} object.}

\item{newE}{Matrix of new values for \code{E} at which predictions are to be
made.}

\item{newG}{Matrix of new values for \code{G} at which predictions are to be
made.}

\item{\dots}{Not used. Other arguments to predict.}
}
\value{
The object returned depends on the \dots{} argument which is passed
on to the \code{predict} method for \code{PTReg} objects.
}
\description{
This function makes predictions from a PTReg model, using the stored \code{"PTReg"} object.
}
\references{
Yaqing Xu, Mengyun Wu, Shuangge Ma, and Syed Ejaz Ahmed.
\emph{Robust gene-environment interaction analysis using penalized trimmed regression. Journal
of Statistical Computation and Simulation, 88(18):3502-3528, 2018.}
}
\seealso{
\code{PTReg}, \code{coef} and \code{plot} methods, and \code{bic.PTReg}.
}
