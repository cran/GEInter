% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PTReg.R
\name{PTReg}
\alias{PTReg}
\title{Robust gene-environment interaction analysis using penalized trimmed regression}
\usage{
PTReg(
  G,
  E,
  Y,
  lambda1,
  lambda2,
  gamma1 = 6,
  gamma2 = 6,
  max_init,
  h = NULL,
  tau = 0.4,
  mu = 2.5,
  family = c("continuous", "survival")
)
}
\arguments{
\item{G}{Input matrix of \code{p} genetic (G) measurements consisting of \code{n} rows. Each
row is an observation vector.}

\item{E}{Input matrix of \code{q} environmental (E) risk factors. Each row is an observation
vector.}

\item{Y}{Response variable. A quantitative vector for \code{family="continuous"}. For
\code{family="survival"}, \code{Y} should be a two-column matrix with the first column being
the log(survival time) and the second column being the censoring indicator. The indicator is a
binary variable, with "1" indicating dead, and "0" indicating right censored.}

\item{lambda1}{A user supplied lambda for MCP accommodating main G effect selection.}

\item{lambda2}{A user supplied lambda for MCP accommodating G-E interaction selecton.}

\item{gamma1}{The regularization parameter of the MCP penalty corresponding to G effects.}

\item{gamma2}{The regularization parameter of the MCP penalty corresponding to G-E
interactions.}

\item{max_init}{The number of initializations.}

\item{h}{The number of the trimmed samples if the parameter \code{mu} is not given.}

\item{tau}{The threshold value used in stability selection.}

\item{mu}{The parameter for screening outliers with extreme absolute residuals if the number
of the trimmed samples \code{h} is not given.}

\item{family}{Response type of \code{Y} (see above).}
}
\value{
An object with S3 class \code{"PTReg"} is returned, which is a list with the following components.
\item{call}{The call that produced this object.}
\item{intercept}{The intercept estimate.}
\item{alpha}{The matrix of the coefficients for main E effects.}
\item{beta}{The matrix of the regression coefficients for all main G effects (the first row)
and interactions.}
\item{df}{The number of nonzeros.}
\item{BIC}{Bayesian Information Criterion.}
\item{select_sample}{Selected samples where samples with extreme absolute residuals are
trimmed.}
\item{family}{The same as input \code{family}.}
}
\description{
Gene-environment interaction analysis using penalized trimmed regression, which is robust to
outliers in both predictor and response spaces. The objective function is based on trimming
technique, where the samples with extreme absolute residuals are trimmed. A decomposition
framework is adopted for accommodating "main effects-interactions" hierarchy, and minimax
concave penalty (MCP) is adopted for regularized estimation and interaction (and main genetic
effect) selection.
}
\examples{
sigmaG<-AR(rho=0.3,p=30)
sigmaE<-AR(rho=0.3,p=3)
set.seed(300)
G=MASS::mvrnorm(150,rep(0,30),sigmaG)
EC=MASS::mvrnorm(150,rep(0,2),sigmaE[1:2,1:2])
ED = matrix(rbinom((150),1,0.6),150,1)
E=cbind(EC,ED)
alpha=runif(3,0.8,1.5)
beta=matrix(0,4,30)
beta[1,1:4]=runif(4,1,1.5)
beta[2,c(1,2)]=runif(2,1,1.5)

\donttest{
#continuous response
y1=simulated_data(G=G,E=E,alpha=alpha,beta=beta,error=c(rnorm(130),
rcauchy(20,0,5)),family="continuous")
fit1<-PTReg(G=G,E=E,y1,lambda1=0.3,lambda2=0.3,gamma1=6,gamma2=6,
max_init=50,h=NULL,tau=0.6,mu=2.5,family="continuous")
coef1=coef(fit1)
y_hat1=predict(fit1,E,G)
plot(fit1)

# survival response
y2=simulated_data(G,E,alpha,beta,rnorm(150,0,1),
family="survival",0.7,0.9)
fit2<-PTReg(G=G,E=E,y2,lambda1=0.3,lambda2=0.3,gamma1=6,gamma2=6,
max_init=50,h=NULL,tau=0.6,mu=2.5,family="survival")
coef2=coef(fit2)
y_hat2=predict(fit2,E,G)
plot(fit2)
}

}
\references{
Yaqing Xu, Mengyun Wu, Shuangge Ma, and Syed Ejaz Ahmed.
\emph{Robust gene-environment interaction analysis using penalized trimmed regression. Journal
of Statistical Computation and Simulation, 88(18):3502-3528, 2018.}
}
\seealso{
\code{coef}, \code{predict}, and \code{plot} methods, and \code{bic.PTReg} method.
}
