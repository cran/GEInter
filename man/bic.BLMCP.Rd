% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bic.BLMCP.R
\name{bic.BLMCP}
\alias{bic.BLMCP}
\title{BIC for BLMCP}
\usage{
bic.BLMCP(
  G,
  E,
  Y,
  weight = NULL,
  lambda1_set = NULL,
  lambda2_set = NULL,
  nlambda1 = 20,
  nlambda2 = 20,
  gamma1 = 6,
  gamma2 = 6,
  max_iter = 200
)
}
\arguments{
\item{G}{Input matrix of \code{p} genetic (G) measurements consisting of \code{n} rows.
Each row is an observation vector.}

\item{E}{Input matrix of \code{q} environmental (E) risk factors. Each row is an observation
vector.}

\item{Y}{Response variable. A quantitative vector for continuous response. For survival response, \code{Y} should be a two-column matrix with the first column being
the log(survival time) and the second column being the censoring indicator. The indicator is a
binary variable, with "1" indicating dead, and "0" indicating right censored.}

\item{weight}{Observation weights.}

\item{lambda1_set}{A user supplied lambda sequence for group minimax concave penalty (MCP),
where each main G effect and its corresponding interactions are regarded as a group.}

\item{lambda2_set}{A user supplied lambda sequence for MCP accommodating interaction selection.}

\item{nlambda1}{The number of lambda1 values.}

\item{nlambda2}{The number of lambda2 values.}

\item{gamma1}{The regularization parameter of the group MCP penalty.}

\item{gamma2}{The regularization parameter of the MCP penalty.}

\item{max_iter}{Maximum number of iterations.}
}
\value{
An object with S3 class \code{"bic.BLMCP"} is returned, which is a list with the ingredients of the BIC fit.
\item{call}{The call that produced this object.}
\item{alpha}{The matrix of the coefficients for main E effects, each column corresponds to one
combination of (lambda1,lambda2).}
\item{beta}{The coefficients for main G effects and G-E interactions, each column corresponds to
one combination of (lambda1,lambda2). For each column, the first element is the first G effect and
the second to (\code{q+1}) elements are the interactions for the first G factor, and so on.}
\item{df}{The number of nonzeros for each value of (lambda1,lambda2).}
\item{BIC}{Bayesian Information Criterion for each value of (lambda1,lambda2).}
\item{alpha_estimate}{Final alpha estimate using Bayesian Information Criterion.}
\item{beta_estimate}{Final beta estimate using Bayesian Information Criterion.}
\item{lambda_combine}{The matrix of (lambda1, lambda2), with the first column being the values of
lambda1, the second being the values of lambda2.}
}
\description{
Selects a point along the regularization path of a fitted \code{BLMCP} object according to
the BIC.
}
\examples{
set.seed(100)
sigmaG=AR(0.3,50)
G=MASS::mvrnorm(150,rep(0,50),sigmaG)
E=matrix(rnorm(150*5),150,5)
E[,2]=E[,2]>0;E[,3]=E[,3]>0
alpha=runif(5,2,3)
beta=matrix(0,5+1,50);beta[1,1:8]=runif(8,2,3)
beta[2:4,1]=runif(3,2,3)
beta[2:3,2]=runif(2,2,3)
beta[5,3]=runif(1,2,3)

# continuous with Normal error
y1=simulated_data(G=G,E=E,alpha=alpha,beta=beta,error=rnorm(150),family="continuous")

# survival with Normal error
y2=simulated_data(G,E,alpha,beta,rnorm(150,0,1),family="survival",0.8,1)

# continuous
fit1<-bic.BLMCP(G,E,y1,weight=NULL,lambda1_set=NULL,lambda2_set=NULL,
nlambda1=10,nlambda2=10,gamma1=6,gamma2=6,max_iter=200)
coef1=coef(fit1)
y1_hat=predict(fit1,E,G)
plot(fit1)

\donttest{
## survival
fit2<-bic.BLMCP(G,E,y2,weight=NULL,lambda1_set=NULL,lambda2_set=NULL,
nlambda1=20,nlambda2=20,gamma1=6,gamma2=6,max_iter=200)
coef2=coef(fit2)
y2_hat=predict(fit2,E,G)
plot(fit2)
}
}
\references{
Mengyun Wu, Yangguang Zang, Sanguo Zhang, Jian Huang, and Shuangge Ma.
\emph{Accommodating missingness in environmental measurements in gene-environment interaction
analysis. Genetic Epidemiology, 41(6):523-554, 2017.}\cr Jin Liu, Jian Huang, Yawei Zhang, Qing
Lan, Nathaniel Rothman, Tongzhang Zheng, and Shuangge Ma.
\emph{Identification of gene-environment interactions in cancer studies using penalization.
Genomics, 102(4):189-194, 2013.}
}
\seealso{
\code{predict}, \code{coef} and \code{plot} methods,
and the \code{BLMCP} function.
}
