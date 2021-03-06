% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Miss.boosting.R
\name{Miss.boosting}
\alias{Miss.boosting}
\title{Robust gene-environment interaction analysis approach via sparse boosting, where the
missingness in environmental measurements is effectively accommodated using multiple
imputation approach}
\usage{
Miss.boosting(
  G,
  E,
  Y,
  im_time = 10,
  loop_time = 500,
  num.knots = c(2),
  Boundary.knots,
  degree = c(2),
  v = 0.1,
  tau,
  family = c("continuous", "survival"),
  knots = NULL,
  E_type
)
}
\arguments{
\item{G}{Input matrix of \code{p} genetic measurements consisting of \code{n} rows. Each row
is an observation vector.}

\item{E}{Input matrix of \code{q} environmental risk factors. Each row is an observation
vector.}

\item{Y}{Response variable. A quantitative vector for \code{family}="continuous". For
\code{family}="survival", \code{Y} should be a two-column matrix with the first column being
the log(survival time) and the second column being the censoring indicator. The indicator is a
binary variable, with "1" indicating dead, and "0" indicating right censored.}

\item{im_time}{Number of imputation for accommodating missingness
in E variables.}

\item{loop_time}{Number of iterations of the sparse boosting.}

\item{num.knots}{Numbers of knots for the B spline basis.}

\item{Boundary.knots}{The boundary of knots for the B spline basis.}

\item{degree}{Degree for the B spline basis.}

\item{v}{The step size used in the sparse boosting process. Default is 0.1.}

\item{tau}{Threshold used in the stability selection at the third step.}

\item{family}{Response type of \code{Y} (see above).}

\item{knots}{List of knots for the B spline basis. Default is NULL and knots can be generated
with the given \code{num.knots, degree} and \code{Boundary.knots}.}

\item{E_type}{A vector indicating the type of each E factor, with "ED" representing discrete E factor, and "EC" representing continuous E factor.}
}
\value{
An object with S3 class \code{"Miss.boosting"} is returned, which is a list with the following components
\item{call}{The call that produced this object.}
\item{alpha0}{A vector with each element indicating whether the corresponding E factor is
selected.}
\item{beta0}{A vector with each element indicating whether the corresponding G factor or G-E
interaction is selected. The first element is the first G effect and the second to
(\code{q+1}) elements are the interactions for the first G factor, and so on.}
\item{intercept}{The intercept estimate.}
\item{unique_variable}{A matrix with two columns that represents the variables that are
selected for the model after removing the duplicates, since the \code{loop_time} iterations of
the method may produce variables that are repeatedly selected into the model. Here, the first
and second columns correspond to the indexes of E factors and G factors. For example, (1, 0)
represents that this variable is the first E factor, and (1,2) represents that the variable is
the interaction between the first E factor and second G factor.}
\item{unique_coef}{Coefficients corresponding to \code{unique_variable}. Here, the coefficients are simple regression coefficients for the linear effect (discrete E factor, G factor, and
their interaction), and B spline coefficients for the nonlinear effect (continuous E factor,
and corresponding G-E interaction).}
\item{unique_knots}{A list of knots corresponding to \code{unique_variable}. Here, when the
type of \code{unique_variable} is discrete E factor, G factor or their interaction, knot will
be NULL, and knots will be B spline otherwise.}
\item{unique_Boundary.knots}{A list of boundary knots corresponding to
\code{unique_variable}.}
\item{unique_vtype}{A vector representing the variable type of \code{unique_variable}.
Here, "EC" stands for continuous E effect, "ED" for discrete E effect,
"G"  for genetic factor variable, "EC-G" for the interaction between "EC" and "G",
and "ED-G" for the interaction between "ED" and "G".}
\item{degree}{Degree for the B spline basis.}
\item{NorM}{The values of B spline basis.}
\item{E_type}{The type of E effects.}
}
\description{
This gene-environment analysis approach includes three steps to accommodate both missingness
in environmental (E) measurements and long-tailed or contaminated outcomes. At the first step,
the  multiple imputation approach based on sparse boosting method is developed to accommodate
missingness in E measurements, where we use \code{NA} to represent those E measurments which
are missing. Here a semiparametric model is assumed to accommodate nonlinear effects, where we
model continuous E factors in a nonlinear way, and discrete E factors in a linear way. For
estimating the nonlinear functions, the B spline expansion is adopted. At the second step, for
each imputed data, we develop \code{RobSBoosting} approach for identifying important main E
and genetic (G) effects, and G-E interactions, where the Huber loss function and Qn estimator are
adopted to accommodate long-tailed distribution/data contamination (see \code{RobSBoosting}).
At the third step, the identification results from Step 2 are combined based on stability
selection technique.
}
\examples{
data(Rob_data)
G=Rob_data[,1:20];E=Rob_data[,21:24]
Y=Rob_data[,25];Y_s=Rob_data[,26:27]
knots=list();Boundary.knots=matrix(0,(20+4),2)
for (i in 1:4){
  knots[[i]]=c(0,1)
  Boundary.knots[i,]=c(0,1)
}
E2=E1=E

##continuous
E1[7,1]=NA
fit1<-Miss.boosting(G,E1,Y,im_time=1,loop_time=100,num.knots=c(2),Boundary.knots,
degree=c(2),v=0.1,tau=0.3,family="continuous",knots=knots,E_type=c("EC","EC","ED","ED"))
y1_hat=predict(fit1,matrix(E1[1,],nrow=1),matrix(G[1,],nrow=1))
plot(fit1)

\donttest{
##survival
E2[4,1]=NA
fit2<-Miss.boosting(G,E2,Y_s,im_time=2,loop_time=200,num.knots=c(2),Boundary.knots,
degree=c(2),v=0.1,tau=0.3,family="survival",knots,E_type=c("EC","EC","ED","ED"))
y2_hat=predict(fit2,matrix(E1[1,],nrow=1),matrix(G[1,],nrow=1))
plot(fit2)
}
}
\references{
Mengyun Wu and Shuangge Ma.
\emph{Robust semiparametric gene-environment interaction analysis using sparse boosting.
Statistics in Medicine, 38(23):4625-4641, 2019.}
}
