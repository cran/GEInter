#' A matrix containing the simulated data for \code{RobSBoosting} and \code{Miss.boosting} methods
#'
#' A matrix containing the simulated genetic (G) effects (the first 20 columns), environmental (E)
#' effects (column 21 to column 24), continuous response (column 25), logarithm of survival time
#' (column 26), and censoring indicator (column 27).
#'
#' @docType data
#'
#' @usage data(Rob_data)
#' @format A matrix with 100 rows and 27 variables.
#' @examples
#' data(Rob_data)
"Rob_data"


#' A data frame containing the TCGA head and neck squamous cell carcinoma (HNSCC) data.
#'
#' A data frame containing the 7 environmental (E)
#' effects (the first 7 columns), 2000 genetic (G) effects (column 8 to column 2007), logarithm of survival time
#' (column 2008), and censoring indicator (column 2009). All of them can be downloaded  from TCGA Provisional using the
#' R package \code{cgdsr}. See details.
#'
#' There are seven E effects, namely alcohol consumption frequency (ACF), smoking pack
#' years (SPY), age, gender, PN, PT, and ICD O3 site. For G effects, 2,000 gene
#' expressions are considered. Among 484 subjects, 343 subjects have missingness in ACF and/or SPY.
#' For G effects, we analyze mRNA gene expressions. A total of 18,409 gene expression
#' measurements are available, then prescreening is conducted using marginal Cox models, finally,
#' the top 2,000 genes with the smallest p-values are selected for downstream analysis.
#'
#' @docType data
#'
#' @usage data(HNSCC)
#' @format A data frame with 484 rows and 2009 variables.
#' @examples
#' data(HNSCC)
#' E=as.matrix(HNSCC[,1:7])
#' G=as.matrix(HNSCC[,8:2007])
#' Y=as.matrix(HNSCC[,2008:2009])
#' \donttest{
#' fit<-Miss.boosting(G,E,Y,im_time=10,loop_time=1000,v=0.25,num.knots=5,degree=3,tau=0.3,
#' family="survival",E_type=c(rep("EC",3),rep("ED",4)))
#' plot(fit)
#' }
"HNSCC"