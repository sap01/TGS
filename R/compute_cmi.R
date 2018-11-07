#' Compute Conditional Mutual Information (CMI) the way it is done in the implementation of the PCA-CMI algo
#'
#' Compute Conditional Mutual Info between random variables 'v1' and 'v2'
#' given random vector 'vcs'.
#' Therefore, when 'vcs' is NULL, this
#' function returns the mutual info between
#' random variables 'v1' and 'v2'.
#' The implementation of this R function
#' is adapted from the mutual info estimator
#' function used in the PCA-CMI algo [1].
#' The original funciton, namely 'cmi()' [2],
#' is written in MATLAB.
#'
#' @import stats
#'
#' @param v1 random variable 1
#' @param v2 random variable 2
#' @param vcs random vector 'vcs'
#'
#' @return Conditional Mutual info between 'v1' and 'v2'
#'
#' @references
#' [1] Zhang X, Zhao X M, He K, et al. Inferring gene
#' regulatory networks from gene expression data by
#' path consistency algorithm based on conditional
#' mutual information[J]. Bioinformatics, 2012, 28(1): 98-104.
#' Companion website:
#' https://sites.google.com/site/xiujunzhangcsb/software/pca-cmi
#'
#' [2] Function declaration: 'function cmiv=cmi(v1,v2,vcs)',
#' Source code file: http://www.comp-sysbio.org/grn/pca_cmi.m .
#'
#' @examples
#' ComputeCmiPcaCmi(c(3,5),c(4,2))
#'
#' @export
ComputeCmiPcaCmi <- function(v1,v2,vcs) {
  cmiv <- 0

  if(base::nargs() == 2) {
    c1 <-  stats::var(v1)
    c2 <-  stats::var(v2)
    v <-base::cbind(v1,v2)
    c3 <- base::det(stats::cov(v))

    if (c3 != 0) {
      cmiv <- 0.5*(base::log(c1*c2/c3))
    }
    # else {
    #   ## Perfectly correlated
    #   cmiv <- .Machine$double.xmax
    # }

  } else if(base::nargs() == 3) {
    v <-base::cbind(v1,vcs)
    c1 <- base::det(stats::cov(v))
    v <-base::cbind(v2,vcs)
    c2 <- base::det(stats::cov(v))

    if (base::ncol(vcs) > 1) {
      c3 <- base::det(stats::cov(vcs))
    } else {
      c3  <- stats::var(vcs)
    }

    v <- base::cbind(v1,v2,vcs)
    c4 <- base::det(stats::cov(v))

    if ((c3*c4) != 0) {
      cmiv <- 0.5 * base::log(base::abs((c1 * c2) / (c3 * c4)))
    }
  }
  return (cmiv)
}
