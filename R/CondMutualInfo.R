#' Compute Conditional Mutual Infortion (CMI)
#'
#' Compute Conditional Mutual Info between random variables 'v1' and 'v2'
#' given random vector 'vcs'.
#'
#' @import stats
#'
#' @param v1 random variable 1
#' @param v2 random variable 2
#' @param vcs random vector 'vcs'
#'
#' @return Conditional Mutual info between 'v1' and 'v2'
#' @examples
#' computeCmi(c(3,5),c(3,5))
#'
#' @export
computeCmi <- function(v1,v2,vcs)
{
  if(nargs() == 2)
  {
    c1 <-  stats::var(v1)
    c2 <-  stats::var(v2)
    v <-base::cbind(v1,v2)
    c3 <- base::det(stats::cov(v))
    cmiv <- 0.5*(base::log(c1*c2/c3))
  }
  else if(nargs() == 3)
  {
    v <-base::cbind(v1,vcs)
    c1=base::det(stats::cov(v))
    v <-base::cbind(v2,vcs)
    c2=base::det(stats::cov(v))
    if(ncol(vcs >1))
       {
         c3=base::det(stats::cov(vcs))
    }
    else
    {
      c3 =stats::var(vcs)
    }
    v <-base::cbind(v1,v2,vcs)
    c4=base::det(stats::cov(v))
    cmiv=0.5*(base::log((c1*c2)/(c3*c4)))
  }
  return (cmiv)
}
