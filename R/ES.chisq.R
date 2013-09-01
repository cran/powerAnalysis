#' Compute effect size of chi-squared test
#'
#' @param ct              a m x n Contingency Table (matrix with m rows and n columes)
#' @param chisq           the value the chi-squared test statistic
#' @param p               p value for the chi-squared test
#' @param df              degree of freedom (e.g., df=(m-1)*(n-1))
#' @param n               total number of observations (sample size)
#' @param mindf           the degrees of freedom for the variable with the smaller number of levels, if m > n, mindf=n-1, otherwise, mindf=m-1
#' @seealso               \code{\link{power.chisq}}
#' @return                es(effect size), chisq, p, df, n
#' @keywords              effect size, chi-squared test, study power, degree of freedom
#' @export
#' @examples
#' counts <- matrix(c(225,125,85,95),nrow=2,byrow=TRUE);ES.chisq(ct=counts)
#' @examples
#' case <- c(225,85,100);control <- c(125,95,125);counts <- cbind(case,control);ES.chisq(ct=counts)
#' @examples
#' p1 <- c(225,85,100);p2 <- c(125,95,125);p3 <- c(175,90,113);counts <- cbind(p1,p2,p3);ES.chisq(ct=counts)
#' @examples
#' ES.chisq(chisq=13.561,n=530,df=1,mindf=1)
#' @examples
#' ES.chisq(p=0.000231,n=530,df=1,mindf=1)
ES.chisq <- function(ct=NULL, chisq=NULL, p=NULL, n=NULL, df=NULL, mindf=NULL){
  es=NULL
  if(!is.null(ct)){
    row = nrow(ct)
    col = ncol(ct)
    for(i in 1:row){
      for(j in 1:col){
        if(ct[i,j]<1){
          stop("numbers in the contingency table must be at least 1\n")
        }
      }
    }
    es$df = (row -1)*(col -1)
    n = sum(ct)
    es$n=n
    mindf <- ifelse(row>col, col-1, row-1)
    chisq=chisq.test(ct)$statistic[[1]]
    es$chisq=chisq
    es$p=chisq.test(ct)$p.value
    es$es=sqrt(chisq/(n*mindf))
  }else if(!is.null(chisq)){
    if(chisq < 0){
      stop("chisq must be at least 0\n")
    }
    if(!is.null(n) && !is.null(df) && !is.null(mindf)){
      if(n<1){
        stop("total number of observations must be at least 1\n")
      }
      if(df < 1){
        stop("df must be at least 1\n")
      }
      if(mindf < 1){
        stop("mindf must be at least 1")
      }
      es$df=df
      es$n=n
      es$chisq=chisq
      es$p=pchisq(chisq,df,lower.tail=FALSE)
      es$es=sqrt(chisq/(n*mindf))
    }else{
      stop("n, df, and mindf are needed to calculate es(effect size)\n")
    }
  }else if(!is.null(p)){
    if(p < 0 || p > 1){
      stop("p must be in [0,1]\n")
    }
    if(!is.null(n) & !is.null(df) & !is.null(mindf)){
        if(n<1){
          stop("total number of observations must be at least 1\n")
        }
        if(df < 1){
          stop("df must be at least 1\n")
        }
        if(mindf < 1){
          stop("mindf must be at least 1\n")
        }
      es$df=df
      es$n=n
      chisq=qchisq(p,df,lower.tail=FALSE)
      es$chisq=chisq
      es$p=p
      es$es=sqrt(chisq/(n*mindf))
    }else{
      stop("n, df, and mindf are needed to calculate es(effect size)\n")
    }
  }else{
    stop("one of ct, chisq and p is needed to calculate es(effect size)\n")
  }
return(es)
}

