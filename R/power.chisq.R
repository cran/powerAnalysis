#' Compute study power of chi-squared test
#'
#' @param es              effect size. A numeric value or output of ES.chisq
#' @param df              degree of freedom
#' @param n               total number of observations (sample size)
#' @param power           power of study
#' @param sig.level       significance level
#' @seealso               \code{\link{ES.chisq}}
#' @return                power, es, df, n, sig.level
#' @keywords              study power, effect size, chi-squared test, sample size, significant level
#' @export
#' @examples
#' counts <- matrix(c(225,125,85,95),nrow=2,byrow=TRUE);
#' power.chisq(es=ES.chisq(ct=counts),sig.level=0.05)
#' @examples
#' power.chisq(es=0.16,df=1,n=530,sig.level=0.05)
#' @examples
#' power.chisq(es=0.16,df=1,n=530,power=0.9576)
#' @examples
#' power.chisq(es=0.16,df=1,power=0.9576,sig.level=0.05)
#' @examples
#' power.chisq(df=1,n=530,power=0.9576,sig.level=0.05)
power.chisq <- function(es=NULL,df=NULL,n=NULL,power=NULL,sig.level=NULL){
  pwr=NULL
  if(is.list(es) & length(es)==5){
    myes=es$es
    df=es$df
    n=es$n
    pwr$es=es$es
    pwr$df=es$df
    pwr$n=es$n
    if (sum(sapply(list(power, sig.level), is.null)) != 1){
      stop("exactly one of power and sig.level must be NULL\n")
    }else if(is.null(power)){
      temp <- qchisq(sig.level, df = df, lower.tail = FALSE)
      pwr$power=pchisq(temp, df = df, ncp = n * myes^2, lower.tail = FALSE)
      pwr$sig.level=sig.level
    }else{
      pwr$power=power
      func1 <- function(sig.level){
        temp <- qchisq(sig.level, df = df, lower.tail = FALSE)
        pchisq(temp, df = df, ncp = n * myes^2, lower.tail = FALSE) - power
      }
      pwr$sig.level=uniroot(func1,lower=1e-15,upper=1-1e-15)$root
    }
  }else{
    if(!is.numeric(df) || df<1){
      stop("df must be at least 1\n")
    }
    if(sum(sapply(list(es, n, power, sig.level), is.null)) != 1){
      stop("exactly one of es, n, power, and sig.level must be NULL\n")
    }
    if(!is.null(es) && (!is.numeric(es) || es<0)){
      stop("es must be positive\n")
    }
    if(!is.null(n) && (!is.numeric(n) || n<1)){
      stop("n must be at least 1\n")
    }
    if(!is.null(power) && (!is.numeric(power) || power>1 || power<0)){
      stop("power must be in [0,1]\n")
    }
    if(!is.null(sig.level) && (!is.numeric(sig.level) || sig.level>1 || sig.level<0)){
      stop("sig.level must be in [0,1]\n")
    }
    
    if(is.null(power)){
      pwr$es=es
      pwr$df=df
      pwr$n=n
      temp <- qchisq(sig.level, df = df, lower.tail = FALSE)
      pwr$power=pchisq(temp, df = df, ncp = n * es^2, lower.tail = FALSE)
      pwr$sig.level=sig.level
    }else if(is.null(n)){
      pwr$es=es
      pwr$df=df
      func2 <- function(n){
        temp <- qchisq(sig.level, df = df, lower.tail = FALSE)
        pchisq(temp, df = df, ncp = n * es^2, lower.tail = FALSE) - power
      }
      pwr$n=uniroot(func2,lower=1+1e-15,upper=1e+15)$root
      pwr$power=power
      pwr$sig.level=sig.level
    }else if(is.null(es)){
      func3 <- function(es){
        temp <- qchisq(sig.level, df = df, lower.tail = FALSE)
        pchisq(temp, df = df, ncp = n * es^2, lower.tail = FALSE) - power
      }
      pwr$es=uniroot(func3,lower=1e-15,upper=1-1e-15)$root
      pwr$df=df
      pwr$n=n
      pwr$power=power
      pwr$sig.level=sig.level
    }else if(is.null(sig.level)){
      pwr$es=es
      pwr$df=df
      pwr$n=n
      pwr$power=power
      func4 <- function(sig.level){
        temp <- qchisq(sig.level, df = df, lower.tail = FALSE)
        pchisq(temp, df = df, ncp = n * es^2, lower.tail = FALSE) - power
      }
      pwr$sig.level=uniroot(func4,lower=1e-15,upper=1-1e-15)$root
    }
  }
return(pwr)
}
