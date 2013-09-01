#' Sample size estimation of chi-squared test
#' 
#' @param es              effect size.
#' @param power           power of study
#' @param df              degree of freedom
#' @param sig.level       significance level
#' @param allele          in genetic association study, whether test allele or genotype
#' @param xlab            a title for the x axis
#' @param ylab            a title for the y axis
#' @param main            an overall title for the plot
#' @param grid            add grid lines or not
#' @param ...             Arguments to be passed to methods
#' @seealso               \code{\link{power.chisq}}
#' @return                power, es, df, n, sig.level
#' @keywords              study power, effect size, chi-squared test, sample size, significant level
#' @export
#' @examples
#' es=seq(from=0.05,to=0.5,by=0.05);
#' power=seq(from=0.7,to=0.9,by=0.1);
#' power.plot.chisq(es=es,power=power,df=1,sig.level=0.05,allele=TRUE)
power.plot.chisq <- function(es=NULL,power=NULL,df=NULL,sig.level=NULL,allele=FALSE,xlab="Effect Size", ylab="Sample Size",main="Sample Size Estimation",grid=FALSE,...){
  rr=NULL
  nes=length(es)
  npower=length(power)
  
  n <- matrix(rep(0,nes*npower), nrow=nes,byrow=T)
  for (i in 1:npower){
    for (j in 1:nes){
      result <- power.chisq(es=es[j],df=df,n=NULL,power=power[i],sig.level=sig.level)
      if(allele==TRUE){
        n[j,i] <- ceiling(result$n/2)
      }else{
        n[j,i] <- ceiling(result$n)
      }
    }
  }
  
  xrange <- range(es)
  yrange <- round(range(n))
  colors <- rainbow(length(power))
  plot=plot(xrange, yrange, type="n",xaxt="n",xlab=xlab,ylab=ylab,main=main)
  
  for (i in 1:npower){
    lines(es, n[,i], type="l", lwd=2, col=colors[i])
  }
  
  if(grid==TRUE){
    abline(v=0, h=seq(0,yrange[2],length.out=21), lty=2, col="grey89")
    abline(h=0, v=seq(xrange[1],xrange[2],0.01), lty=2,col="grey89")
  }

  axis(side=1,at=es,labels=es)
  legend("topright", title="power", as.character(power), fill=colors)
  
  rownames=paste("effect.size=",es,sep="")
  colnames=paste("power=",power,sep="")
  dimnames(n)=list(rownames,colnames)
  rr$sig.level=sig.level
  rr$df=df
  rr$sample.size=n
  return(rr)
}

