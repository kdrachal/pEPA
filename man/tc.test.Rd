\name{tc.test}
\alias{tc.test}

\title{Computes Test for Time Clusters.}

\description{
This function computes test of the equal predictive accuracy for time clusters. The null hypothesis of this test is that the equal predictive accuracy for the two methods holds within each of the time clusters. The test is suitable if either: \eqn{K \ge 2} and \code{significance level} \eqn{\le 0.08326}, or \eqn{2 \le K \leq 14} and \code{significance level} \eqn{\le 0.1}, or \eqn{K = \{ 2,3 \} } and \code{significance level} \eqn{\le 0.2}, where \eqn{K} denotes the number of time clusters.
}

\usage{
tc.test(evaluated1,evaluated2,realized,loss.type="SE",cl)
}

\arguments{
\item{evaluated1}{same as in \code{\link[pEPA]{pool_av.test}}}
\item{evaluated2}{same as in \code{\link[pEPA]{pool_av.test}}}
\item{realized}{same as in \code{\link[pEPA]{pool_av.test}}}
\item{loss.type}{same as in \code{\link[pEPA]{pool_av.test}}}
\item{cl}{\code{\link[base]{vector}} of the beginning indices of each pre-defined blocks of time -- as a result always \code{cl[1]=1}}
}

\value{class \code{htest} object, \code{\link[base]{list}} of 
\item{statistic}{test statistic}
\item{parameter}{\eqn{K}, number of time clusters}
\item{alternative}{alternative hypothesis of the test}
\item{p.value}{p-value}
\item{method}{name of the test}
\item{data.name}{names of the tested data}
}

\examples{
data(forecasts)
y <- t(observed)
f.bsr <- matrix(NA,ncol=ncol(y),nrow=56)
f.dma <- f.bsr
# extract prices predicted by BSR rec and DMA methods
for (i in 1:56)
  {
    f.bsr[i,] <- predicted[[i]][,1]
    f.dma[i,] <- predicted[[i]][,9]
  }
# 3 time clusters: Jun 1996 -- Nov 2007, Dec 2007 -- Jun 2009, Jul 2009 - Aug 2021
# rownames(observed)[1] 
# rownames(observed)[139] 
# rownames(observed)[158] 
t.cl <- c(1,139,158)
t <- tc.test(evaluated1=f.bsr,evaluated2=f.dma,realized=y,loss.type="SE",cl=t.cl)
}

\references{
Qu, R., Timmermann, A., Zhu, Y. 2024. Comparing forecasting performance with panel data. \emph{International Journal of Forecasting} \bold{40}, 918--941. 
}

\seealso{\code{\link[pEPA]{pool_av.test}}}
