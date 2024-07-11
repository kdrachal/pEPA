\name{csc.test}
\alias{csc.test}

\title{Computes Test for Cross-Sectional Clusters.}

\description{
This function computes test of the equal predictive accuracy for cross-sectional clusters. The null hypothesis of this test is that a pair of forecasts have the same expected accuracy among cross-sectional clusters. However, their predictive accuracy can be different across the clusters, but the same among each cluster. The test is suitable if either: \eqn{K \ge 2} and \code{significance level} \eqn{\le 0.08326}, or \eqn{2 \le K \le 14} and \code{significance level} \eqn{\le 0.1}, or \eqn{K = \{ 2,3 \} } and \code{significance level} \eqn{\le 0.2}, where \eqn{K} denotes the number of time clusters.
}

\usage{
csc.test(evaluated1,evaluated2,realized,loss.type="SE",cl,dc=FALSE)
}

\arguments{
\item{evaluated1}{same as in \code{\link[pEPA]{pool_av.test}}, but cross-sections are ordered rowwise}
\item{evaluated2}{same as in \code{\link[pEPA]{pool_av.test}}, but cross-sections are ordered rowwise}
\item{realized}{same as in \code{\link[pEPA]{pool_av.test}}, but cross-sections are ordered rowwise}
\item{loss.type}{same as in \code{\link[pEPA]{pool_av.test}}}
\item{cl}{\code{\link[base]{vector}} of the beginning indices of rows for each pre-defined clusters -- as a result always \code{cl[1]=1}}
\item{dc}{\code{\link[base]{logical}} indicating if apply decorrelating clusters, if not specified \code{dc=FALSE} is used}
}

\value{class \code{htest} object, \code{\link[base]{list}} of 
\item{statistic}{test statistic}
\item{parameter}{\eqn{K}, number of cross-sectional clusters}
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
# 2 cross-sectional clusters: energy commodities and non-energy commodities
cs.cl <- c(1,9)
t <- csc.test(evaluated1=f.bsr,evaluated2=f.dma,realized=y,loss.type="SE",cl=cs.cl)
}

\references{
Qu, R., Timmermann, A., Zhu, Y. 2024. Comparing forecasting performance with panel data. \emph{International Journal of Forecasting} \bold{40}, 918--941. 
}

\seealso{\code{\link[pEPA]{pool_av.test}}}
