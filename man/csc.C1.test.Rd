\name{csc.C1.test}
\alias{csc.C1.test}

\title{Computes Test for Cross-Sectional Clusters.}

\description{
This function computes test of the equal predictive accuracy for cross-sectional clusters. It corresponds to \eqn{C^{(1)}_{nT}} statistic in the referenced paper by Akgun et al. (2024). The null hypothesis of this test is that a pair of forecasts have the same expected accuracy among cross-sectional clusters. However, their predictive accuracy can be different across the clusters, but the same among each cluster. The test is suitable for situations with cross-sectional independence. 
}

\usage{
csc.C1.test(evaluated1,evaluated2,realized,loss.type="SE",cl)
}

\arguments{
\item{evaluated1}{same as in \code{\link[pEPA]{pool_av.test}}, but cross-sections are ordered rowwise}
\item{evaluated2}{same as in \code{\link[pEPA]{pool_av.test}}, but cross-sections are ordered rowwise}
\item{realized}{same as in \code{\link[pEPA]{pool_av.test}}, but cross-sections are ordered rowwise}
\item{loss.type}{same as in \code{\link[pEPA]{pool_av.test}}}
\item{cl}{\code{\link[base]{vector}} of the beginning indices of rows for each pre-defined clusters -- as a result always \code{cl[1]=1}}
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
\donttest{
data(forecasts)
y <- t(observed)
# just to save time
y <- y[,1:40]
f.bsr <- matrix(NA,ncol=ncol(y),nrow=56)
f.dma <- f.bsr
# extract prices predicted by BSR rec and DMA methods
for (i in 1:56)
  {
    f.bsr[i,] <- predicted[[i]][1:40,1]
    f.dma[i,] <- predicted[[i]][1:40,9]
  }
# 2 cross-sectional clusters: energy commodities and non-energy commodities
cs.cl <- c(1,9)
t <- csc.C1.test(evaluated1=f.bsr,evaluated2=f.dma,realized=y,loss.type="SE",cl=cs.cl)
}
}

\references{
Akgun, O., Pirotte, A., Urga, G., Yang, Z. 2024. Equal predictive ability tests based on panel data with applications to OECD and IMF forecasts. \emph{International Journal of Forecasting} \bold{40}, 202--228.
}

\seealso{\code{\link[pEPA]{pool_av.test}}, \code{\link[pEPA]{csc.C3.test}}}
