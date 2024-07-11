\name{pool_av.S3.test}
\alias{pool_av.S3.test}

\title{Computes Test for Overall Equal Predictive Ability.}

\description{
This function computes test of the equal predictive accuracy for the pooled average. It corresponds to \eqn{S^{(3)}_{nT}} statistic in the referenced paper by Akgun et al. (2024). The null hypothesis of this test is that the pooled average loss is equal in expectation for a pair of forecasts from both considered methods. The alternative one is that the differences do not average out across the cross-sectional and time-series dimensions. The test allows for strong cross-sectional dependence.
}

\usage{
pool_av.S3.test(evaluated1,evaluated2,realized,loss.type="SE")
}

\arguments{
\item{evaluated1}{same as in \code{\link[pEPA]{pool_av.test}}}
\item{evaluated2}{same as in \code{\link[pEPA]{pool_av.test}}}
\item{realized}{same as in \code{\link[pEPA]{pool_av.test}}}
\item{loss.type}{same as in \code{\link[pEPA]{pool_av.test}}}
}

\value{class \code{htest} object, \code{\link[base]{list}} of 
\item{statistic}{test statistic}
\item{alternative}{alternative hypothesis of the test}
\item{p.value}{p-value}
\item{method}{name of the test}
\item{data.name}{names of the tested data}
}

\examples{
\donttest{
data(forecasts)
y <- t(observed)
# just to reduce computation time shorten time-series
y <- y[,1:40]
f.bsr <- matrix(NA,ncol=ncol(y),nrow=56)
f.dma <- f.bsr
# extract prices predicted by BSR rec and DMA methods
for (i in 1:56)
  {
    f.bsr[i,] <- predicted[[i]][1:40,1]
    f.dma[i,] <- predicted[[i]][1:40,9]
  }
t <- pool_av.S3.test(evaluated1=f.bsr,evaluated2=f.dma,realized=y,loss.type="SE")
}
}

\references{
Akgun, O., Pirotte, A., Urga, G., Yang, Z. 2024. Equal predictive ability tests based on panel data with applications to OECD and IMF forecasts. \emph{International Journal of Forecasting} \bold{40}, 202--228.
}

\seealso{\code{\link[pEPA]{pool_av.test}}, \code{\link[pEPA]{pool_av.S1.test}}}
