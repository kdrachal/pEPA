\name{pool_av.test}
\alias{pool_av.test}

\title{Computes Test for the Pooled Average.}

\description{
This function computes test of the equal predictive accuracy for the pooled average. The null hypothesis of this test is that the pooled average loss is equal in expectation for a pair of forecasts from both considered methods. The alternative hypothesis can be formulated as the differences do not average out across the cross-sectional and time-series dimensions. 
}

\usage{
pool_av.test(evaluated1,evaluated2,realized,loss.type="SE",J=NULL)
}

\arguments{
\item{evaluated1}{\code{\link[base]{matrix}} of forecasts from the first method, cross-sections are ordered by rows, and time by columns}
\item{evaluated2}{\code{\link[base]{matrix}} of forecasts from the second method, cross-sections are ordered by rows, and time by columns}
\item{realized}{\code{\link[base]{matrix}} of the observed values, cross-sections are ordered by rows, and time by columns}
\item{loss.type}{a method to compute the loss function, \code{loss.type="SE"} applies squared errors, \code{loss.type="AE"} -- absolute errors, \code{loss.type="SPE"} -- squared proportional error (useful if errors are heteroskedastic), \code{loss.type="ASE"} -- absolute scaled error, if \code{loss.type} is specified as some \code{\link[base]{numeric}}, then the function of type \code{exp(loss.type*errors)-1-loss.type*errors} is applied (useful when it is more costly to underpredict \code{realized} than to overpredict), if not specified \code{loss.type="SE"} is used}
\item{J}{\code{\link[base]{numeric}} maximum lag length, if not specified \code{J=round(T^(1/3))} is used, where \code{T=ncol(realized)}}
}

\value{class \code{htest} object, \code{\link[base]{list}} of 
\item{statistic}{test statistic}
\item{parameter}{\code{J}, maximum lag length}
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
t <- pool_av.test(evaluated1=f.bsr,evaluated2=f.dma,realized=y,loss.type="SE")
}

\references{
Hyndman, R.J., Koehler, A.B. 2006. Another look at measures of forecast accuracy. \emph{International Journal of Forecasting} \bold{22}, 679--688.

Qu, R., Timmermann, A., Zhu, Y. 2024. Comparing forecasting performance with panel data. \emph{International Journal of Forecasting} \bold{40}, 918--941. 

Taylor, S. J., 2005. \emph{Asset Price Dynamics, Volatility, and Prediction},  Princeton University Press.

Triacca, U., 2024. \emph{Comparing Predictive Accuracy of Two Forecasts}, \url{https://www.lem.sssup.it/phd/documents/Lesson19.pdf}.
}
