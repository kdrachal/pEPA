\encoding{UTF-8}

\name{observed}
\alias{observed}

\title{Sample Panel of Commodities Spot Prices.}

\description{Observed spot prices of various commodities.}

\usage{data(forecasts)}

\format{
\itemize{\code{observed} is \code{\link[base]{matrix}} object such that its columns correspond to spot prices of selected 56 commodities.
}
}


\details{
They cover the period between 1996 and 2021, and are in monthly freqency. Variables names are the same as in the paper by Drachal and Pawłowski (2024). The observed prices were taken from The World Bank (2022).  
}

\examples{
data(forecasts)
# WTI prices
t1 <- observed[,3]
}

\references{
Drachal, K., Pawłowski, M. 2024. Forecasting selected commodities' prices with the Bayesian symbolic regression. \emph{International Journal of Financial Studies} \bold{12}, 34, \doi{10.3390/ijfs12020034}

The World Bank. 2022. Commodity Markets. \url{https://www.worldbank.org/en/research/commodity-markets} 
}

\seealso{\code{\link[pEPA]{predicted}}}
