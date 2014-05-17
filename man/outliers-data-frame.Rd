\name{outliers}
\alias{outliers}

\title{Define Outliers in a Data Frame}

\description{
This function is an interface to create a data frame defining 
the type, observation and weight of outliers. The output is properly designed to 
be used as input to other functions such as \code{\link{outliers.effects}} 
and \code{\link{outliers.regressors}}.
}

\usage{
outliers(type, ind, weight = NULL)
}

\arguments{
\item{type}{a character vector. The types of outliers ("IO", "AO", "LS", "TC", "SLS").}
\item{ind}{a numeric vector. The observations at which each outlier takes effect.}
\item{weight}{an optional numeric vector. The weights of the outliers. If \code{NULL} the 
weights are set equal to one.}
}

\value{
A data frame.
}

\examples{
outliers(c("AO", "LS"), c(10, 20))
outliers(c("AO", "LS"), c(10, 20), c(2.1, 3.2))
}

\keyword{ts}
