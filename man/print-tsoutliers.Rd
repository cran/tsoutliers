\name{print.tsoutliers}
\alias{print.tsoutliers}

\title{Print \code{tsoutliers} object}

\description{
This function prints the output from function \code{\link{tso}}. 
}

\usage{
\method{print}{tsoutliers}(x, digits = max(3L, getOption("digits") - 3L), call = FALSE, ...)
}

\arguments{
\item{x}{a list of class \code{tsoutliers} as returned by \code{\link{tso}}.}
\item{digits}{integer, the number of significant digits to be used.}
\item{call}{logical, if \code{TRUE} the \code{call} to the fitting method is printed.}
\item{...}{further arguments, currently ignored.}
}

\details{
When \code{\link[stats]{arima}} is run from \code{\link[base]{do.call}}, the latter 
generates a large object for the \code{call} which includes the entire structure 
of the data. This often implies a long output which may not be desired. 
Setting \code{call=FALSE} avoids printint the \code{call}.
}

\value{
None.
}

\seealso{
\code{\link{tso}}.
}

\keyword{print}
