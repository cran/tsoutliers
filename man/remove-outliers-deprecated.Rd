\name{remove.outliers-deprecated}
\alias{remove.outliers}

\title{Stage II of the Procedure: Discard Outliers}

\description{
As of version 0.6-6, \code{remove.outliers} has been renamed as \code{discard.outliers}.
\code{\link{discard.outliers}} should be used.
}

\usage{
remove.outliers(x, y, cval = NULL, 
  method = c("en-masse", "bottom-up"), 
  delta = 0.7, tsmethod.call = NULL, 
  fdiff = NULL, logfile = NULL)
}

\arguments{
\item{x}{a list. The output returned by \code{\link{locate.outliers.oloop}}.}
\item{y}{a time series.}
\item{cval}{a numeric. The critical value to determine the significance of each type of outlier.}
\item{method}{a character. The method to discard/remove outliers. See details.}
\item{delta}{a numeric. Parameter of the temporary change type of outlier.}
\item{tsmethod.call}{an optional \code{call} object. The call to the function used to fit the time 
series model.}
\item{fdiff}{currently ignored.}
\item{logfile}{a character or \code{NULL}. It is the path to the file where 
tracking information is printed. Ignored if \code{NULL}.}
}

\value{
A list.
}

\seealso{
\code{\link{discard.outliers}}.
}

\keyword{ts}
