\name{plot.tsoutliers}
\alias{plot.tsoutliers}

\title{Display Outlier Effects Detected by \code{tsoutliers}}

\description{
This function displays the output from function \code{\link{tso}}. 
}

\usage{
\method{plot}{tsoutliers}(x, 
  args.lines.y = list(col = "gray80"), args.lines.yadj = list(col = "blue"),
  args.lines.effects = list(type = "s", col = "red"),   
  args.points = list(col = "gray80", bg = "red", pch = 21), plot.points = TRUE, 
  args.x.axis = list(at = pretty(time(x$y)), tcl = -0.5, lwd = 0, lwd.ticks = 1),
  args.y.axis = list(at = pretty(x$y), tcl = -0.5, lwd = 0, lwd.ticks = 1),
  args.effects.axis = list(at = pretty(x$effects), tcl = -0.5, lwd = 0, lwd.ticks = 1),
  ...)
}

\arguments{
\item{x}{a list of class \code{tsoutliers} as returned by \code{\link{tso}}.}
\item{args.lines.y}{a list. Arguments passed to \code{\link[graphics]{lines}} to 
customize the line displaying the original series.}
\item{args.lines.yadj}{a list. Arguments passed to \code{\link[graphics]{lines}} to 
customize the line displaying the series adjusted for outliers effects.}
\item{args.lines.effects}{a list. Arguments passed to \code{\link[graphics]{lines}} to 
customize the line displaying the for outliers effects.}
\item{args.points}{a list. Arguments passed to \code{\link[graphics]{lines}} to 
customize the points drawn if \code{plot.points = TRUE}.}
\item{plot.points}{a logical indicating whether the time points of the outliers should 
be drawn as points over the original series.}
\item{args.x.axis}{a list. Arguments to be passed to \code{\link[graphics]{axis}} to customize
the x-axis (time).}
\item{args.y.axis}{a list. Arguments to be passed to \code{\link[graphics]{axis}} to customize 
the y-axis for the original series.}
\item{args.effects.axis}{a list. Arguments to be passed to \code{\link[graphics]{axis}} to customize 
the y-axis for the outliers effects.}
\item{...}{further arguments to be passed to \code{\link[graphics]{par}}.}
}

\details{
Instead of using the ellipsis, \code{...}, arguments passed to other functions are 
defined by means of a list. This approach is taken because there may be a single 
argument name to be used in different parts of the plot with a different value.
For example, the argument \code{col} can be defined in \code{args.lines.y} to indicate 
the color of the original series, e.g. \code{col = "gray80"}; at the same time 
the color for the adjusted series can be defined in the list argument \code{args.lines.yadj}.

For further customizations, the source code of the function can be modified relatively 
easy. Alternatively, a similar plot can be displayed simply as:
\code{plot(cbind(x$y, x$yadj, x$effects), plot.type = "multiple")}. 
In this way, the plot can be fully customized by setting the desired arguments to 
to \code{plot} or to ancillary functions that can be called afterwards.
}

\value{
None.
}

\seealso{
\code{\link{tso}}.
}

\keyword{hplot}
