\name{calendar.effects}
\alias{calendar.effects}

\title{Calendar Effects}

\description{
This function creates regressor variables for trading day, Easter and leap year effects 
over the sample period where the input time series is sampled.
}

\usage{
calendar.effects(x, trading.day = TRUE, easter = 6, 
  leap.year = FALSE, holidays = NULL, easter.date = FALSE)
}

\arguments{
\item{x}{a monthly time series.}
\item{trading.day}{logical. If \code{TRUE}, trading day regressor variable is returned.}
\item{easter}{numeric. The number of days before Easter over which the Easter effect spans. 
If it is set to zero Easter variable is not returned.}
\item{leap.year}{logical. If \code{TRUE}, leap year regressor variable is returned.}
\item{holidays}{an optional numeric vector of the same length as \code{x} containing 
the number of holidays for each time period.}
\item{easter.date}{logical indicating whether the date of Easter should be returned.}
}

\details{
Let \eqn{wd} be the number of working days in a given month 
and \eqn{nwd} the number of non-working days.
The trading day variable at time \eqn{t} is built as follows:

\deqn{wd = wd - holidays}
\deqn{nwd = nwd + holidays}
\deqn{td_t = wd - (5/2) \times nwd}{td_t = wd - (5/2) x nwd}

By default, working days are the days from Monday to Friday and 
non-working days are Saturdays and Sundays. If there are additional 
non-working days they can be defined in argument \code{holidays}. 
For example, if the 1st of February is a local holiday, the user can define 
a variable containing zeros for all periods except for the periods related 
to February where the 1st of February falls within a working day (Monday to Friday); 
these data are set to one so that they are considered as non working days.

Easter effect is defined as the proportion of days before Easter 
(by default \code{easter = 6}) that lie in March and April, 
respectively for each month. It contains zeros for the remaining months.

The leap year is a vector of zeros for all months except February, 
where the variable takes on the value \eqn{0.75} if the year is a leap year
and \eqn{-0.25} otherwise.
}

\value{
A \code{mts} matrix containing the selected calendar effects by columns.

If \code{easter.date} is \code{TRUE} a list is returned containing the 
\code{mts} matrix of calendar effects as well as the dates of Easter for 
each year of the sample.
}

\examples{
# display calendar effects for a sample span period
# no data are actually necessary in the input series 
# since calendar effects are concerned only with the dates
# at which the data are sampled
x <- ts(frequency = 12, start = c(1980, 1), end = c(2000, 12))
ce <- calendar.effects(x, leap.year = TRUE)
colnames(ce)
plot(ce, main = "calendar effects")
# Easter days for each year
calendar.effects(x, easter.date = TRUE)$easter
}

\keyword{ts}
