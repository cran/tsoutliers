\name{outliers.effects}
\alias{outliers.effects}

\title{Create the Pattern of Different Types of Outliers}

\description{
These functions create a unit or weighted impulse for the five types 
of outliers considered in the package.
}

\usage{
outliers.effects(mo, n, weights = FALSE, delta = 0.7, 
  pars = NULL, freq = 12)
}

\arguments{
\item{mo}{a data frame defining the type, location and weight of the outliers to be created.}
%
\item{n}{a numeric. The length of the variable that will contain the outlier.}
%
\item{weights}{logical. If \code{TRUE}, the outliers are weighted by the 
values in column \code{"coefhat"} of the data frame \code{mo}. Otherwise, 
unit weights are considered.}
%
\item{delta}{a numeric. Parameter of the temporary change type of outlier.}
%
\item{pars}{a list containing the parameters of the time series model fitted to the data. 
Only for innovational outlier. See the details section in \code{\link{locate.outliers}}.}
%
\item{freq}{a numeric. The periodicity of the data. Only for seasonal level shift.}
}

\details{
These functions delineate the effect of each type of outlier on the observed data. 
See the example below for a representation of the outliers.

The function \code{outliers.effects} operates directly on the output returned 
by \code{\link[stats]{time}}. The remaining functions are called by 
\code{outliers.effects} and can be used as a simpler interface to define some outliers.
They generate the type of outliers indicated by their names.

The column names of the data frame \code{mo} must follow the same convention as 
the output returned by \code{\link{locate.outliers}}: the column containing the 
observation at which the outlier sparks is named \code{"ind"}; 
the column containing the weights is named \code{"coefhat"} and the type 
of outlier is specified in a column named \code{"type"}. The column \code{"ind"} 
should contain the index time point, not the time point in terms of year and season 
as given by \code{\link[stats]{time}}.
}

\value{
A \code{n} \eqn{\times} \code{nrow(mo)} or \code{n} \eqn{\times} \code{length(ind)} matrix
containing by columns each outlier.
}

\seealso{
\code{\link{locate.outliers}}, \code{\link{remove.outliers}}, 
\code{\link{tso}}.
}

\references{
Chen, C. and Liu, Lon-Mu (1993). 
\sQuote{Joint Estimation of Model Parameters and Outlier Effects in Time Series}.
\emph{Journal of the American Statistical Association},
\bold{88}(421), pp. 284-297.

Gómez, V. and Maravall, A. (1996). 
\emph{Programs TRAMO and SEATS. Instructions for the user}.
Banco de España, Servicio de Estudios. Working paper number 9628.
\url{http://www.bde.es/f/webbde/SES/Secciones/Publicaciones/PublicacionesSeriadas/DocumentosTrabajo/96/Fich/dt9628e.pdf}

Kaiser, R., and Maravall, A. (1999). 
\emph{Seasonal Outliers in Time Series}.
Banco de España, Servicio de Estudios. Working paper number 9915.
\url{http://www.bde.es/f/webbde/SES/Secciones/Publicaciones/PublicacionesSeriadas/DocumentosTrabajo/99/Fic/dt9915e.pdf}
}

\examples{
n <- 30
# innovative outlier based on ARMA(3, 2) model
mo <- outliers("IO", 10)
io <- outliers.effects(mo, n, pars = list(arcoefs = c(0.8, -0.6, 0.2), 
  macoefs = c(-0.5, 0.2)))
plot(c(io[seq.int(10)], rep(NA, 20)), type = "s", ylim = range(io), 
  ylab = "io", main = "IO based on ARMA(3,2)")
lines(c(rep(NA, 9), io[-seq.int(9)]))

# innovative outlier based on Airlines model ARIMA(0,1,1)(0,1,1)
p1 <- c(1, -1)
p2 <- c(1, rep(0, 3), -1)
p1xp2 <- c(1, -1, 0, 0, -1, 1)
p1b <- c(1, -0.6)
p2b <- c(1, rep(0, 3), -0.6)
p2bxp2b <- c(1, -0.6, 0, 0, -0.6, 0.36)
io2 <- outliers.effects(mo, n, pars = list(arcoefs = -p1xp2[-1], 
  macoefs = p2bxp2b[-1]))
plot(c(io2[seq.int(10)], rep(NA, 20)), type = "s", ylim = range(io2), 
  main = "IO based on ARIMA(0,1,1)(0,1,1)", ylab = "io2")
lines(c(rep(NA, 9), io2[-seq.int(9)]))

# additive outlier
mo <- outliers("AO", 10)
ao <- outliers.effects(mo, n)
plot(ao, type = "h", main = "AO: additive outlier")

# level shift
mo <- outliers("LS", 10)
ls <- outliers.effects(mo, n)
plot(ls, type = "s", main = "LS: level shift")

# temporary change
mo <- outliers("TC", 10)
tc <- outliers.effects(mo, n)
plot(c(tc[seq.int(10)], rep(NA, 20)), type = "s", 
  main = "TC: temporary change", ylab = "tc")
lines(c(rep(NA, 9), tc[-seq.int(9)]))

# the temporary change with parameter 0.7 is equivalent to 
# an IO for an AR(1) model with coefficient 0.7
mo <- outliers("IO", 10)
io3 <- outliers.effects(mo, n = n, pars = list(arcoefs = c(0.7), macoefs = c(0)))
all.equal(io3, tc, check.attributes=FALSE)

# the temporary change with parameter 0.7 is equivalent to 
# an IO for an AR(1) model with coefficient 0.7
mo <- outliers("IO", 10)
io4 <- outliers.effects(mo, n = n, pars = list(arcoefs = 1, macoefs = 0))
all.equal(io4, ls, check.attributes=FALSE)

# seasonal level shift (quarterly data)
mo <- outliers("SLS", 10)
sls <- outliers.effects(mo, n, freq = 4)
plot(sls, type = "h", main = "SLS: seasonal level shift")
}

\keyword{ts}
