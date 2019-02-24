\name{locate.outliers.loops}
\alias{locate.outliers.oloop}
\alias{locate.outliers.iloop}

\title{Stage I of the Procedure: Locate Outliers (Loop Around Functions)}

\description{
These functions implement the inner and outer loops based of the procedure 
to locate outliers following the approach described in Chen and Liu (1993) . 
}

\usage{
locate.outliers.oloop(y, fit, types = c("AO", "LS", "TC"), 
  cval = NULL, maxit.iloop = 4, maxit.oloop = 4, 
  delta = 0.7, logfile = NULL)
locate.outliers.iloop(resid, pars, cval = 3.5, 
  types = c("AO", "LS", "TC"), maxit = 4, delta = 0.7, logfile = NULL)
}

\arguments{
\item{y}{a time series.}
%
\item{fit}{an \code{Arima} object. The output from \code{\link[stats]{arima}} or 
\code{\link[forecast]{auto.arima}}}.
%
\item{resid}{a time series. Residuals from a time series model fitted to the data.}
%
\item{pars}{a list containing the parameters of the model fitted to the data. 
See details in \code{\link{locate.outliers}}.}
%
\item{cval}{a numeric. The critical value to determine the significance of each type of outlier.}
%
\item{types}{a character vector indicating the type of outlier to be considered by the 
detection procedure among the following: 
innovational outliers (\code{"IO"}), additive outliers (\code{"AO"}), 
level shifts (\code{"LS"}), temporary changes (\code{"TC"}) 
and seasonal level shifts (\code{"SLS"}).}
%
\item{maxit}{a numeric. The maximum number of iterations in the inner loop.}
%
\item{maxit.iloop}{a numeric. Same as argument as \code{maxit} to be passed to 
\code{locate.outliers.iloop}.}
%
\item{maxit.oloop}{a numeric. The maximum number of iterations in the outer loop.}
%
\item{delta}{a numeric. Parameter of the temporary change type of outlier.}
%
\item{logfile}{a character or \code{NULL}. It is the path to the file where 
tracking information is printed. Ignored if \code{NULL}.}
}

\details{
See also the details section in \code{\link{locate.outliers}}.
%
%Five types of outliers can be considered. 
%By default: \code{"AO"} additive outliers, \code{"LS"} level shifts,
%and \code{"TC"} temporary changes are selected;
%\code{"IO"} innovative outliers and \code{"SLS"} seasonal level shifts 
%can also be selected.

The function \code{\link{locate.outliers.iloop}} iterates around the function 
\code{locate.outliers} until no additional outliers are found or the maximum number 
of iterations is reached. 
After each iteration, the effect of the outliers on the residuals of the fitted model 
is removed and the \eqn{t}-statistics are obtained again for the modified residuals.
No model selection or refit of the model is conducted within this loop.

The function \code{locate.outliers.oloop} is the outer loop of the procedure to locate 
outliers. It iterates around the function \code{\link{locate.outliers.iloop}}.
At the end of each iteration the detected outliers are removed from the original 
data. Then, the time series model is fitted (or selected) again for the adjusted series 
and a new search for outliers is executed. The outer loop stops when no additional outliers 
are detected.

In function \code{locate.outliers.oloop},
if no value is specified for argument \code{cval} a default value based on the sample size 
is used. Let \eqn{n} be the number of observations. If \eqn{n \le 50} then \code{cval} is set 
equal to \eqn{3.0}; If \eqn{n \ge 450} then \code{cval} is set equal to \eqn{4.0};
otherwise \code{cval} is set equal to \eqn{3 + 0.0025 * (n - 50)}.
}

\value{
\code{locate.outliers.iloop} returns a data frame defining by rows each 
detected outlier. The data frame follows the same format as the output from 
\code{\link{locate.outliers}}.

\code{locate.outliers.oloop} returns a list containing the following elements:
\code{fit}: information from the last fitted model that will be required by 
other functions in the automatic procedure (parameter estimates, residuals and 
number of observations);
\code{outliers}: a data frame defining by rows the detected outliers;
\code{iter}: the number of iterations employed by the outer loop.
}

\note{
In \code{locate.outliers.iloop} 
the default critical value, \code{cval}, is set equal to \eqn{3.5} and, 
hence, it is not based on the sample size.
\code{locate.outliers.oloop} uses a default critical value based on the 
sampel size as in \code{\link{tso}}.  
}

\seealso{
\code{\link{outliers.tstatistics}}, \code{\link{tso}}.
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

Gómez, V. and Taguas, D. (1995). 
\emph{Detección y Corrección Automática de Outliers con TRAMO: 
Una Aplicación al IPC de Bienes Industriales no Energéticos}.
Ministerio de Economía y Hacienda. Document number D-95006.
\url{http://www.sepg.pap.minhap.gob.es/sitios/sepg/es-ES/Presupuestos/Documentacion/Documents/DOCUMENTOS\%20DE\%20TRABAJO/D95006.pdf}

Kaiser, R., and Maravall, A. (1999). 
\emph{Seasonal Outliers in Time Series}.
Banco de España, Servicio de Estudios. Working paper number 9915.
}

\examples{
# additional outliers may be detected in the inner or outlier loops
# in this case, the inner does not find further potential outliers 
# and stops in the first iteration, while the outer loop detects 
# a new outlier
data("hicp")
y <- log(hicp[["011600"]])
fit <- arima(y, order = c(1, 1, 0), seasonal = list(order = c(2, 0, 2)))
resid <- residuals(fit)
pars <- coefs2poly(fit)
otypes <- c("AO", "LS", "TC")
mo0 <- locate.outliers(resid, pars, types = otypes)
mo0
mo1 <- locate.outliers.iloop(resid, pars, types = otypes)
mo1
mo2 <- locate.outliers.oloop(y, fit, types = otypes)
mo2$iter
mo2$outliers
}

\keyword{ts}
