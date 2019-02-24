\name{outliers.tstatistics}
\alias{outliers.tstatistics}
%\alias{tstats.arima}
%\alias{tstas.stsm}

\title{Test Statistics for the Significance of Outliers}

\description{
This function computes the \eqn{t}-statistics to assess the significance 
of different types of outliers at every possible time point. 
The statistics can be based either on an ARIMA model, 
\code{\link[stats]{arima}} or \code{\link[forecast]{auto.arima}}.
}

\usage{
outliers.tstatistics(pars, resid, types = c("AO", "LS", "TC"), 
  sigma = NULL, delta = 0.7)
}

\arguments{
\item{pars}{a list containing the parameters of the model. 
See details section in \code{\link{locate.outliers}}.}
\item{resid}{a time series. Residuals of the ARIMA model fitted to the data.}
\item{types}{a character vector indicating the types of outliers to be considered.}
\item{sigma}{a numeric or \code{NULL}. Standard deviation of residuals.}
\item{delta}{a numeric. Parameter of the temporary change type of outlier.}
%\item{n.start}{a numeric. The number of warming observations added to the 
%input passed to the Kalman filter. Only for \code{pars} of class \code{"stsmSS"}.}
%\item{tsmethod}{an optional character indicating whether argument \code{pars} 
%contains the parameters from an ARIMA model (\code{tsmethod = "arima"} or from a 
%structural time series model (\code{tsmethod = "stsm"}). If it is kept as
%\code{NULL} its value is deduced from the elements passed in \code{pars}.}
}

\details{
%\deqn{
%\begin{array}{lll}
%\hbox{IO:}\quad L(B) = \frac{\theta(B)}{\alpha(B)\phi(B)} \,; & 
%\displaystyle
%\hbox{LS:}\quad L(B) = \frac{1}{(1 - B)} \,; &
%\displaystyle
%\hbox{SLS:}\quad L(B) = \frac{1}{(1 - B^s)} \,; \\ [2ex]
%\displaystyle
%\hbox{AO:}\quad L(B) = 1 \,; &
%\displaystyle
%\hbox{TC:}\quad L(B) = \frac{1}{(1 - \delta B)} \,.
%\end{array}
%}{
%IO: L(B) = \theta(B) / (\alpha(B) \phi(B)) \br
%LS: L(B) = 1 / (1 - B) \br
%SLS: L(B) = 1 / (1 - B^s) \br
%AO: L(B) = 1 \br
%TC: L(B) = 1 / (1 - \delta B) \br
%}
Five types of outliers can be considered. 
By default: \code{"AO"} additive outliers, \code{"LS"} level shifts,
and \code{"TC"} temporary changes are selected;
\code{"IO"} innovative outliers and \code{"SLS"} seasonal level shifts 
can also be selected.

The test statistics are based on the second equation defined in \code{\link{locate.outliers}}.

These functions are the called by \code{\link{locate.outliers}}.
The approach described in Chen & Liu (1993) is implemented to compute 
the \eqn{t}-statistics.

By default (\code{sigma = NULL}), the standard deviation of residuals 
is computed as the mean absolute deviation of \code{resid}. 
}

\value{
For each function, a two-column matrix is returned.
The first column contains the estimate of the coefficients related to the type of outlier 
and the second column contains the \eqn{t}-statistics.
The value of these statistics for each time point is stored by rows, thus 
the number of rows is equal to the length of \code{resid}.
}

\seealso{
\code{\link{locate.outliers}}, \code{\link{outliers.regressors}}.
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
\url{http://www.bde.es/f/webbde/SES/Secciones/Publicaciones/PublicacionesSeriadas/DocumentosTrabajo/99/Fic/dt9915e.pdf}
}

\examples{
# given an ARIMA model detect potential outliers
# for a critical value equal to 3.5
data("hicp")
y <- log(hicp[["011600"]])
fit <- arima(y, order = c(1, 1, 0), seasonal = list(order = c(2, 0, 2)))
resid <- residuals(fit)
pars <- coefs2poly(fit)
tstats <- outliers.tstatistics(pars, resid)
# potential observations affected by an additive outliers
which(abs(tstats[,"AO","tstat"]) > 3.5)
# potential observations affected by a temporary change
which(abs(tstats[,"TC","tstat"]) > 3.5)
# potential observations affected by a level shift
which(abs(tstats[,"LS","tstat"]) > 3.5)
}

\keyword{ts}
