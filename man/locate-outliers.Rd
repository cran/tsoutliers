\name{locate.outliers}
\alias{locate.outliers}

\title{Stage I of the Procedure: Locate Outliers (Baseline Function)}

\description{
This function applies the \eqn{t}-statistics for the significance of outliers 
at every time point and selects those that are significant given a critical value.
}

\usage{
locate.outliers(resid, pars, cval = 3.5, types = c("AO", "LS", "TC"), delta = 0.7)
}

\arguments{
\item{resid}{a time series. Residuals from a time series model fitted to the data.}
%
\item{pars}{a list containing the parameters of the model fitted to the data. 
See details below.}
%
\item{cval}{a numeric. The critical value to determine the significance of each 
type of outlier.}
%
\item{types}{a character vector indicating the types of outliers to be considered.}
%
\item{delta}{a numeric. Parameter of the temporary change type of outlier.}
}

\details{
Five types of outliers can be considered. 
By default: \code{"AO"} additive outliers, \code{"LS"} level shifts,
and \code{"TC"} temporary changes are selected;
\code{"IO"} innovative outliers and \code{"SLS"} seasonal level shifts 
can also be selected.

The approach described in Chen & Liu (1993) is followed to locate outliers.
The original framework is based on ARIMA time series models. 
The extension to structural time series models is currently experimental.

Let us define an ARIMA model for the series \eqn{y_t^*} subject 
to \eqn{m} outliers defined as \eqn{L_j(B)} with weights \eqn{w}:

\deqn{y_t^* = \sum_{j=1}^m \omega_j L_j(B) I_t(t_j) + 
\frac{\theta(B)}{\phi(B) \alpha(B)} a_t \,,}{%
y_t^* = \sum_{j=1}^m \omega_j L_j(B) I_t(t_j) + 
\frac{\theta(B)}{\phi(B) \alpha(B)} a_t,
}

where \eqn{I_t(t_j)} is an indicator variable containing the value 
\eqn{1} at observation \eqn{t_j} where the \eqn{j}-th outlier arises;
\eqn{\phi(B)} is an autoregressive polynomial with all roots outside the unit circle;
\eqn{\theta(B)} is a moving average polynomial with all roots outside the unit circle;
and \eqn{\alpha(B)} is an autoregressive polynomial with all roots on the unit circle.

The presence of outliers is tested by means of \eqn{t}-statistics 
applied on the following regression equation:

\deqn{\pi(B) y_t^* \equiv \hat{e}_t = 
\sum_{j=1}^m \omega_j \pi(B) L_j(B) I_t(t_j) + a_t \,.}{%
\pi(B) y_t^* = \hat{e}_t = 
\sum_{j=1}^m \omega_j \pi(B) L_j(B) I_t(t_j) + a_t.}

where \eqn{\pi(B) = \sum_{i=0}^\infty \pi_i B^i}.
%
The regressors of the above equation are created by the functions
\code{outliers.regressors.arima} and the remaining functions described here.

The function \code{locate.outliers} computes all the \eqn{t}-statistics for each type of 
outlier and for every time point. See \code{\link{outliers.tstatistics}}.
Then, the cases where the corresponding \eqn{t}-statistic are (in absolute value) 
below the threshold \code{cval} are removed. Thus, a potential set of outliers is obtained.

Some polishing rules are applied by \code{locate.outliers}:
%
\itemize{
\item If level shifts are found at consecutive time points, only then point with higher 
\eqn{t}-statistic in absolute value is kept.
%
\item If more than one type of outlier exceed the threshold \code{cval} at a given time point, 
the type of outlier with higher \eqn{t}-statistic in absolute value is kept and the others 
are removed.
%
%\item If at a given time point an innovative outlier (IO) and other types of outlier are 
%detected, the IO is removed and the previous rule is applied for the remaining outliers. 
}

The argument \code{pars} is a list containing the parameters of the model.
In the framework of ARIMA models, the coefficients of the ARIMA must be defined in \code{pars} 
as the product of the autoregressive non-seasonal and seasonal polynomials (if any) and the 
differencing filter (if any). The function \code{\link{coefs2poly}} can be used to define 
the argument \code{pars}. 
}

\value{
A data frame defining by rows the potential set of outliers. 
The type of outlier, the observation, the coefficient and the \eqn{t}-statistic
are given by columns respectively for each outlier.
}

\note{
The default critical value, \code{cval}, is set equal to \eqn{3.5} and, 
hence, it is not based on the sample size as in functions \code{\link{tso}} 
or \code{\link{locate.outliers.oloop}}.
%this function is mainly an internal functions but it is exported 
%for debugging and to allow the interested user to more easily track 
%each step of the procedure.
%Same for \code{\link{locate.outliers.oloop}}.

Currently the innovational outlier \code{"SLS"} is not available 
if \code{pars} is related to a structural time series model.
}

\seealso{
\code{\link{locate.outliers.oloop}},
\code{\link{locate.outliers.iloop}}, 
\code{\link{outliers.tstatistics}}, 
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
data("hicp")
y <- log(hicp[["011600"]])
fit <- arima(y, order = c(1, 1, 0), seasonal = list(order = c(2, 0, 2)))
resid <- residuals(fit)
pars <- coefs2poly(fit)
outliers <- locate.outliers(resid, pars)
outliers
}

\keyword{ts}
