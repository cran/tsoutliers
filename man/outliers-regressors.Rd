\name{outliers.regressors}
\alias{outliers.regressors}

\title{Regressor Variables for the Detection of Outliers}

\description{
These functions create regressor variables to be used included in the 
regression where tests for presence will be applied.
}

\usage{
outliers.regressors(pars, mo, n, weights = TRUE,
  delta = 0.7, freq = 12)
}

\arguments{
\item{pars}{a list containing the parameters of the model. 
See details section in \code{\link{locate.outliers}}.}
%\item{arcoefs}{a numeric vector. The coefficients of the product of the all the 
%autorregressive polynomials in the ARIMA model (non-seasonal, seasonal, 
%stationary, non-stationary).
%This argument is used only if the temporary change, \code{"TC"}, outlier is present in \code{"mo"}.}
%\item{macoefs}{a numeric vector. The coefficients of the product of the all the 
%moving average polynomials in the ARIMA model (non-seasonal, seasonal).
%This argument is used only if the temporary change \code{"TC"}, outlier is present in \code{"mo"}.}
\item{mo}{a data frame defining the type, location and weight of the outliers to be created.}
\item{n}{a numeric. The length of the variable that will contain the outlier.}
\item{weights}{logical. If \code{TRUE}, the variables are weighted by the 
values in column \code{"coefhat"} of the data frame \code{mo}. Otherwise, 
unit weights are considered.}
\item{delta}{a numeric. Parameter of the temporary change type of outlier.}
\item{freq}{a numeric. The periodicity of the data. 
Used only for the seasonal level shift, \code{"SLS"}.}
}

\details{
The variables returned by these functions are the regressors that take part in 
the second equation defined in \code{\link{locate.outliers}},
(equation (20) in Chen-Liu (1993), equation (3) in the documentat 
attached to the package).

Regressions are not actually run since the \eqn{t}-statistics
can be obtained more conveniently as indicated in equation (14) in Chen-Liu (1993).
These variables are used in function \code{\link{locate.outliers.iloop}} to  
adjust the residuals at each iteration.

The function \code{\link{outliers}} can be used to easily create the input 
argument \code{mo}.
}

\value{
A matrix containing the regressors by columms.
}

\seealso{
\code{\link{locate.outliers}}, \code{\link{outliers}}, 
\code{\link{outliers.tstatistics}}, \code{\link{tso}}.
}

\references{
Chen, C. and Liu, Lon-Mu (1993). 
\sQuote{Joint Estimation of Model Parameters and Outlier Effects in Time Series}.
\emph{Journal of the American Statistical Association},
\bold{88}(421), pp. 284-297.

Kaiser, R., and Maravall, A. (1999). 
\emph{Seasonal Outliers in Time Series}.
Banco de España, Servicio de Estudios. Working paper number 9915.
\url{http://www.bde.es/f/webbde/SES/Secciones/Publicaciones/PublicacionesSeriadas/DocumentosTrabajo/99/Fic/dt9915e.pdf}
}

\examples{
# regression of the residuals from the ARIMA model 
# on the corresponding regressors for three additive outliers
# at the 5% level, the first AO is not significant, the others are significant
data("hicp")
y <- log(hicp[["011600"]])
fit <- arima(y, order = c(1, 1, 0), seasonal = list(order = c(2, 0, 2)))
resid <- residuals(fit)
pars <- coefs2poly(fit)
mo <- outliers(rep("AO", 3), c(10, 79, 224))
xreg <- outliers.regressors(pars, mo, length(y))
summary(lm(residuals(fit) ~ 0 + xreg)) 
}

\keyword{ts}
