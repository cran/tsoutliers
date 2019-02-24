\name{discard.outliers}
\alias{discard.outliers}

\title{Stage II of the Procedure: Discard Outliers}

\description{
This functions tests for the significance of a given set of outliers 
in a time series model that is fitted including the outliers as regressor variables.
}

\usage{
discard.outliers(x, y, cval = NULL, 
  method = c("en-masse", "bottom-up"), 
  delta = 0.7, tsmethod.call = NULL, 
  fdiff = NULL, logfile = NULL, check.rank = FALSE)
}

\arguments{
\item{x}{a list. The output returned by \code{\link{locate.outliers.oloop}}.}
\item{y}{a time series.}
\item{cval}{a numeric. The critical value to determine the significance of each type of outlier.}
\item{method}{a character. The method to discard/remove outliers. See details.}
\item{delta}{a numeric. Parameter of the temporary change type of outlier.}
%\item{n.start}{a numeric. The number of warming observations added to the 
%input passed to the Kalman filter. Only for \code{tsmethod = "stsm"}.}
\item{tsmethod.call}{an optional \code{call} object. The call to the function used to fit the time 
series model.}
%\item{tsmethod}{a character. The name of the function used to fit the time series model.
%Omitted if \code{tsmethod.call} is not \code{NULL}.}
%\item{args.tsmethod}{a list. Arguments to be passed to the function chosen to fit the time series model. 
%Omitted if \code{tsmethod.call} is not \code{NULL}.}
%\item{fdiff}{an optional function. Differencing filter to be applied on the original data. 
%Used if \code{method="linear-regression"}.}
\item{fdiff}{currently ignored.}
\item{logfile}{a character or \code{NULL}. It is the path to the file where 
tracking information is printed. Ignored if \code{NULL}.}
\item{check.rank}{logical. If \code{TRUE}, variables generating perfect multicollinearity 
are removed (tentative implementation).}
}

\details{
%ARIMA model is chosen again if tsmethod is "auto.arima"

In the regressions involved in this function, the variables included as regressors 
stand for the effects of the outliers on the data. 
These variables are the output returned by \code{\link{outliers.effects}} 
not by \code{\link{outliers.regressors}}, which returns the regressors used in the 
auxiliar regression where outliers are located 
(see second equation defined in \code{\link{locate.outliers}}).

The outliers are defined in input \code{x}. If there are regressor variables 
in \code{tsmethod.call$xreg} they are considered as other regressor variables 
that are included in the regression to test for the significance of outliers.

Given a set of potential outliers detected by \code{\link{locate.outliers}} and 
\code{\link{locate.outliers.oloop}}, three methods are considered in order to 
determine which outliers are not significant after refitting the model 
(including all the potential outliers):
%
\itemize{
\item \code{"en-masse"}: The complete set of outliers is included as regressor variables and the 
model is fitted again. Those outliers that turn out to be not significant for the critical 
value \code{cval} are discarded/removed. The procedure is iterated until all the outliers are significant 
in the final set of outliers.
%
\item \code{"bottom-up"}: First the, the outlier with larger \eqn{t}-statistic is included in the 
model. If it is significant the presence of the outlier is confirmed. Otherwise it is discarded. 
Then, the next outlier with larger \eqn{t}-statistic is included along with the first outlier and 
tested for significance. 
If after including a new outlier, e.g. the \eqn{i}-th outlier, the outliers that have been 
confirmed so far significant become not significant, then the \eqn{i}-th outlier is discarded 
regardless of the value of its \eqn{t}-statistic.
%
%\item \code{"linear-regression"}: The original series and the outlier regressors are 
%differenced by means of the function \code{fdiff}. Then a linear regression of the the 
%differenced data on the differenced regressors is performed. Those outliers that are not significant 
%at the 5\% significant value are discarded. The p-value from the linear regression is used, 
%not the critical value \code{cval}, which is omitted with this method.
}

The option \code{"en-masse"} may be preferred to \code{"bottom-up"} when there are are several outliers,
since it may be hard to fit an ARIMA model with many regressor variables.

%The method \code{"linear-regression"} is an experimental version intended to be used 
%with \code{tsmethod = "stsm"}, where fitting the model with external regressors seems to be
%harder. In this method the critical value \code{cval} is not used, 
%the p-value of the t-statistics is used instead for the \eqn{5\%} significance level.
}

\value{
A list containing the following elements: 
\code{xreg}, the variables used as regressors;
\code{xregcoefs}, the coefficients of the outlier regressors;
\code{xregtstats}, the \eqn{t}-statistics of the outlier regressors;;
\code{iter}, the number of iterations used by method \code{"en-masse"};
\code{fit}, the fitted model;
\code{outliers}, the set of outliers after removing those that were not significant.
}

\seealso{
\code{\link{locate.outliers}}, \code{\link{tso}}.
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
}

\examples{
\dontrun{
data("hicp")
y <- log(hicp[["011600"]])
fit <- arima(y, order = c(1, 1, 0), seasonal = list(order = c(2, 0, 2)))
# initial set of outliers
res <- locate.outliers.oloop(y, fit, types = c("AO", "LS", "TC"))
res$outliers
# given the model fitted above, the effect on the data of some of 
# the outliers is not significant (method = "en-masse")
discard.outliers(res, y, method = "en-masse", 
  tsmethod.call = fit$call)$outliers
# in this case, using method = "bottom-up" the first four 
# outliers with higher t-statistic are kept
discard.outliers(res, y, method = "bottom-up", 
  tsmethod.call = fit$call)$outliers
}
}

\keyword{ts}
