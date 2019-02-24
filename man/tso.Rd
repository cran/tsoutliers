\name{tso}
\alias{tso}
\alias{tso0}

\title{Automatic Procedure for Detection of Outliers}

\description{
These functions are the interface to the automatic detection 
procedure provided in this package.
}

\usage{
tso(y, xreg = NULL, cval = NULL, delta = 0.7, 
  types = c("AO", "LS", "TC"), 
  maxit = 1, maxit.iloop = 4, maxit.oloop = 4, cval.reduce = 0.14286, 
  discard.method = c("en-masse", "bottom-up"), discard.cval = NULL, 
  remove.method, remove.cval,
  tsmethod = c("auto.arima", "arima"), 
  args.tsmethod = NULL, logfile = NULL, check.rank = FALSE)

tso0(x, xreg = NULL, cval = 3.5, delta = 0.7, 
  types = c("AO", "LS", "TC"), maxit.iloop = 4, maxit.oloop = 4,
  discard.method = c("en-masse", "bottom-up"), discard.cval = NULL,   
  tsmethod = c("auto.arima", "arima"), args.tsmethod = NULL, 
  logfile = NULL, check.rank = FALSE)
}

\arguments{
\item{y}{a time series where outliers are to be detected.}
\item{x}{a time series object.}
\item{xreg}{an optional matrix of regressors with the same number of rows as \code{y}.}
\item{cval}{a numeric. The critical value to determine the significance of each type of outlier.}
\item{delta}{a numeric. Parameter of the temporary change type of outlier.}
\item{types}{a character vector indicating the type of outlier to be considered by the 
detection procedure: innovational outliers (\code{"IO"}), additive outliers (\code{"AO"}), 
level shifts (\code{"LS"}), temporary changes (\code{"TC"}) 
and seasonal level shifts (\code{"SLS"}).}
\item{maxit}{a numeric. The maximum number of iterations.}
\item{maxit.iloop}{a numeric. The maximum number of iterations in the inner loop. 
See \code{\link{locate.outliers}}.}
\item{maxit.oloop}{a numeric. The maximum number of iterations in the outer loop.}
\item{cval.reduce}{a numeric. Factor by which \code{cval} is reduced 
if the procedure is run on the adjusted series, if \code{maxit > 1}.}
\item{discard.method}{a character. The method used in the second stage of the procedure. 
See \code{\link{discard.outliers}}.}
\item{discard.cval}{a numeric. The critical value to determine the significance of each 
type of outlier in the second stage of the procedure (discard outliers). By default it is 
set equal to \code{cval}. See details.}
\item{remove.method}{deprecated, argument \code{discard.method} should be used.}
\item{remove.cval}{deprecated, argument \code{discard.cval} should be used.}
\item{tsmethod}{a character. The framework for time series modelling. It basically is the name 
of the function to which the arguments defined in \code{args.tsmethod} are referred to.}
\item{args.tsmethod}{an optional list containing arguments to be passed to 
the function invoking the method selected in \code{tsmethod}.}
%\item{args.tsmodel}{an optional list containing the arguments to be passed to 
%\code{\link[stsm]{stsm.model}}. Only for \code{tsmethod = "stsm"}}.
\item{logfile}{a character or \code{NULL}. It is the path to the file where 
tracking information is printed. Ignored if \code{NULL}.}
\item{check.rank}{logical. If \code{TRUE} the regressors are checked for perfect collinearity.
The variables related to coefficients that turn out to be \code{NA} due to possible 
perfect collinearity are discarded.}
}

\details{
Five types of outliers can be considered. 
By default: \code{"AO"} additive outliers, \code{"LS"} level shifts,
and \code{"TC"} temporary changes are selected;
\code{"IO"} innovative outliers and \code{"SLS"} seasonal level shifts 
can also be selected.

\code{tso0} is mostly a wrapper function around the functions 
\code{\link{locate.outliers}} and \code{\link{discard.outliers}}.

\code{tso} iterates around \code{tso0} first for the original series 
and then for the adjusted series. The process stops if no additional outliers 
are found in the current iteration or if \code{maxit} iterations are reached.

\code{tso0} is an auxiliar function (it is the workhorse for \code{tso} 
but it is not intended to be called directly by the user,
use \code{tso(maxit = 1, ...)} instead.
\code{tso0} does not check the arguments since they are assumed to be passed 
already checked by \code{tso}; the default value for \code{cval} is not based 
on the sample size.
For the time being, \code{tso0} is exported in the NAMESPACE since it 
is convenient for debugging.

If no value is specified for argument \code{cval} a default value based on the sample size 
is used. Let \eqn{n} be the number of observations. If \eqn{n \le 50} then \code{cval} is set 
equal to \eqn{3.0}; If \eqn{n \ge 450} then \code{cval} is set equal to \eqn{4.0};
otherwise \code{cval} is set equal to \eqn{3 + 0.0025 * (n - 50)}.

If \code{tsmethod} is \code{NULL}, the following default arguments are used in the 
function selected in \code{tsmethod}:
\code{tsmethod = "auto.arima"}: \code{allowdrift = FALSE}, \code{ic = "bic"};
\code{tsmethod = "arima"} = \code{order = c(0, 1, 1)} \code{seasonal = list(order = c(0, 1, 1))}.

If \code{args.tsmethod} is \code{NULL}, the following lists are used by default,
respectively for each method:
\code{auto.arima}: \code{list(allowdrift = FALSE, ic = "bic")};
\code{arima}: \code{list(order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1)))}.

\code{xreg} must be a matrix with time series attributes, \code{\link[stats]{tsp}}, 
that must be the same same as \code{tsp(x)}. Column names are also compulsory.
If there is only one regressor it may still have non-null dimension, i.e. 
it must be a one-column matrix.
  
The external regressors (if any) should be defined in the argument \code{xreg}.
However, they may be also defined as an element in \code{args.tsmethod} since this list 
is passed to function that fits the model. The function \code{tso} 
deals with this possibility and returns a warning if \code{"xreg"} is defined 
twice with different values. No checks are done in \code{tso0}.
  
If \code{maxit = 1} the procedure is run only once on the original series.
If \code{maxit > 1} the procedure is run iteratively, first for the original series 
and then for the adjusted series. The critical value used for the adjusted series 
may be reduced by the factor \code{cval.reduce}, equal to \eqn{0.14286} by default.
The new critical value is defined as \eqn{cval * (1 - cval.reduce)}.

By default, the same critical value is used in the first stage of the procedure 
(location of outliers) and in the second stage (discard outliers). 
Under the framework of structural time series models I noticed that 
the default critical value based on the sample size is too high, since all the 
potential outliers located in the first stage were discarded in the second stage 
(even in simulated series with known location of outliers). 
In order to investigate this issue, the argument \code{discard.cval} has been added. 
In this way a different critical value can be used in the second stage. 
Alternatively, the argument \code{discard.cval} could be omitted and simply choose 
a lower critical value, \code{cval}, to be used in both stages.
However, using the argument \code{discard.cval} is more convenient since it avoids 
locating too many outliers in the first stage.
%If a lower critical value is chosen for the whole procedure, \code{cval}, 
%then potential set of outliers located in the first stage is much larger, although 
%most of them are discarded in the second stage. 
\code{discard.cval} is not affected by \code{cval.reduce}.
}

\value{
A list of class \code{tsoutliers}.
}

\seealso{
\code{\link{locate.outliers}}, \code{\link{discard.outliers}}, 
\code{\link{plot.tsoutliers}}, \code{\link{print.tsoutliers}}.
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
}

\examples{
\dontrun{
data("hicp")
tso(y = log(hicp[[1]]))
}
}

\keyword{ts}
