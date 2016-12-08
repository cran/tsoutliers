\name{tsoutliers-package}
\docType{package}
\alias{tsouliers-package}

\title{Automatic Detection of Outliers in Time Series}

\description{
This package implements a procedure based on the approach described in 
Chen and Liu (1993) for automatic detection of outliers in time series. 
Innovational outliers, additive outliers, level shifts, temporary changes 
and seasonal level shifts are considered.
}

\details{
Time series data often undergo sudden changes that alter
the dynamics of the data transitory or permanently.
These changes are typically non systematic and cannot be captured by 
standard time series models. 
That's why they are known as exogenous or outlier effects. 
Detecting outliers is important because they have an impact on 
the selection of the model, the estimation of parameters and, 
consequently, on forecasts.

Following the approach described in Chen & Liu (1993), 
an automatic procedure for detection of outliers in time series 
is implemented in the package \pkg{tsoutliers}.
The procedure may in turn be run along with the automatic ARIMA model 
selection strategy available in the package \pkg{forecast}.

The function \code{\link{tso}} is the main interface for the 
automatic procedure. The functions \code{\link{locate.outliers.oloop}}
and \code{\link{remove.outliers}} implement respectively the first and 
second stages of the procedure. In practice, the user may stick to use 
the function \code{\link{tso}}. 

Although the purpose of the package is to provide an automatic procedure, 
the implementation allows the user to do a manual inspection of each 
step of the procedure. Thus, the package is also useful to track the behaviour
of the procedure and come up with ideas for possible improvements or 
enhancements.

%For ease of maintenance and debugging the procedure is split in several 
%functions specialised in a stage or substage of the procedure. 
The functions \code{\link{locate.outliers.oloop}} and \code{\link{remove.outliers}} 
implement the major steps of the procedure.
\code{\link{tso}} is the main interface to the automatic procedure. 
All the options at any stage of the procedure can be defined through the 
arguments passed to \code{\link{tso}}.
Despite the user may stick to use the function \code{\link{tso}},
other functions called by this main interface are exported in the namespace of 
the package. They are helpful for debugging and allow the interested user to more 
easily track each step of the procedure.

Information supplemental to these help pages is given in the document that is 
provided with the package (\file{tsoutliers/inst/doc/tsoutliers.pdf} in the source files).
}

\author{
Javier López-de-Lacalle \email{javlacalle@yahoo.es}
}

\references{
Chen, C. and Liu, Lon-Mu (1993). 
\sQuote{Joint Estimation of Model Parameters and Outlier Effects in Time Series}.
\emph{Journal of the American Statistical Association},
\bold{88}(421), pp. 284-297.
\doi{10.2307/2290724}

Gómez, V. and Maravall, A. (1996). 
\emph{Programs TRAMO and SEATS. Instructions for the user}.
Banco de España, Servicio de Estudios. Working paper number 9628.
\url{http://www.bde.es/f/webbde/SES/Secciones/Publicaciones/PublicacionesSeriadas/DocumentosTrabajo/96/Fich/dt9628e.pdf}

Gómez, V. and Taguas, D. (1995). 
\emph{Detección y Corrección Automática de Outliers con TRAMO: 
Una Aplicación al IPC de Bienes Industriales no Energéticos}.
Ministerio de Economía y Hacienda. Document number D-95006.
\url{http://www.sepg.pap.minhap.gob.es/sitios/sepg/es-ES/Presupuestos/Documentacion/Documents/DOCUMENTOS\%20DE\%20TRABAJO/D95006.pdf}

Hyndman, R.J. and Khandakar, Y. (2008). 
\sQuote{Automatic Time Series Forecasting: The \pkg{forecast} Package for \R}.
\emph{Journal of Statistical Software}, \bold{27}(3), pp. 1-22.
\url{http://www.jstatsoft.org/v27/i03}

Hyndman, R.J. with contributions from 
George Athanasopoulos, Slava Razbash, Drew Schmidt, Zhenyu Zhou, Yousaf Khan, 
Christoph Bergmeir and Earo Wang (2014). 
\sQuote{\pkg{forecast}: Forecasting functions for time series and linear models}.
\R package version 5.4. 
\url{https://CRAN.R-project.org/package=forecast}
}

\keyword{package, ts}
