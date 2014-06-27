\name{JarqueBera.test}
\alias{JarqueBera.test}
\alias{print.mhtest}
%\alias{JarqueBera.test.default}
%\alias{JarqueBera.test.Arima}

\title{Jarque-Bera Test for Normality}

\description{
This function applies the test for normality proposed in Jarque and Bera (1980).
}

\usage{
JarqueBera.test(x, fc = 3.5, ...)
}

\arguments{
\item{x}{a time series of residuals or an object of class \code{Arima}.}
\item{fc}{a numeric. Factor to asses whether the first residual observations 
are to be omitted. Ignored if \code{x} is not an \code{Arima} object. See details.}
\item{...}{further arguments. Currently omitted.}
}

\details{
This function is based on function \code{\link[tseries]{jarque.bera.test}} 
available in package \pkg{tseries}.
Here, the results are split in a test for the null hypothesis that the 
skewness is \eqn{0}, the null that the kurtosis is \eqn{3} and the overall 
Jarque-Bera test.

The input can be a time series of residuals, \code{jarque.bera.test.default},
or an \code{Arima} object, \code{jarque.bera.test.Arima} from which the residuals 
are extracted.
In the former case the whole input series of residuals is used. 
In the latter case,
the first \eqn{n0} (defined below) residuals are omitted if they are are equal to zero
or if any of them are in absolute value larger than \code{fc} times 
the standard deviation of the remaining residuals.
\eqn{n0} is set equal to \code{x$arma[6] + x$arma[5] * x$arma[7]}, i.e. 
the number of regular differences times the periodicity of the data times 
the number of seasonal differences. If \eqn{n0} happens to be equal to \eqn{1} 
it is set to \eqn{2}.

If the latter trimming operation is not desired, 
the argument \code{fc} can be set to a high value to ensure the complete 
series of residuals in considered; or the function can be called 
as \code{jarque.bera.test(residuals(x))}.

Missing observations are omitted.
}

\value{
A list containing one \code{htest} object for the null hypothesis that 
the kurtosis is \eqn{3}, the skewness is \eqn{0} and a test combining 
both the kurtosis and the skewness to test for the normality of the input data.
}

\seealso{
\code{\link{print.mhtest}}.
}

\references{
Jarque, C. M. and Bera, A. K. (1980). 
\sQuote{Efficient test for normality, homoscedasticity and 
serial independence of residuals}.
\emph{Economic Letters}, \bold{6}(3), pp. 255-259.
}

\examples{
# fit an ARIMA model to the HICP 011600 series
# ARIMA(0,1,0)(2,0,1) was chosen by forecast::auto.arima(ic = "bic")
# normality of the residuals is rejected at the 5% significance level
# due to an excess of kurtosis
data("hicp")
y <- log(hicp[["011600"]])
fit1 <- arima(y, order = c(0, 1, 0), seasonal = list(order = c(2, 0, 1)))
JarqueBera.test(fit1)
JarqueBera.test(residuals(fit1))

# fit ARIMA model for the same series including outliers that were
# detected by "tso" and for the model chosen by "auto.arima"
# normality of the residuals is not rejected at the 5% significance level
# after including the presence of outliers
mo <- outliers(c("AO", "AO", "LS", "LS"), c(79, 210, 85, 225))
xreg <- outliers.effects(mo, length(y))
fit2 <- arima(y, order = c(1, 1, 0), seasonal = list(order = c(2, 0, 2)), 
  xreg = xreg)
JarqueBera.test(fit2)
}

\keyword{ts, htest}
