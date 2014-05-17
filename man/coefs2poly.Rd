\name{coefs2poly}
\alias{coefs2poly}
%\alias{coefs2poly.default}
%\alias{coefs2poly.Arima}

\title{Product of the Polynomials in an ARIMA Model}

\description{
This function collapses the polynomials of an ARIMA model 
into two polynomials: the product of the autoregressive polynomials 
and the product of the moving average polynomials.
}

\usage{
coefs2poly(x, morder, add = TRUE, ...)
}

\arguments{
\item{x}{a numeric vector containing the coefficients from the 
fitted ARIMA model or an \code{Arima} object returned by \code{\link[stats]{arima}}.}
%
\item{morder}{a numeric vector. The order of the fitted model as returned 
by \code{\link[stats]{arima}} in element \code{arma}. Ignored if \code{x} is an 
\code{Arima} object.}
%
\item{add}{logical. If \code{TRUE}, the polynomial of the differencing filter 
(if present in the model) is multiplied byt the stationary autoregressive 
polynomial. Otherwise only the coefficients of the product of the 
stationary polynomials is returned.}
%
\item{...}{Further arguments to be passed to other functions. Currently ignored.}
}

\details{
In practice, the version \code{coefs2poly.Arima} may be more convenient since it 
requires passing only one argument defining the model. However, since only the 
coefficients of the model and the order of the model is required by this function, 
there is no need to pass the complete \code{Arima} object.
}

\value{
A list containing the elements:
\code{arcoefs}, the coefficients of the product of the 
autoregressive polynomials;
\code{macoefs}, the coefficients of the product of the 
moving average polynomials. 
This list is of class \code{"ArimaPars"} so that it can be recognized by 
\code{\link{outliers.tstatistics}}.
}

\note{
If the ARIMA model contains regressor variables, 
the names of those variables cannot match the following 
regular expressions:
\dQuote{^ar\\d+$}, \dQuote{^sar\\d+$}, \dQuote{^ma\\d+$}, 
and \dQuote{^sma\\d+$}. Otherwise the AR and MA coefficients 
would be misled by the coefficients of the external regressors.
}

\seealso{
\code{\link[polynom]{polynomial}}, 
\code{\link[polynom]{Ops.polynomial}}.
}

\examples{
# ARIMA(0,1,1)(0,1,1) model
fit <- arima(log(AirPassengers), order = c(0,1,1), 
  seasonal = list(order = c(0,1,1)))
coefs <- coef(fit)

# "coefs2poly" returns the coefficients of the product of 
# the non-seasonal and the seasonal moving average polynomials
pma <- polynom::polynomial(c(1, coefs[1]))
psma <- polynom::polynomial(c(1, rep(0, 11), coefs[2]))
coef(pma * psma)[-1]
coefs2poly(coef(fit), fit$arma)$macoefs

# since the model does not contain an autoregressive part
# the product of the regular and the seasonal differencing 
# filter is returned if "add = TRUE"
coefs2poly(coef(fit), fit$arma)$arcoefs
# an empty set nothing is returned if "add = FALSE"
coefs2poly(coef(fit), fit$arma, add = FALSE)$arcoefs

# in a model with non-seasonal part and no differencing filter 
# no multiplication of polynomials are involved and 
# the output coincides with "coef"
fit <- arima(log(AirPassengers), order = c(1,0,1))
coef(fit)
coefs2poly(coef(fit), fit$arma)
}

\keyword{ts, math, symbolmath}
