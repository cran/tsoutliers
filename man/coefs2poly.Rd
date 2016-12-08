\name{coefs2poly}
\alias{coefs2poly}

\title{Product of the Polynomials in an ARIMA Model}

\description{
This function collapses the polynomials of an ARIMA model 
into two polynomials: the product of the autoregressive polynomials 
and the product of the moving average polynomials.
}

\usage{
coefs2poly(x, add = TRUE)
}

\arguments{
\item{x}{an object of class \code{Arima}, as returned by \code{\link[stats]{arima}}.}
%
\item{add}{logical. If \code{TRUE}, the polynomial of the differencing filter 
(if present in the model) is multiplied by the stationary autoregressive 
polynomial. Otherwise only the coefficients of the product of the 
stationary polynomials is returned.}
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

\examples{
# ARIMA(0,1,1)(0,1,1) model
fit <- arima(log(AirPassengers), order = c(0,1,1), 
  seasonal = list(order = c(0,1,1)))
coefs <- coef(fit)

# "coefs2poly" returns the coefficients of the product of 
# the non-seasonal and the seasonal moving average polynomials
a1 <- convolve(c(1, coefs[1]), rev(c(1, rep(0, 11), coefs[2])), type="open")[-1]
a2 <- coefs2poly(fit)$macoefs
a2
all.equal(a1, a2, check.names=FALSE)

# since the model does not contain an autoregressive part
# the product of the regular and the seasonal differencing 
# filter is returned if "add = TRUE"
coefs2poly(fit)$arcoefs
# an empty set is returned if "add = FALSE"
coefs2poly(fit, add = FALSE)$arcoefs

# in a model with non-seasonal part and no differencing filter 
# no multiplication of polynomials are involved and 
# the coefficients are the same as those returned by "coef"
fit <- arima(log(AirPassengers), order = c(1,0,1))
coef(fit)
coefs2poly(fit)
}

\keyword{ts, math, symbolmath}
