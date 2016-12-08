
coefs2poly <- function(x, add = TRUE)
{
  arcoefs <- x$model$phi

  if (add && x$arma[6] > 0)
  {
    #multiply by (1-L) x$arma[6] times
    #p <- 1
    #for (i in seq_len(x$arma[6])) 
    #  p <- c(0, p) - c(p, 0)
    #arcoefs <- -convolve(c(1, -arcoefs), p, type="open")[-1]
    arcoefs <- c(1, -arcoefs)
    for (i in seq_len(x$arma[6]))
      arcoefs <- c(arcoefs[1], diff(arcoefs), -arcoefs[length(arcoefs)])
    arcoefs <- -arcoefs[-1]
  }

  if (add && x$arma[7] > 0)
  {
    # multiply by (1-L^S) x$arma[7] times (at most 2 times)
    #
    #arcoefs <- -convolve(c(1, -arcoefs), c(-1, rep(0, x$arma[5]-1), 1), type="open")[-1]
    #if (x$arma[7] == 2) {
    #  arcoefs <- -convolve(c(1, -arcoefs), c(-1, rep(0, x$arma[5]-1), 1), type="open")[-1]
    #} else 
    #if (x$arma[7] > 2)
    #  stop("unsupported model, ", sQuote("x$arma[7] > 2"))
    arcoefs <- c(1, -arcoefs)
    arcoefs <- c(arcoefs[1], diff(c(rep(0, x$arma[5]-1), arcoefs, rep(0, x$arma[5]-1)), x$arma[5]), 
      -arcoefs[length(arcoefs)])
    if (x$arma[7] == 2) {
      arcoefs <- c(arcoefs[1], diff(c(rep(0, x$arma[5]-1), arcoefs, rep(0, x$arma[5]-1)), x$arma[5]), 
        -arcoefs[length(arcoefs)])
    } else 
    if (x$arma[7] > 2)
      stop("unsupported model, ", sQuote("x$arma[7] > 2"))
    arcoefs <- -arcoefs[-1]
  }

  #if (any(x$arma[c(6,7)] > 0))
  #  arcoefs[abs(arcoefs) < .Machine$double.eps] <- 0

  structure(list(arcoefs=arcoefs, 
    macoefs=x$model$theta[seq_len(x$arma[2] + x$arma[5]*x$arma[4])]), 
    class = "ArimaPars")
}
