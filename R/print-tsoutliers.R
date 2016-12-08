
print.tsoutliers <- function(x, digits = max(3L, getOption("digits") - 3L), call = FALSE, ...)
{
  #avoid long output returned when 'arima' is called from 'do.call()'
  #element 'method' is required below by 'print', so it must be kept
  if (!call)
    x$fit$call <- list(method=x$fit$call$method)

  print(x$fit)

  if (nrow(x$outliers) > 0)
  {
    cat("\nOutliers:\n")
    print(format(x$outliers, digits = digits)) #...
  } else
    cat("\nNo outliers were detected.\n")
  
  invisible(x)
}
