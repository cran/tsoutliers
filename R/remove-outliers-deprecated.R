
remove.outliers <- function(x, y, cval = NULL, 
  method = c("en-masse", "bottom-up"), 
  delta = 0.7, #n.start = 50, 
  tsmethod.call = NULL, fdiff = NULL, logfile = NULL)
{
  # as of version 0.6-6, 'remove.outliers' has been renamed as 'discard.outliers'
  .Deprecated("discard.outliers")

  discard.outliers(x, y, cval, method, delta, #n.start, 
    tsmethod.call, fdiff, logfile)
}
