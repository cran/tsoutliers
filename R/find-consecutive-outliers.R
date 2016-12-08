
find.consecutive.outliers <- function(x, type)
{
  #NOTE this function is intended for internal use; 
  #no check is done for correctness of the arguments passed as input

  xnms0 <- rownames(x) # not used if 'nls==0' but mus be done here
  #apparently "%in%" is slightly faster than "=="; "subset" seems slower
  x <- x[x[,"type"] %in% type,] # subset(x, type == type)

  nls <- nrow(x)
  rmid <- NULL

  if (nls > 1)  
  {
    x <- x[order(x[,"ind"]),]
    tmp <- diff(x[,"ind"])
    if (any(tmp == 1))
    {
      ids <- c(0, which(tmp != 1), nls)
      #apparently here 'lapply' is slightly faster than 'sapply'
      aux <- lapply(as.list(seq_along(ids[-1])), 
        function(i, x) x[seq.int(ids[i]+1, ids[i+1])], x = seq_len(nls))
      #ignore cases not related to consecutive time points, (length(aux[[i]]) == 1)
      aux[which(lapply(aux, length) == 1)] <- NULL
      for (i in seq_along(aux))
      {
        #if (length(aux[[i]]) == 1)
        #  next 
        id <- aux[[i]][-which.max(abs(x[aux[[i]],"tstat"]))]
# #debug
#stopifnot(length(id) > 0)
        rmid <- c(rmid, id)
      }
      #rmid <- unlist(lapply(as.list(seq_along(ids[-1])), FUN=function(i, snls) {
      #  id <- snls[seq.int(ids[i]+1, ids[i+1])] 
      #  if (length(id)>1) id[-which.max(abs(x[id,"tstat"]))]}, snls=seq_len(nls)))
    }
    
    rmid <- pmatch(rownames(x)[rmid], xnms0)
  }

  #rownames(x)[rmid]
  rmid
}
