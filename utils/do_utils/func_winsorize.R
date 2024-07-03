winsorize <- function(x){
  if (!is.null(ncol(x))){ # if x is a dataframe with more than 1 column
    # winzorize all the numeric columns in a dataframe object
    for (i in which(sapply(x, is.numeric))) {
      quantiles <- quantile( x[,i], c(.001, .999 ), na.rm =TRUE)
      x[,i] = ifelse(x[,i] < quantiles[1] , quantiles[1], x[,i])
      x[,i] = ifelse(x[,i] > quantiles[2] , quantiles[2], x[,i])
    }
  }else{ # if x is vector or only one column of numeric variable in dataframe 
    quantiles <- quantile(x, c(.001, .999 ), na.rm =TRUE)
    x = ifelse(x < quantiles[1] , quantiles[1], x)
    x = ifelse(x > quantiles[2] , quantiles[2], x)
  }
  return(x)
}
