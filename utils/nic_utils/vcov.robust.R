#' Clustered statistics from glm with link function of logit
#' 
#' @param mdl a model object from glm package for logistic regression.
#' @param cluster a column vector for cluster id
#' @param type type of model lr = logistic regression, lm linear model
#' @returns vcov a numeric variance covariance matrix




vcov.robust <- function(mdl, cluster, family=c("binomial","gaussian")[1]){
  # initiate variance covariance matrix to return
  vcov <- NULL
  
  # get x, y, b from the model object
  if(any(class(mdl)%in%c("glm", "lm") ) ){
    # find b
    b <- mdl$coefficients
    # find x
    mdl$data$`(Intercept)` <- 1
    x <- mdl$data[,names(b)]
    # find y 
    y <- mdl$y
    # prepare c
    c <- cluster
    if(length(c)!=length(y)){
      stop("cluster length and number of rows in modeling data differ!")
    }
  }
  
  if(family=="binomial"){
    vcov <- vcov.robust.xy(x,y,b,c)$cHScov
  }else if(family=="gaussian"){
    vcov <- vcov.robust.xy_lm(x,y,b,c)$cHScov
  }
  
  
  return(list(vcov = vcov))
}
