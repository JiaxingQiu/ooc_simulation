#' Clustered statistics from x, y, coefficient and cluster matrices
#' 
#' @param x matrix of predictors, of dimension \eqn{n \times p}{n*p}; each row
#' is an observation vector.
#' @param y response vector of length n. This argument should be quantitative for
#' regression (least squares), and a two-level factor for classification
#' (logistic model, huberized SVM, squared SVM).
#' @param b estimated coefficient of length p.
#' @param c cluster ids vector of length n.
#' @returns vcov a numeric variance covariance matrix


vcov.robust.xy <- function(x,y,b,c){
  
  
  # setup, error checking
  # check dim
  x <- as.matrix(x) # "x must be a n*p matrix"
  nobs <- as.integer(dim(x)[1]) # number of observations
  nvar <- as.integer(dim(x)[2]) # number of coefficients
  
  if(length(y) != nobs) stop("x and y have different number of rows")
  y <- matrix(y, nrow = nobs)
  if(dim(y)[1]!=nobs|dim(y)[2]!=1) stop("y must be a n*1 matrix")
  if(!all(as.numeric(y)%in%c(0,1))) stop("y can only take value 0 or 1")# check value
  
  if(length(c) != nobs) stop("x and c have different number of rows")
  c <- matrix(c, nrow = nobs)
  if(dim(c)[1]!=nobs|dim(c)[2]!=1) stop("c must be a n*1 matrix")
  
  if(length(b) != nvar) stop("x and b have different number of columns")
  b <- matrix(b, nrow = nvar)
  if(dim(b)[1]!=nvar|dim(b)[2]!=1) stop("b must be a p*1 matrix")
  
  
  # predicted probability p [n*1]
  p <- exp(x%*%b)/(1+exp(x%*%b))
  
  # calculate the gradient g [1*p]
  g <- t(y-p) %*% x
  
  # calculate information matrix J [p*p]
  J <- t(x * as.numeric(p*(1-p))) %*% x
  # J <- t(x) %*% diag(as.numeric(p*(1-p))) %*% x # very slow
  
  # Deviance = -2 * log likelihood
  deviance <- -2 * sum(log(p)*y + log(1-p)*(1-y),na.rm = T)
  
  # -------------------- unclustered ------------------------
  # the derivative of log likelihood per component gi [n*p]
  gi <- x
  for(i in 1:dim(x)[2]) gi[,i] <- (y-p)*x[,i]
  
  # calculate fisher information matrix K [p*p]
  K <- t(gi)%*%gi
  
  # Huber Sandwich Estimator matrix without clustering
  HScov <- solve(J)%*%K%*%solve(J)
  
  # approximation of the number of parameters (effective d.f. of the model (counting intercept terms))
  dof <- length(diag(K%*%solve(J))) # sum(diag(K%*%solve(J)))
  
  # AIC
  aic <- deviance + 2*dof
  
  # NIC(AIC)
  nic <- deviance + 2*sum(diag(K%*%solve(J)))
  
  
  # ---------------------- Clustered ------------------------
  # the derivative of log likelihood per cluster cgi [n*p]
  cgi <- NULL
  for(cid in unique(c)){
    if(sum(c==cid)>1){
      gi_c <- colSums(gi[c==cid,]) # gi cluster of id
    }else if(sum(c==cid)==1){
      gi_c <- gi[c==cid,]
    }
    if(is.null(cgi)){
      cgi <- matrix(gi_c, ncol=nvar) 
    }else{
      cgi <- rbind(cgi, matrix(gi_c,ncol=nvar))
    }
  }
  
  # calculate clustered fisher information matrix K [p*p]
  cK <- t(cgi)%*%cgi
  
  # Huber Sandwich Estimator matrix with clustering
  cHScov <- solve(J)%*%cK%*%solve(J)
  
  # estimated number of parameters
  cdof <- sum(diag(cK%*%solve(J)))
  
  # NIC
  nicc <- deviance + 2*cdof
  
  # calculate NIC
  return(list(deviance = deviance,
              K = K,
              cK = cK,
              J = J,
              HScov = HScov,
              dof = dof,
              cHScov = cHScov,
              cdof = cdof,
              aic = aic,
              nic = nic,
              nicc = nicc))
  
}
