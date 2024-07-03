# robust covariance matrix for linear model

vcov.robust.xy_lm <- function(x,y,b,c){
  
  # setup, error checking
  # check dim
  x <- as.matrix(x) # "x must be a n*p matrix"
  nobs <- as.integer(dim(x)[1]) # number of observations
  nvar <- as.integer(dim(x)[2]) # number of coefficients
  
  if(length(y) != nobs) stop("x and y have different number of rows")
  y <- matrix(y, nrow = nobs)
  if(dim(y)[1]!=nobs|dim(y)[2]!=1) stop("y must be a n*1 matrix")
  
  if(length(c) != nobs) stop("x and c have different number of rows")
  c <- matrix(c, nrow = nobs)
  if(dim(c)[1]!=nobs|dim(c)[2]!=1) stop("c must be a n*1 matrix")
  
  if(length(b) != nvar) stop("x and b have different number of columns")
  b <- matrix(b, nrow = nvar)
  if(dim(b)[1]!=nvar|dim(b)[2]!=1) stop("b must be a p*1 matrix")
  
  
  # predicted response
  p <- x%*%b
  
  # sigma estimation
  s <- sqrt(1/nobs * sum( (y-p)^2 )) # check with sigma(mdl)
  
  # calculate information matrix J [p*p] (- Hessian matrix)
  J <- (1/s^2) * t(x) %*% x
  
  # MLE Deviance = -2 * log likelihood
  deviance <- nobs*log(2*pi*s^2) + 1/s^2*sum( (y-p)^2 )  # check with AIC(mdl)

  # -------------------- unclustered ------------------------
  # the gradient of log likelihood per observation per variable i
  gi <- x
  for(i in 1:dim(x)[2]) gi[,i] <- 1/s^2 * (y - p)*x[,i]
  
  # calculate fisher information matrix K [p*p]
  K <- t(gi)%*%gi
  
  # Huber Sandwich Estimator matrix without clustering
  HScov <- solve(J)%*%K%*%solve(J) # check with sandwich::vcovHC(mdl, type = "HC")
  # vcov(mdl) check with sandwich::vcovHC(mdl, type = "const")
  
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
  cHScov <- solve(J)%*%cK%*%solve(J)  # check with sandwich::vcovCL(mdl, cluster = mdl$c)
  
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
