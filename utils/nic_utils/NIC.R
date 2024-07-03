NIC <- function(mdl, family=c("binomial","gaussian")[1]){
  tryCatch({
    b <- coef(mdl)[!is.na(coef(mdl))]
    x <- as.matrix(model.matrix(mdl))[,names(b)]
    y <- mdl$y
    c <- mdl$c
    if(family == "binomial") res <- vcov.robust.xy(x, y, b, c)
    if(family == "gaussian") res <- vcov.robust.xy_lm(x, y, b, c)
    
    return(list("dev"=res$deviance,
                "aic"=res$aic,
                "nic" = res$nic,
                "nicc"=res$nicc,
                "robcov" = res$cHScov,
                "cov" = res$HScov ))
    
  },error=function(e){
    print(e)
  })
}
