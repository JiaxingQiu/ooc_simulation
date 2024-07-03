gglasso_cv_logit <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
  typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss", auc = "AUC")
  if (pred.loss == "default") 
    pred.loss <- "misclass"
  if (!match(pred.loss, c("misclass", "loss", "auc"), FALSE)) {
    warning("Only 'misclass', 'loss', 'auc' available for logistic regression; 'misclass' used")
    pred.loss <- "misclass"
  }
  prob_min <- 1e-05
  fmax <- log(1/prob_min - 1)
  fmin <- -fmax
  ###Turn y into c(0,1)
  y <- as.factor(y)
  y <- c(-1, 1)[as.numeric(y)]
  nfolds <- max(foldid)
  predmat <- matrix(NA, length(y), length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    which <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
    nlami <- length(outlist[[i]]$lambda)
    predmat[which, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  predmat <- pmin(pmax(predmat, fmin), fmax)
  if(pred.loss %in%c("loss", "misclass")){
    cvraw <- switch(pred.loss, 
                    loss = 2 * log(1 + exp(-y * predmat)), 
                    misclass = (y != ifelse(predmat > 0, 1, -1)))
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))
  }
  if(pred.loss == "auc"){
    C_fun <- function(y_prob, y_true){
      y_true <- ifelse(y_true>0,1,0)
      AUROC <- ifelse(is.na(MLmetrics::AUC(y_prob, y_true)), round(pROC::auc(pROC::roc(y_true, round(y_prob,6))),6), MLmetrics::AUC(y_prob, y_true))
      AUROC[which(AUROC<0.5)] <- 0.5#1-AUROC[which(AUROC<0.5)]
      return(AUROC)
    }
    cvm <- apply(predmat, 2, C_fun, y_true=y)
    cvsd <- rep(0, length(cvm))
  }
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}
