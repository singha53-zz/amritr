#' GLM single biomarker
#'
#' build a glm biomarker
#' @param X nxp matrix - training dataset
#' @param Y categorical variables
#' @param X.test nxp matrx  - test dataset
#' @param Y.test lablel of test dataset
#' @export
glmPanel = function(X, Y, X.test = NULL, Y.test = NULL){
  library(pROC); library(OptimalCutpoints);
  fit <- glm(Y ~ ., data = as.data.frame(X), family = binomial(link='logit'))
  if(!is.null(X.test)){
    probs <- predict(fit, as.data.frame(X.test))
    perfTest <- amritr::tperformance(weights = as.numeric(as.matrix(probs)),
      trueLabels = Y.test)
  } else {
    perfTest <- NA
  }

  return(list(X = X, Y =Y, X.test = X.test, Y.test = Y.test, perfTest = perfTest))
}

#' interal function (glm cross-validation)
#'
#' Estimate test error of glm biomarker
#' @param X nxp matrix - training dataset
#' @param Y categorical variables
#' @param M = # of folds
#' @param folds list of length M specifiying which samples are in which fold
#' @param progressBar (TRUE/FALSE)
#' @export
glmCV = function(X, Y, M, folds, progressBar){
  library(pROC); library(OptimalCutpoints)

  probs <- predictResponseList <- enet.panel <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    omit = folds[[i]]
    X.train = X[-omit, , drop = FALSE]
    Y.train = Y[-omit]
    X.test = X[omit, colnames(X.train), drop = FALSE]
    fit <- glm(Y.train ~ ., data = as.data.frame(X.train), family = binomial(link='logit'))
    probs[[i]] <- predict(fit, as.data.frame(X.test))
  }
  Y.test <- Y[unlist(folds)]
  perf <- amritr::tperformance(weights = as.numeric(as.matrix(unlist(probs))), trueLabels = Y.test)
  return(list(perf = perf))
}

#' cross-validation function for glm
#'
#' Estimate test error of glm biomarker
#' @param object = glmPanel object
#' @param validation - Mfold or LOOCV
#' @param M = # of folds
#' @param iter = number of times to repeat cross-valiation
#' @param threads - number of cpus (each cross-valiation scheme performed on a separate node)
#' @param progressBar - show progressbar (TRUE/FALE)
#' @export
perf.glm = function(object, validation = c("Mfold", "loo"), M = 5, iter = 10, threads = 4, progressBar = TRUE){
  library(dplyr); library(tidyr)
  X = object$X
  Y = object$Y
  n = nrow(X)

  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) createFolds(Y, k = M))
    require(parallel)
    cl <- parallel::makeCluster(mc <- getOption("cl.cores",
      threads))
    parallel::clusterExport(cl, varlist = c("glmPanel", "glmCV",
      "X", "Y", "M", "folds", "progressBar"), envir = environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi,
      X, Y, M, progressBar) {
      glmCV(X = X, Y = Y, M = M, folds = foldsi, progressBar = progressBar)
    }, X, Y, M, progressBar) %>% amritr::zip_nPure()
    parallel::stopCluster(cl)
    perf <- do.call(rbind, cv$perf) %>% as.data.frame %>%
      gather(ErrName, Err) %>% dplyr::group_by(ErrName) %>%
      dplyr::summarise(Mean = mean(Err), SD = sd(Err))
  }
  else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- glmCV(X, Y, M, folds, progressBar)
    perf <- data.frame(Mean = cv$perf) %>% mutate(ErrName = rownames(.))
    perf$SD <- NA
  }
  result = list()
  result$perf = perf
  method = "glm.mthd"
  result$meth = "glm.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}
