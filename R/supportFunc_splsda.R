#' splsda model after tuning of the number of variables
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param X nxp dataset
#' @param Y vector of phenotype labels with names(Y) == rownames(X)
#' @param ncomp number of components
#' @param keepXgrid sequence of integers (# of variables to select per component)
#' @param validatoin "Mfold" or "loo"
#' @param M Number of folds in the cross-validation
#' @export
tuned.spslda = function(X, Y, ncomp, keepXgrid, validation, M){
  library(dplyr);
  resultPerf <- lapply(keepXgrid, function(i){
    result <- mixOmics::splsda(X, Y, keepX = rep(i, ncomp), ncomp = ncomp)
    cv <- as.data.frame(mixOmics::perf(result, validation = validation, folds = M, progressBar = FALSE)$error.rate)
    cv$Comp <- paste("Comp1", 1:ncomp, sep = "-")
    cv %>% tidyr::gather(Type, ErrorRate, -Comp) %>%
      tidyr::separate(col = Type, into = c("Type", "Method"), sep = "\\.") %>%
      dplyr::mutate(keepX = rep(i, nrow(.)))
  })

  err <- do.call(rbind, resultPerf)  %>% as.data.frame %>%
    dplyr::filter(Comp == paste("Comp1", ncomp, sep = "-"), Type == "BER") %>%
    dplyr::group_by(Type, Method) %>%
    dplyr::filter(ErrorRate == min(ErrorRate)) %>% slice(1)

  return(list(X = X, Y = Y, ncomp = ncomp, keepXgrid = keepXgrid, validation = validation, M = M, err = err))
}

#' Runs a cross-validation scheme a given number (iter) of times
#'
#' estimate the error rate and AUC of tuned.splsda panel
#' @param object tuned.splsda model
#' @param validatoin "Mfold" or "loo"
#' @param M Number of folds in the cross-validation
#' @param iter Number of times to repeat cross-validation
#' @param progressBar = show progress bar or not
#' @export
perf.tuned.splsda=function (object, validation = "Mfold", M = 5, iter = 5, progressBar = TRUE){
  if(is.null(Y))
    stop("names(Y) should equal rownams(X)")
  for(i in 1 : nlevels(Y)){
    if(levels(Y)[i] == "ER")
      stop("levels of Y include ER or BER")
  }

  library(dplyr)
  library(tidyr)
  X = object$X
  Y = object$Y
  ncomp = object$ncomp
  keepXgrid = object$keepXgrid
  validation = object$validation
  M = object$M
  n = nrow(X)
  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) createFolds(Y, k = M))
    cvList <- list()
    for (i in 1:length(folds)) {
      cvList[[i]] <- tuned.spsldaCV(X = X, Y = Y, folds = folds[[i]],
        progressBar = progressBar, ncomp = ncomp, keepXgrid = keepXgrid,
        validation = validation, M = M)
    }
    cv <- cvList %>% amritr::zip_nPure()
    errRes <- cv$errorRate %>% amritr::zip_nPure() %>% lapply(.,
      function(j) {
        do.call(rbind, j)
      })
    errorRate <- do.call(rbind, errRes) %>% as.data.frame %>%
      mutate(Method = rep(names(errRes), each = nrow(errRes[[1]]))) %>%
      tidyr::gather(Index, Value, -Method) %>% dplyr::group_by(Index,
        Method) %>% dplyr::summarise(Mean = mean(Value),
          SD = sd(Value)) %>% ungroup %>% dplyr::group_by(Method) %>%
      dplyr::filter(Index == "BER")
    aucRes <- cv$auc %>% amritr::zip_nPure() %>% lapply(.,
      function(j) {
        do.call(rbind, j)
      })
    auc <- do.call(rbind, aucRes) %>% as.data.frame %>% mutate(Method = rep(names(errRes),
      each = nrow(errRes[[1]]))) %>% tidyr::gather(Index,
        Value, -Method) %>% dplyr::group_by(Index, Method) %>%
      dplyr::summarise(Mean = mean(Value), SD = sd(Value)) %>%
      ungroup
  } else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- tuned.spsldaCV(X = X, Y = Y, folds = folds, progressBar = progressBar,
      ncomp = ncomp, keepXgrid = keepXgrid, validation = validation,
      M = M)
    errorRate <- do.call(rbind, cv$errorRate) %>% as.data.frame %>%
      dplyr::mutate(Method = rownames(.)) %>%
      tidyr::gather(Index, Value, -Method) %>%
      dplyr::group_by(Index, Method) %>%
      dplyr::summarise(Mean = mean(Value), SD = sd(Value)) %>%
      ungroup
    auc <- do.call(rbind, cv$auc) %>% as.data.frame %>%
      dplyr::mutate(Method = rownames(.)) %>%
      tidyr::gather(Index, Value, -Method) %>%
      dplyr::group_by(Index, Method) %>%
      dplyr::summarise(Mean = mean(Value), SD = sd(Value)) %>%
      ungroup
  }
  result = list()
  result$folds = folds
  result$errorRate = errorRate
  result$auc = auc
  method = "tune.splsda.mthd"
  result$meth = "tune.splsda.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}

#' Runs a cross-validation scheme for a tuned.splsda model once
#'
#' estimate the error rate and AUC of tuned.splsda panel
#' @param X nxp dataset
#' @param Y vector of phenotype labels with names(Y) == rownames(X)
#' @param ncomp number of components
#' @param folds list of length M where each element contains the indices of the samples in that fold
#' @param progressBar = TRUE/FALSE
#' @param keepXgrid sequence of integers (# of variables to select per component)
#' @param validatoin "Mfold" or "loo"
#' @param M Number of folds in the cross-validation
#' @export
tuned.spsldaCV=function (X = X, Y = Y, ncomp = ncomp, folds = folds, progressBar = FALSE,
  keepXgrid = keepXgrid, validation = validation, M = M)
{
  library(mixOmics)
  probs <- predictResponseList <- panel <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    omit = folds[[i]]
    if (length(omit) == 1)
      stop.user = TRUE
    X.train = X[-omit, ]
    Y.train = Y[-omit]
    X.test = matrix(X[omit, ], nrow = length(omit))
    X.test = X[omit, , drop = FALSE]
    spls.res0 = tuned.spslda(X.train, Y.train, ncomp, keepXgrid,
      validation, M)
    testErr <- lapply(spls.res0$err$Method, function(z) {
      keepX <- rep(subset(spls.res0$err, Method == z)$keepX,
        ncomp)
      spls.res = mixOmics::splsda(X.train, Y.train, ncomp,
        keepX, mode = "regression")
      if (!is.null(spls.res$nzv$Position))
        X.test = X.test[, -spls.res$nzv$Position]
      predictResponseList <- predict(object = spls.res,
        newdata = X.test, dist = paste(z, "dist", sep = "."))$class
      probs <- predict(object = spls.res, newdata = X.test,
        dist = paste(z, "dist", sep = "."))$predict
      return(list(predictResponseList = predictResponseList,
        probs = probs))
    }) %>% zip_nPure()
    predictResponseList[[i]] <- testErr$predictResponseList
    probs[[i]] <- testErr$probs
  }
  errorRate <- predictResponseList %>% zip_nPure() %>% lapply(.,
    function(i) {
      predictedLabels <- do.call(rbind, lapply(i, function(j) j[[1]]))[,
        paste("comp", ncomp, sep = " ")]
      trueLabels = Y[names(predictedLabels)]
      mat <- table(trueLabels, predictedLabels)
      mat2 <- mat
      diag(mat2) <- 0
      classError <- colSums(mat2)/colSums(mat)
      er <- sum(mat2)/sum(mat)
      ber <- mean(classError)
      perf <- c(classError, er, ber)
      names(perf) <- c(names(classError), "ER", "BER")
      perf
    })
  if(validation == "Mfold"){
    auc <- probs %>% zip_nPure() %>% lapply(., function(i) {
      pb <- do.call(rbind, lapply(i, function(j) {
        j[, , paste("dim", ncomp, sep = " ")]
      }))
      prob <- unlist(pb[, ncomp])
      trueLabels = Y[rownames(pb)]
      library(pROC)
      library(OptimalCutpoints)
      amritr::tperformance(weights = prob, trueLabels = trueLabels)
    })
  } else {
    auc <- probs %>% zip_nPure() %>% lapply(., function(i) {
      pb <- do.call(rbind, lapply(i, function(j) {
        j[, , paste("dim", ncomp, sep = " ")]
      }))
      rownames(pb) <- rownames(X)
      prob <- unlist(pb[, ncomp])
      trueLabels = Y[rownames(pb)]
      library(pROC)
      library(OptimalCutpoints)
      amritr::tperformance(weights = prob, trueLabels = trueLabels)
    })
  }

  names(errorRate) <- names(auc) <- spls.res0$err$Method
  return(list(errorRate = errorRate, auc = auc))
}
