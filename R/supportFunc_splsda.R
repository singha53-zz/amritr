#' sPLSDA function
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param X.train n1xp dataset (discovery cohort)
#' @param Y.train vector of phenotype labels with names(Y.train) == rownames(X.train)
#' @param keepX (# of variables to select per component)
#' @param ncomp number of components
#' @param X.test n2xp dataset (validaton cohort)
#' @param Y.test vector of phenotype labels with names(Y.test) == rownames(X.test)
#' @param filter apply no filter "none" or a p.value filter "p.value
#' @param topranked select the top significant variables based on limma
#' @export
sPLSDA = function(X.train, Y.train, keepX, ncomp, X.test = NULL, Y.test = NULL,
  filter = "p.value", topranked = 50){
  if (filter == "none") {
    X.train1 <- X.train
  }
  if (filter == "p.value") {
    design <- model.matrix(~Y.train)
    fit <- eBayes(lmFit(t(X.train), design))
    top <- topTable(fit, coef = 2, adjust.method = "BH",
      n = nrow(fit))
    X.train1 <- X.train[, rownames(top)[1:topranked], drop = FALSE]
  }
  fit <- mixOmics::splsda(X = X.train1, Y = Y.train, keepX = keepX,
    ncomp = ncomp)
  panel <- unique(unlist(lapply(1:ncomp, function(i){
    names(which(fit$loadings$X[, i] != 0))
  })))
  if (!is.null(X.test)) {
    predMod <- predict(object = fit, newdata = X.test[, colnames(X.train1), drop = FALSE],
      dist = "all")
    if(length(Y.test) > 1){
      predictResponse <- lapply(predMod$class, function(i) {
        factor(i[, ncomp], levels = levels(Y.test))
      })
      errorRate <- lapply(predMod$class, function(i) {
        predictResponse <- factor(i[, ncomp], levels = levels(Y.test))
        trueLabels = Y.test
        mat <- table(trueLabels, predictResponse)
        mat2 <- mat
        diag(mat2) <- 0
        classError <- colSums(mat2)/colSums(mat)
        er <- sum(mat2)/sum(mat)
        ber <- mean(classError)
        perf <- c(classError, er, ber)
        names(perf) <- c(names(classError), "ER", "BER")
        perf
      }) %>% do.call(rbind, .) %>% as.data.frame %>% mutate(Method = rownames(.)) %>%
        tidyr::gather(Index, Value, -Method)
      probs <- predMod$predict[, , ncomp][, levels(Y)[2]]
      perfTest <- amritr::tperformance(weights = as.numeric(as.matrix(probs)),
        trueLabels = Y.test)
    } else {
      predictResponse <- lapply(predMod$class, function(i) {
        factor(i[, ncomp], levels = levels(Y.test))
      })
      probs <- predMod$predict[, , ncomp][levels(Y)[2]]
      perfTest <- errorRate <- NA
    }
  } else {
    predictResponse <- probs <- perfTest <- errorRate <- NA
  }
  return(list(X.train = X.train, Y.train = Y.train, fit = fit,
    panel = panel, predictResponse = predictResponse, probs = probs,
    perfTest = perfTest, ncomp = ncomp, keepX = keepX, errorRate = errorRate,
    filter = filter, topranked = topranked))
}

#' CV function for a sPLSDA model
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param X nxp dataset
#' @param Y vector of phenotype labels with names(Y) == rownames(X)
#' @param keepX (# of variables to select per component)
#' @param ncomp number of components
#' @param M Number of folds in the cross-validation
#' @param folds list of elements (sample indices) for the M folds
#' @param progressBar display progress bar or not (TRUE/FALSE)
#' @param filter apply no filter "none" or a p.value filter "p.value
#' @param topranked select the top significant variables based on limma
#' @export
sPLSDACV = function (X, Y, keepX, ncomp, M, folds, progressBar, filter,
  topranked){
  library(pROC)
  library(limma)
  library(mixOmics)
  library(OptimalCutpoints)
  probsList <- predictResponseList <- panelList <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    omit = folds[[i]]
    X.train = X[-omit, , drop = FALSE]
    Y.train = Y[-omit]
    if (filter == "none") {
      X.train1 <- X.train
    }
    if (filter == "p.value") {
      design <- model.matrix(~Y.train)
      fit <- eBayes(lmFit(t(X.train), design))
      top <- topTable(fit, coef = 2, adjust.method = "BH",
        n = nrow(fit))
      X.train1 <- X.train[, rownames(top)[1:topranked],
        drop = FALSE]
    }
    X.test1 = X[omit, colnames(X.train1), drop = FALSE]
    fit <- mixOmics::splsda(X = X.train1, Y = Y.train, keepX = keepX,
      ncomp = ncomp)
    panelList[[i]] <- unique(as.character(as.matrix(apply(fit$loadings$X,
      2, function(i) names(which(i != 0))))))
    predMod <- predict(object = fit, newdata = X.test1, dist = "all")
    predictResponseList[[i]] <- lapply(predMod$class, function(i) i[,
      ncomp, drop = FALSE])
    if (M == nrow(X)) {
      probsList[[i]] <- predMod$predict[, , ncomp][levels(Y)[2]]
    }
    else {
      probsList[[i]] <- predMod$predict[, , ncomp][, levels(Y)[2]]
    }
  }
  err <- predictResponseList %>% zip_nPure %>% lapply(., function(i) {
    predictResponse <- factor(unlist(i), levels = levels(Y))
    trueLabels = Y[unlist(folds)]
    mat <- table(trueLabels, predictResponse)
    mat2 <- mat
    diag(mat2) <- 0
    classError <- colSums(mat2)/colSums(mat)
    er <- sum(mat2)/sum(mat)
    ber <- mean(classError)
    perf <- c(classError, er, ber)
    names(perf) <- c(names(classError), "ER", "BER")
    perf
  })
  errorRate <- do.call(rbind, err) %>% as.data.frame %>% mutate(Method = rownames(.)) %>%
    tidyr::gather(Index, Value, -Method)
  probs <- unlist(probsList)
  trueLabels = Y[unlist(folds)]
  library(pROC)
  library(OptimalCutpoints)
  perf <- amritr::tperformance(weights = probs, trueLabels = trueLabels)
  panelFreq <- unlist(panelList)
  return(list(predictResponseList = predictResponseList, probsList = probsList,
    panelFreq = panelFreq, errorRate = errorRate, perf = perf))
}

#' repeated CV function for a sPLSDA model
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param object sPLSDA model object
#' @param valdation type of cross-validation; Mfold "Mfold" or LOOCV "loo"
#' @param M Number of folds in the cross-validation
#' @param iter number of times to repeat the cross-validation
#' @param threads number of cores to use to parrallel the CV
#' @param progressBar display progress bar or not (TRUE/FALSE)
#' @export
perf.sPLSDA = function (object, validation = c("Mfold", "loo"), M = 5, iter = 10,
  threads = 4, progressBar = TRUE)
{
  library(dplyr)
  library(tidyr)
  X = object$X.train
  Y = object$Y.train
  n = nrow(X)
  keepX = object$keepX
  ncomp = object$ncomp
  filter = object$filter
  topranked = object$topranked
  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) createFolds(Y, k = M))
    require(parallel)
    cl <- parallel::makeCluster(mc <- getOption("cl.cores",
      threads))
    parallel::clusterExport(cl, varlist = c("sPLSDACV", "X",
      "Y", "keepX", "ncomp", "M", "folds", "progressBar",
      "filter", "topranked"), envir = environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi,
      X, Y, keepX, ncomp, M, progressBar, filter, topranked) {
      library(mixOmics)
      library(dplyr)
      library(amritr)
      sPLSDACV(X = X, Y = Y, keepX = keepX, ncomp = ncomp,
        M = M, folds = foldsi, progressBar = progressBar,
        filter = filter, topranked = topranked)
    }, X, Y, keepX, ncomp, M, progressBar, filter, topranked) %>%
      amritr::zip_nPure()
    parallel::stopCluster(cl)
    panelFreq <- table(unlist(cv$panelFreq))[order(table(unlist(cv$panelFreq)),
      decreasing = TRUE)]/(M * iter)
    errorRate <- do.call(rbind, cv$errorRate) %>% dplyr::group_by(Method,
      Index) %>% dplyr::summarise(Mean = mean(Value), SD = sd(Value))
    perf <- do.call(rbind, cv$perf) %>% as.data.frame %>%
      gather(ErrName, Err) %>% dplyr::group_by(ErrName) %>%
      dplyr::summarise(Mean = mean(Err), SD = sd(Err))
  }
  else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- sPLSDACV(X = X, Y = Y, keepX = keepX, ncomp = ncomp,
      M = M, folds = folds, progressBar = progressBar,
      filter = filter, topranked = topranked)
    panelFreq <- table(unlist(cv$panelFreq))[order(table(unlist(cv$panelFreq)),
      decreasing = TRUE)]/M
    errorRate <- cv$errorRate %>% as.data.frame %>% dplyr::group_by(Method,
      Index) %>% dplyr::summarise(Mean = mean(Value), SD = sd(Value))
    perf <- data.frame(ErrName = names(cv$perf), Mean = as.numeric(cv$perf),
      SD = NA)
  }
  if (nlevels(Y) == 2) {
    if (validation == "Mfold") {
      colPalette <- c("#66C2A5", "#FC8D62")
      probDat <- lapply(cv$probsList, unlist) %>% lapply(.,
        function(i) i[names(Y)]) %>% do.call(rbind, .) %>%
        as.data.frame %>% gather(Sample, Prob) %>% mutate(phenotype = Y[Sample]) %>%
        group_by(phenotype, Sample) %>% dplyr::summarise(Mean = mean(Prob),
          SD = sd(Prob)) %>% arrange(Mean) %>% mutate(Sample = factor(as.character(Sample),
            levels = as.character(Sample)))
      probPlot <- probDat %>% ggplot(aes(x = Sample, y = Mean,
        color = phenotype)) + geom_point() + geom_rect(aes(xmin = -Inf,
          xmax = Inf, ymin = -Inf, ymax = subset(perf,
            ErrName == "cutoff")$Mean), fill = colPalette[1],
          alpha = 0.01) + geom_rect(aes(xmin = -Inf, xmax = Inf,
            ymin = subset(perf, ErrName == "cutoff")$Mean,
            ymax = Inf), fill = colPalette[2], alpha = 0.01) +
        geom_point() + geom_errorbar(aes(ymin = Mean -
            SD, ymax = Mean + SD)) + geom_hline(yintercept = subset(perf,
              ErrName == "cutoff")$Mean, linetype = "dashed") +
        customTheme(sizeStripFont = 10, xAngle = 90,
          hjust = 1, vjust = 0.5, xSize = 10, ySize = 10,
          xAxisSize = 10, yAxisSize = 10) + ylab(paste0("AUC Mean+/-SD ",
            "(", iter, "x", M, "-fold CV)")) + xlab("Observations") +
        ggtitle(paste(paste(perf$ErrName, round(perf$Mean,
          2), sep = "="), collapse = ", ")) + scale_color_manual(values = colPalette,
            name = "True labels")
    } else {
      colPalette <- c("#66C2A5", "#FC8D62")
      probDat <- data.frame(Mean = unlist(cv$probsList),
        Sample = names(Y)[unlist(folds)], SD = rep(NA,
          length(folds)), phenotype = Y[unlist(folds)]) %>%
        arrange(phenotype, Mean) %>% mutate(Sample = factor(as.character(Sample),
          levels = as.character(Sample)))
      probPlot <- probDat %>% ggplot(aes(x = Sample, y = Mean,
        color = phenotype)) + geom_point() + geom_rect(aes(xmin = -Inf,
          xmax = Inf, ymin = -Inf, ymax = subset(perf,
            ErrName == "cutoff")$Mean), fill = colPalette[1],
          alpha = 0.01) + geom_rect(aes(xmin = -Inf, xmax = Inf,
            ymin = subset(perf, ErrName == "cutoff")$Mean,
            ymax = Inf), fill = colPalette[2], alpha = 0.01) +
        geom_point() + geom_errorbar(aes(ymin = Mean -
            SD, ymax = Mean + SD)) + geom_hline(yintercept = subset(perf,
              ErrName == "cutoff")$Mean, linetype = "dashed") +
        customTheme(sizeStripFont = 10, xAngle = 90,
          hjust = 1, vjust = 0.5, xSize = 10, ySize = 10,
          xAxisSize = 10, yAxisSize = 10) + ylab("AUC (LOOCV)") + xlab("Observations") +
        ggtitle(paste(paste(perf$ErrName, round(perf$Mean,
          2), sep = "="), collapse = ", ")) + scale_color_manual(values = colPalette,
            name = "True labels")
    }
  }
  result = list()
  result$predictResponseList = cv$predictResponseList
  result$probsList = cv$probsList
  result$probsList = cv$probsList
  result$panelFreq = panelFreq
  result$errorRate = errorRate
  result$perf = perf
  if (nlevels(Y) == 2) {
    result$probPlot = probPlot
  }
  method = "sPLSDA.mthd"
  result$meth = "sPLSDA.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}

#' sPLSDA model after tuning of the number of variables
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param X nxp dataset
#' @param Y vector of phenotype labels with names(Y) == rownames(X)
#' @param ncomp number of components
#' @param keepXgrid sequence of integers (# of variables to select per component)
#' @param validatoin "Mfold" or "loo"
#' @param M Number of folds in the cross-validation
#' @export
tuned.sPLSDA = function(X.train, Y.train, keepXgrid, ncomp, X.test = X.test,
  Y.test = Y.test, filter = filter, topranked = topranked,
  validation = validation, M = M, iter = iter, threads = threads,
  progressBar = progressBar, optimal = optimal, errorMethod = errorMethod){
  indx <- unlist(lapply(keepXgrid, function(i) {
    result <- sPLSDA(X.train, Y.train, keepX = rep(i, ncomp),
      ncomp, X.test = X.test, Y.test = Y.test, filter = filter,
      topranked = topranked)
    cv <- perf.sPLSDA(object = result, validation = validation,
      M = M, iter = iter, threads = threads, progressBar = progressBar)

    if(optimal == "auc"){
      subset(cv$perf, ErrName == "AUC")$Mean
    } else {
      subset(cv$errorRate, Method == errorMethod & Index == "BER")$Mean
    }
  }))
  if(optimal == "auc"){
    keepX <- rep(keepXgrid[which(indx == max(indx))][1], 2)
  } else {
    keepX <- rep(keepXgrid[which(indx == min(indx))][1], 2)
  }

  fit <- sPLSDA(X.train, Y.train, keepX = keepX, ncomp, X.test = X.test,
    Y.test = Y.test, filter = filter, topranked = topranked)
  panel <- fit$panel
  predictResponse <- fit$predictResponse
  probs <- fit$probs
  perfTest <- fit$perfTest
  errorRate <- fit$errorRate
  return(list(X.train = X.train, Y.train = Y.train, keepXgrid = keepXgrid,
    panel = panel, predictResponse = predictResponse, probs = probs,
    perfTest = perfTest, ncomp = ncomp, keepX = keepX, errorRate = errorRate,
    filter = filter, topranked = topranked, validation = validation,
    M = M, threads = threads, progressBar = progressBar,
    iter = iter, optimal = optimal, errorMethod = errorMethod))
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
tuned.sPLSDACV = function(object, folds){
  X <- object$X.train
  Y <- object$Y.train
  keepXgrid <- object$keepXgrid
  ncomp <- object$ncomp
  filter = object$filter
  topranked = object$topranked
  validation = object$validation
  M = length(folds)
  iter = object$iter
  threads = object$threads
  progressBar = object$progressBar
  optimal = object$optimal
  errorMethod = object$errorMethod
  library(pROC)
  library(limma)
  library(mixOmics)
  library(OptimalCutpoints)
  probsList <- predictResponseList <- panelList <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    omit = folds[[i]]
    X.train = X[-omit, , drop = FALSE]
    Y.train = Y[-omit]
    if (filter == "none") {
      X.train1 <- X.train
    }
    if (filter == "p.value") {
      design <- model.matrix(~Y.train)
      fit <- eBayes(lmFit(t(X.train), design))
      top <- topTable(fit, coef = 2, adjust.method = "BH",
        n = nrow(fit))
      X.train1 <- X.train[, rownames(top)[1:topranked],
        drop = FALSE]
    }
    X.test1 = X[omit, colnames(X.train1), drop = FALSE]
    fit <- tuned.sPLSDA(X.train = X.train1, Y.train = Y.train,
      keepXgrid, ncomp, X.test = X.test1, Y.test = Y[omit],
      filter = filter, topranked = topranked, validation = validation,
      M = M, iter = iter, threads = threads, progressBar = progressBar,
      optimal = optimal, errorMethod = errorMethod)
    panelList[[i]] = fit$panel
    predictResponseList[[i]] <- fit$predictResponse
    probsList[[i]] <- fit$probs
  }
  err <- predictResponseList %>% zip_nPure %>% lapply(., function(i) {
    predictResponse <- factor(unlist(i), levels = levels(Y))
    trueLabels = Y[unlist(folds)]
    mat <- table(trueLabels, predictResponse)
    mat2 <- mat
    diag(mat2) <- 0
    classError <- colSums(mat2)/colSums(mat)
    er <- sum(mat2)/sum(mat)
    ber <- mean(classError)
    perf <- c(classError, er, ber)
    names(perf) <- c(names(classError), "ER", "BER")
    perf
  })
  errorRate <- do.call(rbind, err) %>% as.data.frame %>% mutate(Method = rownames(.)) %>%
    tidyr::gather(Index, Value, -Method)
  probs <- unlist(probsList)
  trueLabels = Y[unlist(folds)]
  library(pROC)
  library(OptimalCutpoints)
  perf <- amritr::tperformance(weights = probs, trueLabels = trueLabels)
  if (validation == "Mfold") {
    panelFreq <- table(unlist(panelList))/(M * iter)
    panelFreq <- panelFreq[order(panelFreq, decreasing = TRUE)]
  } else {
    panelFreq <- table(unlist(panelList))/(M)
    panelFreq <- panelFreq[order(panelFreq, decreasing = TRUE)]
  }
  return(list(predictResponseList = predictResponseList, probsList = probsList,
    panelFreq = panelFreq, errorRate = errorRate, perf = perf))
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
perf.tuned.sPLSDA = function(object){
  library(dplyr)
  library(tidyr)
  Y <- object$Y.train
  n <- length(Y)
  validation <- object$validation
  iter <- object$iter
  if (validation == "Mfold") {
    M <- object$M
    folds <- lapply(1:iter, function(i) createFolds(Y, k = M))
    cv <- list()
    for (i in 1:iter) {
      cv[[i]] <- tuned.sPLSDACV(object = object, folds = folds[[i]])
    }
    cv <- cv %>% amritr::zip_nPure()
    panelFreq <- table(unlist(cv$panelFreq))[order(table(unlist(cv$panelFreq)),
      decreasing = TRUE)]/(M * iter)
    errorRate <- do.call(rbind, cv$errorRate) %>% dplyr::group_by(Method,
      Index) %>% dplyr::summarise(Mean = mean(Value), SD = sd(Value))
    perf <- do.call(rbind, cv$perf) %>% as.data.frame %>%
      gather(ErrName, Err) %>% dplyr::group_by(ErrName) %>%
      dplyr::summarise(Mean = mean(Err), SD = sd(Err))
  } else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- tuned.sPLSDACV(object = object, folds = folds)
    panelFreq <- cv$panelFreq
    errorRate <- cv$errorRate
    perf <- data.frame(ErrName = names(cv$perf), Mean = as.numeric(cv$perf),
      SD = NA)
  }
  if (nlevels(Y) == 2) {
    if (validation == "Mfold") {
      colPalette <- c("#66C2A5", "#FC8D62")
      probDat <- lapply(cv$probsList, unlist) %>% lapply(.,
        function(i) i[names(Y)]) %>% do.call(rbind, .) %>%
        as.data.frame %>% gather(Sample, Prob) %>% mutate(phenotype = Y[Sample]) %>%
        group_by(phenotype, Sample) %>% dplyr::summarise(Mean = mean(Prob),
          SD = sd(Prob)) %>% arrange(Mean) %>% mutate(Sample = factor(as.character(Sample),
            levels = as.character(Sample)))
      probPlot <- probDat %>% ggplot(aes(x = Sample, y = Mean,
        color = phenotype)) + geom_point() + geom_rect(aes(xmin = -Inf,
          xmax = Inf, ymin = -Inf, ymax = subset(perf,
            ErrName == "cutoff")$Mean), fill = colPalette[1],
          alpha = 0.01) + geom_rect(aes(xmin = -Inf, xmax = Inf,
            ymin = subset(perf, ErrName == "cutoff")$Mean,
            ymax = Inf), fill = colPalette[2], alpha = 0.01) +
        geom_point() + geom_errorbar(aes(ymin = Mean -
            SD, ymax = Mean + SD)) + geom_hline(yintercept = subset(perf,
              ErrName == "cutoff")$Mean, linetype = "dashed") +
        customTheme(sizeStripFont = 10, xAngle = 90,
          hjust = 1, vjust = 0.5, xSize = 10, ySize = 10,
          xAxisSize = 10, yAxisSize = 10) + ylab(paste0("AUC Mean+/-SD ",
            "(", iter, "x", M, "-fold CV)")) + xlab("Observations") +
        ggtitle(paste(paste(perf$ErrName, round(perf$Mean,
          2), sep = "="), collapse = ", ")) + scale_color_manual(values = colPalette,
            name = "True labels")
    } else {
      colPalette <- c("#66C2A5", "#FC8D62")
      probDat <- data.frame(Mean = unlist(cv$probsList),
        Sample = names(Y)[unlist(folds)], SD = rep(NA,
          length(folds)), phenotype = Y[unlist(folds)]) %>%
        arrange(phenotype, Mean) %>% mutate(Sample = factor(as.character(Sample),
          levels = as.character(Sample)))
      probPlot <- probDat %>% ggplot(aes(x = Sample, y = Mean,
        color = phenotype)) + geom_point() + geom_rect(aes(xmin = -Inf,
          xmax = Inf, ymin = -Inf, ymax = subset(perf,
            ErrName == "cutoff")$Mean), fill = colPalette[1],
          alpha = 0.01) + geom_rect(aes(xmin = -Inf, xmax = Inf,
            ymin = subset(perf, ErrName == "cutoff")$Mean,
            ymax = Inf), fill = colPalette[2], alpha = 0.01) +
        geom_point() + geom_errorbar(aes(ymin = Mean -
            SD, ymax = Mean + SD)) + geom_hline(yintercept = subset(perf,
              ErrName == "cutoff")$Mean, linetype = "dashed") +
        customTheme(sizeStripFont = 10, xAngle = 90,
          hjust = 1, vjust = 0.5, xSize = 10, ySize = 10,
          xAxisSize = 10, yAxisSize = 10) + ylab("AUC (LOOCV)") +
        xlab("Observations") + ggtitle(paste(paste(perf$ErrName,
          round(perf$Mean, 2), sep = "="), collapse = ", ")) +
        scale_color_manual(values = colPalette, name = "True labels")
    }
  }
  result = list()
  result$predictResponseList = cv$predictResponseList
  result$probsList = cv$probsList
  result$probsList = cv$probsList
  result$panelFreq = panelFreq
  result$errorRate = errorRate
  result$perf = perf
  if (nlevels(Y) == 2) {
    result$probPlot = probPlot
  }
  method = "sPLSDA.mthd"
  result$meth = "sPLSDA.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}
