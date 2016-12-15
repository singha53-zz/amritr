#' CV function for a splsda model
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
splsdaCV = function(X, Y, keepX, ncomp, M, folds, progressBar, filter, topranked){
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

    ## fit model to train dataset
    fit <- mixOmics::splsda(X = X.train1, Y = Y.train, keepX = keepX, ncomp = ncomp)
    panelList[[i]] <- unique(as.character(as.matrix(apply(fit$loadings$X, 2, function(i) names(which(i != 0))))))
    predMod <- predict(object = fit, newdata = X.test1, dist = "all")

    predictResponseList[[i]] <- lapply(predMod$class, function(i) i[, ncomp, drop = FALSE])
    if(M == nrow(X)){
      probsList[[i]] <- predMod$predict[,, ncomp][levels(Y)[2]]
    } else {
      probsList[[i]] <- predMod$predict[,, ncomp][, levels(Y)[2]]
    }
  }

  ## Error Rate per distance method
  err <- predictResponseList %>% zip_nPure %>% lapply(., function(i){
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

  ## AUC per distance method
  probs <- unlist(probsList)
  trueLabels = Y[unlist(folds)]
  library(pROC)
  library(OptimalCutpoints)
  perf <- amritr::tperformance(weights = probs, trueLabels = trueLabels)

  panelFreq <- unlist(panelList)

  return(list(predictResponseList = predictResponseList, probsList = probsList, panelFreq = panelFreq,
    errorRate = errorRate, perf = perf))
}


#' repeated CV function for a splsda model
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param object splsda model object
#' @param valdation type of cross-validation; Mfold "Mfold" or LOOCV "loo"
#' @param M Number of folds in the cross-validation
#' @param iter number of times to repeat the cross-validation
#' @param threads number of cores to use to parrallel the CV
#' @param progressBar display progress bar or not (TRUE/FALSE)
#' @export
perf.splsda = function (object, validation = c("Mfold", "loo"), M = 5, iter = 10, threads = 4, progressBar = TRUE){
  library(dplyr)
  library(tidyr)
  X = object$X
  Y = object$Y
  n = nrow(X)
  keepX = object$keepX
  ncomp = object$ncomp
  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) createFolds(Y, k = M))
    require(parallel)
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", threads))
    parallel::clusterExport(cl, varlist = c("splsdaCV", "X", "Y", "keepX", "ncomp",
      "M", "folds", "progressBar", "filter", "topranked"), envir = environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi,
      X, Y, keepX, ncomp, M, progressBar, filter, topranked) {
      library(mixOmics); library(dplyr); library(amritr);
      splsdaCV(X = X, Y = Y, keepX = keepX, ncomp = ncomp, M = M, folds = foldsi,
        progressBar = progressBar, filter = filter, topranked = topranked)
    }, X, Y, keepX, ncomp, M, progressBar, filter, topranked) %>% amritr::zip_nPure()
    parallel::stopCluster(cl)

    ## CV panels
    panelFreq <- table(unlist(cv$panelFreq))[order(table(unlist(cv$panelFreq)), decreasing = TRUE)]/(M*iter)

    ## error rate
    errorRate <- do.call(rbind, cv$errorRate) %>% dplyr::group_by(Method, Index) %>%
      dplyr::summarise(Mean = mean(Value), SD = sd(Value))

    ## classification performance measures
    perf <- do.call(rbind, cv$perf) %>% as.data.frame %>%
      gather(ErrName, Err) %>% dplyr::group_by(ErrName) %>%
      dplyr::summarise(Mean = mean(Err), SD = sd(Err))
  } else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- splsdaCV(X = X, Y = Y, keepX = keepX, ncomp = ncomp, M = M, folds = folds,
      progressBar = progressBar, filter = filter, topranked = topranked)
    ## CV panels
    panelFreq <- table(unlist(cv$panelFreq))[order(table(unlist(cv$panelFreq)), decreasing = TRUE)]/M

    ## error rate
    errorRate <- cv$errorRate %>% as.data.frame %>%
      dplyr::group_by(Method, Index) %>%
      dplyr::summarise(Mean = mean(Value), SD = sd(Value))

    ## classification performance measures
    perf <- data.frame(ErrName = names(cv$perf),
      Mean = as.numeric(cv$perf), SD = NA)
  }

  if(nlevels(Y) == 2){
    colPalette <- c("#66C2A5", "#FC8D62")
    ## Generate summarise plots
    ## plot prob vs. samples plot
    probDat <- lapply(cv$probsList, unlist) %>% lapply(., function(i) i[names(Y)]) %>%
      do.call(rbind, .) %>% as.data.frame %>% gather(Sample, Prob) %>%
      mutate(phenotype = Y[Sample]) %>%
      group_by(phenotype, Sample) %>% dplyr::summarise(Mean = mean(Prob), SD = sd(Prob)) %>%
      arrange(Mean) %>%
      mutate(Sample = factor(as.character(Sample), levels = as.character(Sample)))
    probPlot <- probDat %>% ggplot(aes(x = Sample, y = Mean, color = phenotype)) + geom_point() +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = subset(perf, ErrName == "cutoff")$Mean), fill = colPalette[1], alpha = 0.01) +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = subset(perf, ErrName == "cutoff")$Mean, ymax = Inf), fill = colPalette[2], alpha = 0.01) +
      geom_point() +
      geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD)) +
      geom_hline(yintercept = subset(perf, ErrName == "cutoff")$Mean, linetype = "dashed") +
      customTheme(sizeStripFont = 10, xAngle = 90, hjust = 1, vjust = 0.5, xSize = 10,
        ySize = 10, xAxisSize = 10, yAxisSize = 10) +
      ylab(paste0("AUC Mean+/-SD ", "(", iter, "x", M, "-fold CV)")) +
      xlab("Observations") +
      ggtitle(paste(paste(perf$ErrName, round(perf$Mean, 2), sep = "="),
        collapse=", ")) + scale_color_manual(values=colPalette, name = "True labels")
  }

  result = list()
  result$predictResponseList = cv$predictResponseList
  result$probsList  = cv$probsList
  result$probsList = cv$probsList
  result$panelFreq = panelFreq
  result$errorRate = errorRate
  result$perf = perf
  if(nlevels(Y) == 2){result$probPlot = probPlot}
  method = "splsda.mthd"
  result$meth = "splsda.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}





















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
