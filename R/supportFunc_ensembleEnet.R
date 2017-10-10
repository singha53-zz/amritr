#' Build ensemble classification panel
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param X.train - list of training datasets (nxpi); i number of elements
#' @param Y.train - n-vector of class labels
#' @param alpha = list of alpha values
#' @param lambda = list of lambda values
#' @param X.test - list of test datasets (nxpi); i number of elements
#' @param Y.test - n-vector of class labels
#' @export
ensembleEnet = function (X.trainList, Y.train, alphaList, lambdaList, family = "binomial",
  X.testList, Y.test, filter, topranked, keepVarList){

  ## load libraries
  library(glmnet); library(limma); library(pROC);
  library(OptimalCutpoints); library(tidyverse);

  ## perform pre-filtering (none, p-value, and keep certain variables)
  if (filter == "none") {
    X1.trainList <- X.trainList
  }
  if (filter == "p.value") {
    X1.trainList <- lapply(X.trainList, function(i){
      design <- model.matrix(~Y.train)
      fit <- eBayes(lmFit(t(i), design))
      top <- topTable(fit, coef = 2, adjust.method = "BH", n = nrow(fit))
      i[, rownames(top)[1:topranked]]
    })
  }
  if (is.null(keepVarList)) {
    penalty.factorList <- lapply(X1.trainList, function(i){rep(1, ncol(i))})
    X2.trainList  <- X1.trainList
  } else {
    X2.trainList <- mapply(function(x, x1, y){
      X1 <- x1[, setdiff(colnames(x1), y)]
      X2 <- as.matrix(cbind(X1, x[, y]))
      colnames(X2) <- c(colnames(X1), y)
      X2
    }, x = X.trainList, x1 = X1.trainList, y = keepVarList)

    penalty.factorList <- mapply(function(x, y){
      c(rep(1, ncol(X1)), rep(0, length(keepVar)))
    }, x = X1.trainList, y = keepVarList)
  }

  ## build glmnet classifier
  model <- mapply(function(X, alpha, lambda, penalty.factor){
    if(family == "binomial") {
      fit <- glmnet(X, Y.train, family = "binomial", alpha = alpha,
        penalty.factor = penalty.factor)
      cv.fit <- cv.glmnet(X, Y.train, family = "binomial")
    } else {
      fit <- glmnet(X, Y.train, family = "multinomial", alpha = alpha,
        type.multinomial = "grouped", penalty.factor = penalty.factor)
      cv.fit <- cv.glmnet(X, Y.train, family = "multinomial")
    }
    if(is.null(lambda)) {lambda = cv.fit$lambda.min} else {lambda = lambda}
    Coefficients <- coef(fit, s = lambda)
    if(family == "binomial"){
      Active.Index <- which(Coefficients[, 1] != 0)
      Active.Coefficients <- Coefficients[Active.Index, ]
    } else {
      Active.Index <- which(Coefficients[[1]][, 1] != 0)
      Active.Coefficients <- Coefficients[[1]][Active.Index, ]
    }
    enet.panel <- names(Active.Coefficients)[-1]
    enet.panel.length <- length(enet.panel)
    return(list(fit=fit, Coefficients=Coefficients, Active.Index=Active.Index, lambda = lambda,
      Active.Coefficients=Active.Coefficients, enet.panel=enet.panel, enet.panel.length=enet.panel.length))
  }, X = X2.trainList, alpha = alphaList, lambda = lambdaList, penalty.factor = penalty.factorList, SIMPLIFY = FALSE)

  ## Test performance in test dataset
  if(!is.null(X.testList)){
    if(!all(sapply(1 : length(X.trainList), function(i) any(colnames(X.trainList[[i]]) == colnames(X.testList[[i]])))))
      stop("features of the train and test datasets are not in the same order")
    if(!any(levels(Y.train) == levels(Y.test)))
      stop("levels of Y.train and Y.test are not in the same order")

    testPerf <- mapply(function(mod, test){
      predictResponse <- unlist(predict(mod$fit, newx = test[, rownames(mod$Coefficients)[-1]], s = mod$lambda, type = "class"))
      probs <- predict(mod$fit, newx = test[, rownames(mod$Coefficients)[-1]], s = mod$lambda, type = "response") %>% as.numeric
      names(probs) <- rownames(predictResponse)
      predictResponse <- as.character(predictResponse)
      names(predictResponse) <- names(probs)

      ## compute error rate
      mat <- table(pred=factor(as.character(predictResponse), levels = levels(Y.train)), truth=Y.test)
      mat2 <- mat
      diag(mat2) <- 0
      classError <- colSums(mat2)/colSums(mat)
      er <- sum(mat2)/sum(mat)
      ber <- mean(classError)
      error <- c(classError, er, ber) %>% matrix
      rownames(error) <- c(names(classError), "Overall", "Balanced")
      colnames(error) <- "Error_0.5"

      ## compute AUROC
      if(length(Y.test) > 1) {
        if(nlevels(Y.train) == 2){
          Y.test <- factor(as.character(Y.test), levels(Y.train))
          perfTest <- amritr::tperformance(weights = as.numeric(as.matrix(probs)), trueLabels = Y.test) %>% as.matrix
          colnames(perfTest) <- paste(levels(Y.train), collapse = "_vs_")
        } else {
          perfTest <- NA
        }
      } else {
        perfTest <- NA
      }

      return(list(probs=probs, predictResponse=predictResponse, error=error, perfTest=perfTest))
    }, mod = model, test = X.testList, SIMPLIFY = FALSE)
  } else {testPerf <- NA}

  return(list(model=model, testPerf=testPerf, X.trainList=X.trainList,
    Y.train=Y.train, alphaList=alphaList, lambdaList=lambdaList, family=family, X.testList=X.testList,
    Y.test=Y.test, filter=filter, topranked=topranked, keepVarList=keepVarList))
}


#' Estimate test error using repeated cross-validation
#'
#'
#' @param object - ensembleEnet object
#' @param validation = Mfold or loo
#' @param M - # of folds
#' @param iter - Number of iterations of cross-validation
#' @param threads - # of cores, running each iteration on a separate node
#' @param progressBar = TRUE (show progress bar or not)
#' @export
perfEnsemble = function(object, validation = "Mfold", M = 5, iter = 5, threads = 5, progressBar = TRUE){
  library(tidyverse)
  X.trainList=object$X.trainList
  Y.train=object$Y.train
  alphaList=object$alphaList
  lambdaList=object$lambdaList
  family=object$family
  filter=object$filter
  topranked=object$topranked
  keepVarList=object$keepVarList

  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) createFolds(Y.train, k = M))
    require(parallel)
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", threads))
    parallel::clusterExport(cl, varlist = c("ensembleCV",
      "ensembleEnet", "X.trainList", "Y.train", "alphaList", "lambdaList",
      "family", "filter", "topranked", "keepVarList",
      "M", "folds", "progressBar"),
      envir = environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi, X.trainList, Y.train, alphaList, lambdaList, family, filter, topranked, keepVarList, M, progressBar) {
      ensembleCV(X.trainList=X.trainList, Y.train=Y.train, alphaList=alphaList, lambdaList=lambdaList, family=family, filter=filter, topranked=topranked, keepVarList=keepVarList, M=M, folds=foldsi, progressBar=progressBar)
    }, X.trainList, Y.train, alphaList, lambdaList, family, filter, topranked,
      keepVarList, M, progressBar) %>%
      amritr::zip_nPure()
    parallel::stopCluster(cl)
    error <- do.call(rbind, cv$error) %>% as.data.frame %>%
      mutate(ErrName = factor(rownames(.), unique(rownames(.)))) %>%
      dplyr::group_by(ErrName) %>%
      dplyr::summarise(Mean = mean(Error_0.5), SD = sd(Error_0.5))
    perfTest <- do.call(rbind, cv$perfTest) %>% as.data.frame %>%
      mutate(ErrName = factor(rownames(.), unique(rownames(.)))) %>%
      dplyr::group_by(ErrName) %>%
      dplyr::summarise(Mean = mean(perf), SD = sd(perf))
  }
  else {
    n <- length(Y.train)
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- ensembleCV(X.trainList, Y.train, alphaList, lambdaList, family, filter, topranked,
      keepVarList, M, folds, progressBar)
    error <- cv$error
    perfTest <- cv$perfTest
  }
  result = list()
  result$error = error
  result$perfTest = perfTest
  method = "enetEnsemble.mthd"
  result$meth = "enetEnsemble.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}

#' Estimate cross-validation using cross-validation
#'
#'
#' @param X - list of training datasets (nxpi); i number of elements
#' @param Y - n-vector of class labels
#' @param alpha = list of alpha values
#' @param lambda = list of lambda values
#' @param M - # of folds
#' @param folds - list of length M, where each element contains the indices for samples for a given fold
#' @param progressBar (TRUE/FALSE) - show progress bar or not
#' @param filter - "none" or "p.value"
#' @param topranked - # of significant features to use to build classifier
#' @export
ensembleCV = function (X.trainList, Y.train, alphaList, lambdaList, family, filter, topranked,
  keepVarList, M, folds, progressBar){
  J <- length(X.trainList)
  assign("X.training", NULL, pos = 1)
  assign("Y.training", NULL, pos = 1)
  X.training = lapply(folds, function(x) {
    lapply(1:J, function(y) {
      X.trainList[[y]][-x, , drop = FALSE]
    })
  })
  Y.training = lapply(folds, function(x) {
    Y.train[-x]
  })
  X.test = lapply(folds, function(x) {
    lapply(1:J, function(y) {
      X.trainList[[y]][x, , drop = FALSE]
    })
  })
  Y.test = lapply(folds, function(x) {
    Y.train[x]
  })
  avgProbList <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    ## build ensemble panel
    result <- ensembleEnet(X.trainList=X.training[[i]], Y.train=Y.training[[i]],
      alphaList, lambdaList, family = family,
      X.testList=X.test[[i]], Y.test=Y.test[[i]], filter, topranked, keepVarList)
    # combine predictions using average probability
    avgProbList[[i]] <- do.call(cbind, lapply(result$testPerf,
      function(i) {
        i$probs
      })) %>% rowMeans
  }

  probs <- unlist(avgProbList)
  ## Error and AUROC
  predictResponse <- rep(levels(Y.train)[1], length(probs))
  predictResponse[probs >= 0.5] <- levels(Y.train)[2]

  ## compute error rate
  truth <- sapply(strsplit(names(unlist(Y.test)), "\\."), function(i) i[2])
  if(!all(names(probs) == truth))
    stop("predicted probability is not in the same order as the test labels")
  mat <- table(pred=factor(as.character(predictResponse), levels = levels(Y.train)), truth=unlist(Y.test))
  mat2 <- mat
  diag(mat2) <- 0
  classError <- colSums(mat2)/colSums(mat)
  er <- sum(mat2)/sum(mat)
  ber <- mean(classError)
  error <- c(classError, er, ber) %>% matrix
  rownames(error) <- c(names(classError), "Overall", "Balanced")
  colnames(error) <- "Error_0.5"

  ## compute AUROC
  if(nlevels(Y.train) == 2){
    perfTest <- amritr::tperformance(weights = probs, trueLabels = unlist(Y.test)) %>% as.matrix
    colnames(perfTest) <- "perf"
  } else {
    perfTest <- NA
  }

  return(list(error = error, perfTest = perfTest))
}
