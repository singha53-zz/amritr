#' Build ensemble enet classification panel
#'
#' @param X.trainList - list of training datasets (nxpi); i number of elements
#' @param y.train - n-vector of class labels (must be a factor)
#' @param alphaList = list of alpha values
#' @param lambdaList = list of lambda values
#' @param family - can be "binomial" or "multinomial"
#' @param X.testList - list of test datasets (nxpi); i number of elements
#' @param y.test - n-vector of class labels (must be a factor)
#' @param filter - pre-filtering of initial datasets - "none" or "p.value"
#' @param topranked - Number of topranked features based on differential expression to use to build classifer
#' @param keepVarList - which variables to keep and not omit (set to NULL if no variables are forced to be kept)
#' @return model
#' @return testPerf
#' @return X.trainList
#' @return y.train
#' @return alphaList
#' @return lambdaList
#' @return family
#' @return X.testList
#' @return y.test
#' @return filter
#' @return topranked
#' @return keepVarList
#' @export
ensembleEnet = function(X.trainList, y.train, alphaList, lambdaList, family = "binomial",
  X.testList=NULL, y.test=NULL, filter="none", topranked=50, keepVarList=NULL){
  if (class(y.train) == "character")
    stop("y.train is not a factor")
  ## load libraries
  library(glmnet); library(limma); library(pROC);
  library(OptimalCutpoints); library(tidyverse);

  ## perform pre-filtering (none, p-value, and keep certain variables)
  if (filter == "none") {
    X1.trainList <- X.trainList
  }
  if (filter == "p.value") {
    X1.trainList <- lapply(X.trainList, function(i){
      design <- model.matrix(~y.train)
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
      fit <- glmnet(X, y.train, family = "binomial", alpha = alpha,
        penalty.factor = penalty.factor)
      cv.fit <- cv.glmnet(X, y.train, family = "binomial")
    } else {
      fit <- glmnet(X, y.train, family = "multinomial", alpha = alpha,
        type.multinomial = "grouped", penalty.factor = penalty.factor)
      cv.fit <- cv.glmnet(X, y.train, family = "multinomial")
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
    if(!any(levels(y.train) == levels(y.test)))
      stop("levels of y.train and y.test are not in the same order")

    testPerf <- mapply(function(mod, test){
      predictResponse <- unlist(predict(mod$fit, newx = test[, rownames(mod$Coefficients)[-1]], s = mod$lambda, type = "class"))
      probs <- predict(mod$fit, newx = test[, rownames(mod$Coefficients)[-1]], s = mod$lambda, type = "response") %>% as.numeric
      names(probs) <- rownames(predictResponse)
      predictResponse <- as.character(predictResponse)
      names(predictResponse) <- names(probs)

      ## compute error rate
      mat <- table(pred=factor(as.character(predictResponse), levels = levels(y.train)), truth=y.test)
      mat2 <- mat
      diag(mat2) <- 0
      classError <- colSums(mat2)/colSums(mat)
      er <- sum(mat2)/sum(mat)
      ber <- mean(classError)
      error <- c(classError, er, ber) %>% matrix
      rownames(error) <- c(names(classError), "Overall", "Balanced")
      colnames(error) <- "Error_0.5"

      ## compute AUROC
      if(length(y.test) > 1) {
        if(nlevels(y.train) == 2){
          y.test <- factor(as.character(y.test), levels(y.train))
          perfTest <- amritr::tperformance(weights = as.numeric(as.matrix(probs)), trueLabels = y.test) %>% as.matrix
          colnames(perfTest) <- paste(levels(y.train), collapse = "_vs_")
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
    y.train=y.train, alphaList=alphaList, lambdaList=lambdaList, family=family, X.testList=X.testList,
    y.test=y.test, filter=filter, topranked=topranked, keepVarList=keepVarList))
}



#' Estimate classification performance using repeated cross-validation using an elastic net classifier
#'
#'
#' @param object - ensembleEnet object
#' @param validation = "Mfold" or "loo"
#' @param M - # of folds
#' @param iter - Number of iterations of cross-validation
#' @param threads - # of cores, running each iteration on a separate node
#' @param progressBar = TRUE (show progress bar or not)
#' @export
perf.ensembleEnet = function(object, validation = "Mfold", M = 5, iter = 5, threads = 5, progressBar = TRUE){
  library(tidyverse)
  X.trainList=object$X.trainList
  y.train=object$y.train
  alphaList=object$alphaList
  lambdaList=object$lambdaList
  family=object$family
  filter=object$filter
  topranked=object$topranked
  keepVarList=object$keepVarList

  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) createFolds(y.train, k = M))
    require(parallel)
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", threads))
    parallel::clusterExport(cl, varlist = c("ensembleEnetCV",
      "ensembleEnet", "X.trainList", "y.train", "alphaList", "lambdaList",
      "family", "filter", "topranked", "keepVarList",
      "M", "folds", "progressBar"),
      envir = environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi, X.trainList, y.train, alphaList, lambdaList, family, filter, topranked, keepVarList, M, progressBar) {
      ensembleEnetCV(X.trainList=X.trainList, y.train=y.train, alphaList=alphaList, lambdaList=lambdaList, family=family, filter=filter, topranked=topranked, keepVarList=keepVarList, M=M, folds=foldsi, progressBar=progressBar)
    }, X.trainList, y.train, alphaList, lambdaList, family, filter, topranked,
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
    n <- length(y.train)
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- ensembleEnetCV(X.trainList, y.train, alphaList, lambdaList, family, filter, topranked,
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


#' Estimate classification performance using cross-validation using an elastic net classifier
#'
#' @param X.trainList list of training datasets (nxpi); i number of elements
#' @param y.train n-vector of class labels (must be a factor)
#' @param alphaList list of alpha values
#' @param lambdaList list of lambda values
#' @param family can be "binomial" or "multinomial"
#' @param filter pre-filtering of initial datasets - "none" or "p.value"
#' @param topranked Number of topranked features based on differential expression to use to build classifer
#' @param keepVarList which variables to keep and not omit (set to NULL if no variables are forced to be kept)
#' @param M # of folds
#' @param folds list of length M, where each element contains the indices for samples for a given fold
#' @param progressBar (TRUE/FALSE) - show progress bar or not
#' @return error computes error rate (each group, overall and balanced error rate)
#' @return perfTest classification performance measures
#' @export
ensembleEnetCV = function(X.trainList, y.train, alphaList, lambdaList, family="binomial", filter="none", topranked=50,
  keepVarList=NULL, M=5, folds=5, progressBar=FALSE){
  J <- length(X.trainList)
  assign("X.training", NULL, pos = 1)
  assign("y.training", NULL, pos = 1)
  X.training = lapply(folds, function(x) {
    lapply(1:J, function(y) {
      X.trainList[[y]][-x, , drop = FALSE]
    })
  })
  y.training = lapply(folds, function(x) {
    y.train[-x]
  })
  X.test = lapply(folds, function(x) {
    lapply(1:J, function(y) {
      X.trainList[[y]][x, , drop = FALSE]
    })
  })
  y.test = lapply(folds, function(x) {
    y.train[x]
  })
  avgProbList <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    ## build ensemble panel
    result <- ensembleEnet(X.trainList=X.training[[i]], y.train=y.training[[i]],
      alphaList, lambdaList, family = family,
      X.testList=X.test[[i]], y.test=y.test[[i]], filter, topranked, keepVarList)
    # combine predictions using average probability
    avgProbList[[i]] <- do.call(cbind, lapply(result$testPerf,
      function(i) {
        i$probs
      })) %>% rowMeans
  }

  probs <- unlist(avgProbList)
  ## Error and AUROC
  predictResponse <- rep(levels(y.train)[1], length(probs))
  predictResponse[probs >= 0.5] <- levels(y.train)[2]

  ## compute error rate
  truth <- sapply(strsplit(names(unlist(y.test)), "\\."), function(i) i[2])
  if(!all(names(probs) == truth))
    stop("predicted probability is not in the same order as the test labels")
  mat <- table(pred=factor(as.character(predictResponse), levels = levels(y.train)), truth=unlist(y.test))
  mat2 <- mat
  diag(mat2) <- 0
  classError <- colSums(mat2)/colSums(mat)
  er <- sum(mat2)/sum(mat)
  ber <- mean(classError)
  error <- c(classError, er, ber) %>% matrix
  rownames(error) <- c(names(classError), "Overall", "Balanced")
  colnames(error) <- "Error_0.5"

  ## compute AUROC
  if(nlevels(y.train) == 2){
    perfTest <- amritr::tperformance(weights = probs, trueLabels = unlist(y.test)) %>% as.matrix
    colnames(perfTest) <- "perf"
  } else {
    perfTest <- NA
  }

  return(list(error = error, perfTest = perfTest))
}
