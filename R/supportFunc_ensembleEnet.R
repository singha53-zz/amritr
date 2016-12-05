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
ensembleEnet = function(X.train, Y.train, alphaList, lambdaList, X.test, Y.test, filter, topranked){
  if(is.null(X.test)){  ## run if no test set is provided (only build ensemble models)
    result <- mapply(function(X.train, alpha, lambda) {
      enet(X = X.train, Y = Y.train, alpha = alpha, lambda = lambda, family = "multinomial", X.test = NULL,
        Y.test = NULL, filter = filter, topranked = topranked)
    }, X.train = X.train, alpha = alphaList, lambda = lambdaList, SIMPLIFY = FALSE)
    Y.vote <- error <- NA
  } else { ## run if test set is provided
    if(length(X.train) == length(X.test)){   ## if the same numbers of train and test datasets are available
      result <- mapply(function(X.train, X.test, alpha, lambda) {
        enet(X = X.train, Y = Y.train, alpha = alpha, lambda = lambda, family = "multinomial", X.test = X.test,
          Y.test = Y.test, filter = filter, topranked = topranked)
      }, X.train = X.train, X.test = X.test, alpha = alphaList, lambda = lambdaList, SIMPLIFY = FALSE)
      predConcat <- do.call(cbind, lapply(result, function(i) {
        i$predictResponse
      }))
      Y.vote <- apply(predConcat, 1, function(z) {
        temp = table(z)
        if (length(names(temp)[temp == max(temp)]) > 1) {
          "zz"
        }
        else {
          names(temp)[temp == max(temp)]
        }
      })
      temp = table(pred = factor(Y.vote, levels = c(levels(Y.test), "zz")), truth = unlist(Y.test))
      diag(temp) <- 0
      error = c(colSums(temp)/summary(Y.test), sum(temp)/length(Y.test),
        mean(colSums(temp)/summary(Y.test)))
      names(error) <- c(names(error)[1:nlevels(Y.test)], "ER", "BER")
    } else { ## if the different numbers of train and test datasets are available
      comDataset <- intersect(names(X.train), names(X.test))
      X.train <- X.train[comDataset]
      X.test <- X.test[comDataset]
      alphaList <- alphaList[comDataset]
      lambdaList <- lambdaList[comDataset]

      result <- mapply(function(X.train, X.test, alpha, lambda) {
        enet(X = X.train, Y = Y.train, alpha = alpha, lambda = lambda, family = "multinomial", X.test = X.test,
          Y.test = Y.test, filter = filter, topranked = topranked)
      }, X.train = X.train, X.test = X.test, alpha = alphaList, lambda = lambdaList, SIMPLIFY = FALSE)
      predConcat <- do.call(cbind, lapply(result, function(i) {
        i$predictResponse
      }))
      Y.vote <- apply(predConcat, 1, function(z) {
        temp = table(z)
        if (length(names(temp)[temp == max(temp)]) > 1) { "zz"
        }
        else {
          names(temp)[temp == max(temp)]
        }
      })
      temp = table(pred = factor(Y.vote, levels = c(levels(Y.test), "zz")), truth = unlist(Y.test))
      diag(temp) <- 0
      error = c(colSums(temp)/summary(Y.test), sum(temp)/length(Y.test),
        mean(colSums(temp)/summary(Y.test)))
      names(error) <- c(names(error)[1:nlevels(Y.test)], "ER", "BER")
    }
  }

  return(list(result = result, Y.vote = Y.vote, perfTest = error, X.train = X.train, Y.train = Y.train, alphaList = alphaList, lambdaList = lambdaList, filter = filter, topranked = topranked))
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
perfEnsemble = function(object, validation = "Mfold", M = M, iter = iter, threads = threads,
  progressBar = TRUE){
  library(dplyr)
  X <- object$X.train
  Y = object$Y.train
  n = length(Y)
  alpha = object$alphaList
  lambda = object$lambdaList
  filter = object$filter
  topranked = object$topranked
  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) caret::createFolds(Y, k = M))
    require(parallel)
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", threads))
    parallel::clusterExport(cl, varlist = c("ensembleCV",
      "ensembleEnet", "enet", "X", "Y", "alpha", "lambda",
      "M", "folds", "progressBar", "filter", "topranked"), envir = environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi,
      X, Y, alpha, lambda, M, progressBar, filter, topranked) {
      ensembleCV(X = X, Y = Y, alpha = alpha, lambda = lambda,
        M = M, folds = foldsi, progressBar = progressBar, filter = filter, topranked = topranked)
    }, X, Y, alpha, lambda, M, progressBar, filter, topranked) %>% amritr::zip_nPure()
    parallel::stopCluster(cl)
    perf <- do.call(rbind, cv$error) %>% as.data.frame %>%
      tidyr::gather(ErrName, Err) %>%
      dplyr::group_by(ErrName) %>%
      dplyr::summarise(Mean = mean(Err), SD = sd(Err))
  }
  else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- ensembleCV(X, Y, alpha, lambda, M, folds, progressBar, filter, topranked)
    perf <- data.frame(Mean = cv$error) %>% mutate(ErrName = rownames(.))
    perf$SD <- NA
  }
  result = list()
  result$perf = perf
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
ensembleCV = function(X, Y, alpha, lambda, M, folds, progressBar, filter, topranked){
  J <- length(X)
  assign("X.training", NULL, pos = 1)
  assign("Y.training", NULL, pos = 1)
  X.training = lapply(folds, function(x) {
    lapply(1:J, function(y) {
      X[[y]][-x, ]
    })
  })
  Y.training = lapply(folds, function(x) {
    Y[-x]
  })
  X.test = lapply(folds, function(x) {
    lapply(1:J, function(y) {
      X[[y]][x, , drop = FALSE]
    })
  })
  Y.test = lapply(folds, function(x) {
    Y[x]
  })
  predConcatList <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    result <- mapply(function(X.train, X.test, alpha, lambda) {
      enet(X.train, Y.training[[i]], alpha = alpha, lambda = lambda,
        X.test = X.test, Y.test = Y.test[[i]], family = "multinomial", filter = filter, topranked = topranked)
    }, X.train = X.training[[i]], X.test = X.test[[i]], alpha = alpha, lambda = lambda, SIMPLIFY = FALSE)
    predConcatList[[i]] <- do.call(cbind, lapply(result,
      function(i) {
        i$predictResponse
      }))
  }
  predConcat <- do.call(rbind, predConcatList)
  Y.vote <- apply(predConcat, 1, function(z) {
    temp = table(z)
    if (length(names(temp)[temp == max(temp)]) > 1) {
      "zz"
    }
    else {
      names(temp)[temp == max(temp)]
    }
  })
  temp = table(pred = factor(Y.vote, levels = c(levels(Y),
    "zz")), truth = unlist(Y)[unlist(folds)])
  diag(temp) <- 0
  error = c(colSums(temp)/summary(Y), sum(temp)/length(Y),
    mean(colSums(temp)/summary(Y)))
  names(error) <- c(names(error)[1:nlevels(Y)], "ER", "BER")
  return(list(error = error))
}
