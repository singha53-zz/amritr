#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
ensembleSvm = function(X.train, Y.train, X.test, Y.test){
  result <- mapply(function(X.train, X.test){
    supportVectorMachine(X.train, Y.train, X.test=X.test, Y.test=Y.test)
  }, X.train = X.train, X.test = X.test, SIMPLIFY = FALSE)

  predConcat <- do.call(cbind, lapply(result, function(i){
    i$predictResponse
  }))

  Y.vote <- apply(predConcat, 1, function(z){
    temp = table(z)
    if (length(names(temp)[temp == max(temp)]) > 1){
      NA
    } else {
      names(temp)[temp == max(temp)]
    }
  })

  ## Compare error rates
  temp=table(pred=Y.vote, truth=Y.test)
  diag(temp) <- 0
  error = c(colSums(temp)/summary(Y.test), sum(temp)/length(Y.test), mean(colSums(temp)/summary(Y.test)))
  names(error) <- c(names(error)[1:nlevels(Y.test)], "ER", "BER")

  return(list(result=result, Y.vote=Y.vote, error=error, X.train=X.train, Y.train=Y.train, alpha=alpha))
}

#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
ensembleSvmCV = function(X, Y, M, folds, progressBar){
  J <- length(X)
  ### Start: Training samples (X.training and Y.training) and Test samples (X.test / Y.test)
  assign("X.training", NULL, pos = 1); assign("Y.training", NULL, pos = 1)
  X.training = lapply(folds, function(x){lapply(1:J, function(y) {X[[y]][-x, ]})})
  Y.training = lapply(folds, function(x) {Y[-x]});

  X.test = lapply(folds, function(x){lapply(1:J, function(y) {X[[y]][x, , drop = FALSE]})})
  Y.test = lapply(folds, function(x) {Y[x]});
  ### End: Training samples (X.training and Y.training) and Test samples (X.test / Y.test)

  predConcatList <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)

    result <- mapply(function(X.train, X.test){
      supportVectorMachine(X.train, Y.training[[i]], X.test=X.test, Y.test=Y.test[[i]])
    }, X.train = X.training[[i]], X.test = X.test[[i]], SIMPLIFY = FALSE)

    predConcatList[[i]] <- do.call(cbind, lapply(result, function(i){
      i$predictResponse
    }))
  }

  ## estimate error rate
  predConcat <- do.call(rbind, predConcatList)

  Y.vote <- apply(predConcat, 1, function(z){
    temp = table(z)
    if (length(names(temp)[temp == max(temp)]) > 1){
      NA
    } else {
      names(temp)[temp == max(temp)]
    }
  })

  ## Compare error rates
  temp=table(pred=Y.vote, truth=unlist(Y.test))
  diag(temp) <- 0
  error = c(colSums(temp)/summary(Y), sum(temp)/length(Y), mean(colSums(temp)/summary(Y)))
  names(error) <- c(names(error)[1:nlevels(Y)], "ER", "BER")

  return(list(error=error))
}

#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
perfEnsembleSvm = function(object, validation = "Mfold", M = M, iter = iter, threads = threads, progressBar = TRUE){

  library(dplyr)
  X <- object$X.train
  Y = object$Y.train
  n = length(Y)
  alpha = object$alpha
  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) caret::createFolds(Y, k = M))
    require(parallel)
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", threads))
    parallel::clusterExport(cl, varlist=c("ensembleSvmCV", "ensembleSvm", "supportVectorMachine", "X", "Y", "M", "folds", "progressBar"), envir=environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi, X, Y, M, progressBar){
      ensembleSvmCV(X=X, Y=Y, M=M, folds = foldsi, progressBar=progressBar)
    }, X, Y, M, progressBar) %>% amritr::zip_nPure()
    parallel::stopCluster(cl)

  } else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- ensembleSvmCV(X, Y, alpha, M, folds, progressBar)
  }

  ## Summarise performance results
  perf <- do.call(rbind, cv$error) %>% as.data.frame %>%
    gather(ErrName, Err) %>% dplyr::group_by(ErrName) %>%
    dplyr::summarise(Mean = mean(Err), SD = sd(Err))

  result = list()
  result$folds = folds
  result$probs = cv$probs
  result$trueLabels = cv$trueLabels
  result$perf = perf
  method = "enet.mthd"
  result$meth = "enet.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}
