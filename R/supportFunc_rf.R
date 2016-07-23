#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
rforest = function(X, Y, X.test=NULL, Y.test=NULL){
  library(randomForest)

  # Run random forest classifier
  colnames(X) <- paste("Feature", 1:ncol(X), sep=".")
  y.train0 <- Y[rownames(X)]
  fit = randomForest(y.train0 ~ ., data=X, importance=TRUE, proximity=TRUE)


  if(!is.null(X.test)){
    colnames(X.test) <- paste("Feature", 1:ncol(X), sep=".")
    ## assess panel performance using test data if available
    predictResponse <- predict(fit, as.matrix(X.test), type='response')

    mat <- table(Y.test, predictResponse)
    mat2 <- mat
    diag(mat2) <- 0
    classError <- colSums(mat2)/colSums(mat)
    er <- sum(mat2)/sum(mat)
    ber <- mean(classError)
    perfTest <- c(classError, er, ber)
    names(perfTest) <- c(names(classError), "ER", "BER")
  } else {
    perfTest <- predictResponse <- NA
  }

  return(list(X = X, Y = Y, fit = fit, perfTest=perfTest, predictResponse=predictResponse))
}


#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
rfCV = function(X, Y, M, folds, progressBar){
  probs <- predictResponseList <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    omit = folds[[i]]
    X.train = X[-omit, ]
    Y.train = Y[-omit]
    X.test = matrix(X[omit, ], nrow = length(omit))
    colnames(X.test) <- paste("Feature", 1:ncol(X.test), sep=".")
    rf.res = suppressWarnings(rforest(X=X.train, Y=Y.train))
    probs[[i]] <- predict(rf.res$fit, X.test, type='response')
    predictResponseList[[i]] <- predict(rf.res$fit, X.test, type='class')
  }
  predictResponse <- unlist(predictResponseList)

  ## Error rate
  trueLabels = Y[unlist(folds)]
  mat <- table(trueLabels, predictResponse)
  mat2 <- mat
  diag(mat2) <- 0
  classError <- colSums(mat2)/colSums(mat)
  er <- sum(mat2)/sum(mat)
  ber <- mean(classError)
  perf <- data.frame(Error = c(classError, er, ber))
  perf$ErrorRate <- c(names(classError), "ER", "BER")

  perf
}

#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
perf.rf = function (object, validation = c("Mfold", "loo"), M = 5, iter = 10,
  threads = 4, progressBar = TRUE)
{
  require(caret, quietly = TRUE)
  library(dplyr)
  X = object$X
  Y = object$Y
  n = nrow(X)
  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) createFolds(Y, k = M))
    require(parallel)
    cl <- parallel::makeCluster(mc <- getOption("cl.cores", threads))
    parallel::clusterExport(cl, varlist=c("rfCV", "rf", "X", "Y", "M", "folds", "progressBar"), envir=environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi, X, Y, M, progressBar){
      rfCV(X=X, Y=Y, M=M, folds = foldsi, progressBar=progressBar)
    }, X, Y, M, progressBar) %>% amritr::zip_nPure()
    parallel::stopCluster(cl)

  } else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- rfCV(X, Y, M, folds, progressBar)
  }

  ## Summarise performance results
  perf0 <- do.call(rbind, cv$Error)
  colnames(perf0) <- cv$ErrorRate[[1]]
  perf <- perf0 %>% as.data.frame() %>%
    gather(ErrName, Err) %>% dplyr::group_by(ErrName) %>%
    dplyr::summarise(Mean = mean(Err), SD = sd(Err))

  result = list()
  result$perf = perf
  method = "rf.mthd"
  result$meth = "rf.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}


