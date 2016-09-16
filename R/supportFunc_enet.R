#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
## Elastic net
enet = function(X, Y, alpha, lambda=NULL, family, X.test=NULL, Y.test=NULL){
  library(glmnet)

  # Run Elastic net classifier
  if(family == "binomial"){
    fit <- glmnet(X, Y, family="binomial", alpha = alpha)
    cv.fit <- cv.glmnet(X, Y, family="binomial")
    if(is.null(lambda)){
      set.seed(1)
      lambda = cv.fit$lambda.min
    } else {
      lambda = lambda
    }
    Coefficients <- coef(fit, s = lambda)
    Active.Index <- which(Coefficients[, 1] != 0)
    Active.Coefficients  <- Coefficients[Active.Index,]
    enet.panel <- names(Active.Coefficients)[-1]
    enet.panel.length <- length(enet.panel)
  } else {
    fit <- glmnet(X, Y, family="multinomial", alpha = alpha, type.multinomial = "grouped")
    cv.fit <- cv.glmnet(X, Y, family="multinomial")
    if(is.null(lambda)){
      set.seed(1)
      lambda = cv.fit$lambda.min
    } else {
      lambda = lambda
    }
    Coefficients <- coef(fit, s = lambda)
    Active.Index <- which(Coefficients[[1]][, 1] != 0)
    Active.Coefficients  <- Coefficients[[1]][Active.Index,]
    enet.panel <- names(Active.Coefficients)[-1]
    enet.panel.length <- length(enet.panel)
  }

  ## assess panel performance using test data if available
  if(!is.null(X.test)){
    library(pROC)
    library(OptimalCutpoints)
    probs <- predict(fit, newx = X.test, s = lambda, type = "response")
    predictResponse <- unlist(predict(fit, newx = X.test, s = lambda, type = "class"))

    if(family == "binomial"){
      perfTest <- amritr::tperformance(weights = probs, trueLabels = Y.test, direction = "auto")
    } else {
      mat <- table(Y.test, predictResponse)
      mat2 <- mat
      diag(mat2) <- 0
      classError <- colSums(mat2)/colSums(mat)
      er <- sum(mat2)/sum(mat)
      ber <- mean(classError)
      perfTest <- c(classError, er, ber)
      names(perfTest) <- c(names(classError), "ER", "BER")
    }
  } else {
    perfTest <- predictResponse <- NA
  }

  return(list(X = X, Y = Y, fit = fit, enet.panel = enet.panel, lambda = lambda, alpha = alpha, family = family, perfTest=perfTest, predictResponse=predictResponse))
}


#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
runCV = function(X, Y, alpha, M, folds, progressBar, family){
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
    enet.res = suppressWarnings(enet(X.train, Y.train, alpha = alpha, lambda = NULL, family = family))
    probs[[i]] <- predict(enet.res$fit, newx=X.test, s = enet.res$lambda, type='response')
    predictResponseList[[i]] <- predict(enet.res$fit, newx=X.test, s = enet.res$lambda, type='class')
  }
  predictResponse <- unlist(predictResponseList)

  ## Error rate, AUC
  if(family == "binomial"){
    probs <- unlist(probs)
    trueLabels = Y[unlist(folds)]
    library(pROC); library(OptimalCutpoints);
    perf <- amritr::tperformance(weights = probs, trueLabels = trueLabels, direction = "auto")

  } else {
    ## Error rate
    trueLabels = Y[unlist(folds)]
    mat <- table(trueLabels, predictResponse)
    mat2 <- mat
    diag(mat2) <- 0
    classError <- colSums(mat2)/colSums(mat)
    er <- sum(mat2)/sum(mat)
    ber <- mean(classError)
    perf <- c(classError, er, ber)
    names(perf) <- c(names(classError), "ER", "BER")
  }

  return(list(probs=probs, trueLabels=trueLabels, perf=perf, predictResponse=predictResponse))
}

#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
perf.enet = function (object, validation = c("Mfold", "loo"), M = 5, iter = 10,
  threads = 4, progressBar = TRUE)
{
  require(caret, quietly = TRUE)
  library(dplyr)
  X = object$X
  Y = object$Y
  n = nrow(X)
  alpha = object$alpha
  family = object$family
  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) createFolds(Y, k = M))
    require(parallel)
    cl <- parallel::makeCluster(mc <- getOption("cl.cores",
      threads))
    parallel::clusterExport(cl, varlist = c("runCV", "enet",
      "X", "Y", "alpha", "M", "folds", "progressBar", "family"),
      envir = environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi,
      X, Y, alpha, M, progressBar, family) {
      runCV(X = X, Y = Y, alpha = alpha, M = M, folds = foldsi,
        progressBar = progressBar, family = family)
    }, X, Y, alpha, M, progressBar, family) %>% amritr::zip_nPure()
    parallel::stopCluster(cl)

    perf <- do.call(rbind, cv$perf) %>% as.data.frame %>% gather(ErrName,
      Err) %>% dplyr::group_by(ErrName) %>% dplyr::summarise(Mean = mean(Err),
        SD = sd(Err))
  }
  else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- runCV(X, Y, alpha, M, folds, progressBar, family = family)
    perf <- cv$perf
  }
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
