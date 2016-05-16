#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
## Elastic net
enet = function(X, Y, alpha, lambda=NULL){
  library(glmnet)

  # Run Elastic net classifier
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
  return(list(X = X, Y = Y, fit = fit, enet.panel = enet.panel, lambda = lambda, alpha = alpha))
}

#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
runCV = function(X, Y, alpha, M, folds, progressBar){
  probs <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    omit = folds[[i]]
    X.train = X[-omit, ]
    Y.train = Y[-omit]
    X.test = matrix(X[omit, ], nrow = length(omit))
    enet.res = suppressWarnings(enet(X.train, Y.train, alpha = alpha, lambda = NULL))
    probs[[i]] <- predict(enet.res$fit, newx=X.test, s = enet.res$lambda, type='response')
  }
  probs <- unlist(probs)
  trueLabels = Y[unlist(folds)]
  library(pROC); library(OptimalCutpoints);
  perf <- amritr::tperformance(weights = probs, trueLabels = trueLabels, direction = "auto")
  return(list(probs=probs, trueLabels=trueLabels, perf=perf))
}

#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
perf.enet = function (object, validation = c("Mfold", "loo"), M = 5, iter = 10, threads = 4, progressBar = TRUE){
  require(caret, quietly=TRUE)   ## require to make cv folds
  require(snow, quietly=TRUE)   ## require to parallelize
  library(dplyr)
  X = object$X
  Y = object$Y
  n = nrow(X)
  alpha = object$alpha

  if (validation == "Mfold") {
   folds <- lapply(1 : iter, function(i) createFolds(1:n, k = M))

   # parallelize each iteration
   cl <- makeCluster(threads, type = "SOCK")
   clusterExport(cl, c("runCV", "enet"))
   clusterExport(cl, c("X", "Y", "alpha", "M", "folds", "progressBar"))
   cv <- parLapply(cl, folds, function(x){
     runCV(X, Y, alpha, M, folds=x, progressBar)
   }) %>% amritr::zip_nPure()
   stopCluster(cl)
  }
  else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- runCV(X, Y, alpha, M, folds, progressBar)
  }

  result = list()
  result$folds = folds
  result$probs = cv$probs
  result$trueLabels = cv$trueLabels
  result$perf = cv$perf
  method = "enet.mthd"
  result$meth = "enet.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}



#library(mixOmics)
#data(srbct)
#X <- srbct$gene[srbct$class %in% c("EWS", "RMS"),]
#Y <- droplevels(srbct$class[srbct$class %in% c("EWS", "RMS")])

#model <- enet(X, Y, alpha=0.5, lambda=NULL)
#cv <- perf.enet(model, validation = "Mfold", M = 5, iter = 10, progressBar = TRUE)

