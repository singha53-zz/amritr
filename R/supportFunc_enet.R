#' Elastic net classification panel
#'
#' build an elastic net classification panel
#' @param X nxp matrix - training dataset
#' @param Y categorical variables
#' @param alpha = 1 (lasso), alpha = 0 (ridge), 0 < alpha < 1 (elastic net penalty)
#' @param lambda = strength of elastic net penalty
#' @param family "binomial" or "multinomial"
#' @param X.test nxp matrx  - test dataset
#' @param filter = "none" or "p.value"
#' @param topranked = 50 (top number of features to select and build a classifier)
#' @export
enet = function (X, Y, alpha, lambda = NULL, family, X.test = NULL,
  Y.test = NULL, filter = "p.value", topranked = 50)
{
  library(glmnet); library(limma); library(pROC); library(OptimalCutpoints);

  if (filter == "none") {
    X1 <- X
  }
  if (filter == "p.value") {
    design <- model.matrix(~Y)
    fit <- eBayes(lmFit(t(X), design))
    top <- topTable(fit, coef = 2, adjust.method = "BH",
      n = nrow(fit))
    X1 <- X[, rownames(top)[1:topranked]]
  }

  if (family == "binomial") {
    fit <- glmnet(X1, Y, family = "binomial", alpha = alpha)
    cv.fit <- cv.glmnet(X1, Y, family = "binomial")
    if (is.null(lambda)) {
      lambda = cv.fit$lambda.min
    }
    else {
      lambda = lambda
    }
    Coefficients <- coef(fit, s = lambda)
    Active.Index <- which(Coefficients[, 1] != 0)
    Active.Coefficients <- Coefficients[Active.Index, ]
    enet.panel <- names(Active.Coefficients)[-1]
    enet.panel.length <- length(enet.panel)
  }
  if (family == "multinomial") {
    fit <- glmnet(X1, Y, family = "multinomial", alpha = alpha,
      type.multinomial = "grouped")
    cv.fit <- cv.glmnet(X1, Y, family = "multinomial")
    if (is.null(lambda)) {
      lambda = cv.fit$lambda.min
    }
    else {
      lambda = lambda
    }
    Coefficients <- coef(fit, s = lambda)
    Active.Index <- which(Coefficients[[1]][, 1] != 0)
    Active.Coefficients <- Coefficients[[1]][Active.Index,
      ]
    enet.panel <- names(Active.Coefficients)[-1]
    enet.panel.length <- length(enet.panel)
  }
  if (!is.null(X.test)) {
    library(pROC)
    library(OptimalCutpoints)
    probs <- predict(fit, newx = X.test, s = lambda, type = "response")
    predictResponse <- unlist(predict(fit, newx = X.test,
      s = lambda, type = "class"))
    if (family == "binomial") {
      perfTest <- amritr::tperformance(weights = as.numeric(as.matrix(probs)),
        trueLabels = Y.test)
    }
    else {
      mat <- table(factor(as.character(predictResponse),
        levels = levels(Y.test)), Y.test)
      mat2 <- mat
      diag(mat2) <- 0
      classError <- colSums(mat2)/colSums(mat)
      er <- sum(mat2)/sum(mat)
      ber <- mean(classError)
      perfTest <- c(classError, er, ber)
      names(perfTest) <- c(names(classError), "ER", "BER")
    }
  }
  else {
    perfTest <- predictResponse <- probs <- NA
  }
  return(list(X = X, Y = Y, fit = fit, enet.panel = enet.panel,
    lambda = lambda, alpha = alpha, family = family, probs = probs,
    Active.Coefficients = Active.Coefficients, perfTest = perfTest,
    predictResponse = predictResponse, filter = filter, topranked = topranked))
}

#' interal function (enet cross-validation)
#'
#' Estimate test error of elastic net panel
#' @param X nxp matrix - training dataset
#' @param Y categorical variables
#' @param M = # of folds
#' @param folds list of length M specifiying which samples are in which fold
#' @param family "binomial" or "multinomial"
#' @param filter = "none" or "p.value"
#' @param topranked = 50 (top number of features to select and build a classifier)
#' @export
enetCV = function (X, Y, alpha, lambda, M, folds, progressBar, family,
  filter, topranked)
{
  library(pROC); library(limma); library(glmnet); library(OptimalCutpoints)

  probs <- predictResponseList <- enet.panel <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    omit = folds[[i]]
    X.train = X[-omit, ]
    Y.train = Y[-omit]

    if (filter == "none") {
      X.train1 <- X.train
    }
    if (filter == "p.value") {
      design <- model.matrix(~Y.train)
      fit <- eBayes(lmFit(t(X.train), design))
      top <- topTable(fit, coef = 2, adjust.method = "BH", n = nrow(fit))
      X.train1 <- X.train[, rownames(top)[1:topranked]]
    }
    X.test1 = X[omit, colnames(X.train1), drop = FALSE]

    if (family == "binomial") {
      fit <- glmnet(X.train1, Y.train, family = "binomial", alpha = alpha)
      cv.fit <- cv.glmnet(X.train1, Y.train, family = "binomial")
      if (is.null(lambda)) {
        lambda = cv.fit$lambda.min
      }
      else {
        lambda = lambda
      }
      Coefficients <- coef(fit, s = lambda)
      Active.Index <- which(Coefficients[, 1] != 0)
      Active.Coefficients <- Coefficients[Active.Index, ]
      enet.panel[[i]] <- names(Active.Coefficients)[-1]
    }
    if (family == "multinomial") {
      fit <- glmnet(X.train1, Y.train, family = "multinomial", alpha = alpha,
        type.multinomial = "grouped")
      cv.fit <- cv.glmnet(X.train1, Y.train, family = "multinomial")
      if (is.null(lambda)) {
        lambda = cv.fit$lambda.min
      }
      else {
        lambda = lambda
      }
      Coefficients <- coef(fit, s = lambda)
      Active.Index <- which(Coefficients[[1]][, 1] != 0)
      Active.Coefficients <- Coefficients[[1]][Active.Index,
        ]
      enet.panel[[i]] <- names(Active.Coefficients)[-1]
    }

    probs[[i]] <- predict(fit, newx = X.test1, s = lambda, type = "response")
    predictResponseList[[i]] <- predict(fit, newx = X.test1, s = lambda, type = "class")
  }
  predictResponse <- unlist(predictResponseList)
  if (family == "binomial") {
    probs <- unlist(probs)
    trueLabels = Y[unlist(folds)]
    library(pROC)
    library(OptimalCutpoints)
    perf <- amritr::tperformance(weights = probs, trueLabels = trueLabels)
  }
  else {
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
  return(list(probs = probs, trueLabels = trueLabels, perf = perf,
    enet.panel = enet.panel, predictResponse = predictResponse))
}

#' cross-validation function for elastic net panel
#'
#' Estimate test error of elastic net panel
#' @param X nxp matrix - training dataset
#' @param Y categorical variables
#' @param alpha = 1 (lasso), alpha = 0 (ridge), 0 < alpha < 1 (elastic net penalty)
#' @param family "binomial" or "multinomial"
#' @param lambda = strength of elastic net penalty
#' @param M = # of folds
#' @param iter = number of times to repeat cross-valiation
#' @param threads - number of cpus (each cross-valiation scheme performed on a separate node)
#' @param progressBar - show progressbar (TRUE/FALE)
#' @export
perf.enet = function (object, validation = c("Mfold", "loo"), M = 5, iter = 10,
  threads = 4, progressBar = TRUE)
{
  library(dplyr)
  library(tidyr)
  X = object$X
  Y = object$Y
  n = nrow(X)
  alpha = object$alpha
  family = object$family
  lambda = object$lambda
  filter = object$filter
  topranked = object$topranked
  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) caret::createFolds(Y,
      k = M))
    require(parallel)
    cl <- parallel::makeCluster(mc <- getOption("cl.cores",
      threads))
    parallel::clusterExport(cl, varlist = c("enetCV", "enet",
      "X", "Y", "alpha", "lambda", "M", "folds", "progressBar",
      "family", "filter", "topranked"), envir = environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi,
      X, Y, alpha, lambda, M, progressBar, family, filter,
      topranked) {
      enetCV(X = X, Y = Y, alpha = alpha, lambda = lambda,
        M = M, folds = foldsi, progressBar = progressBar,
        family = family, filter = filter, topranked = topranked)
    }, X, Y, alpha, lambda, M, progressBar, family, filter,
      topranked) %>% amritr::zip_nPure()
    parallel::stopCluster(cl)
    perf <- do.call(rbind, cv$perf) %>% as.data.frame %>%
      gather(ErrName, Err) %>% dplyr::group_by(ErrName) %>%
      dplyr::summarise(Mean = mean(Err), SD = sd(Err))
  }
  else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- enetCV(X, Y, alpha, lambda, M, folds, progressBar,
      family, filter, topranked)
    perf <- cv$perf
  }
  result = list()
  result$folds = folds
  result$probs = cv$probs
  result$trueLabels = cv$trueLabels
  result$panels = cv$enet.panel
  result$perf = perf
  method = "enet.mthd"
  result$meth = "enet.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}
