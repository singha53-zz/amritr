library(glmnet)

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

## cross-validation Enet
perf.enet = function (object, validation = c("Mfold", "loo"), folds = 10,
  progressBar = TRUE, near.zero.var = FALSE)
{
  X = object$X
  Y = object$Y
  n = nrow(X)
  alpha = object$alpha
  nzv = nearZeroVar(X)
  if (length(nzv$Position > 0)) {
    warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
    X = X[, -nzv$Position, drop = TRUE]
    if (ncol(X) == 0) {
      stop("No more predictors after Near Zero Var has been applied!")
    }
  }

  if (validation == "Mfold") {
    if (is.list(folds)) {
      if (length(folds) < 2 | length(folds) > n)
        stop("Invalid number of folds.")
      if (length(unique(unlist(folds))) != n)
        stop("Invalid folds.")
      M = length(folds)
    }
    else {
      if (is.null(folds) || !is.numeric(folds) || folds <
          2 || folds > n)
        stop("Invalid number of folds.")
      else {
        M = round(folds)
        folds = split(sample(1:n), rep(1:M, length = n))
      }
    }
  }
  else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
  }

  probs <- list()
  stop.user = FALSE
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
    enet.res = suppressWarnings(enet(X.train, Y.train, alpha = alpha, lambda = NULL))
    probs[[i]] <- predict(enet.res$fit, newx=X.test, s = enet.res$lambda, type='response')
  }
  probs <- unlist(probs)
  if (stop.user == TRUE & validation == "Mfold")
    stop("The folds value was set too high to perform cross validation. Choose validation = \"loo\" or set folds to a lower value")
  if (progressBar == TRUE)
    cat("\n")

  result = list()
  result$nzvX = nzv$Position
  result$probs = probs
  result$folds = folds
  method = "enet.mthd"
  result$meth = "enet.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}
