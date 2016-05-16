#' @title Prediction fuction for splsda
#'
#' @description predict scores for newdata
#' @param object splsda model
#' @param method.predict type of distance measure to use: "max.dist", "centroids.dist", "mahalanobis.dist" or "all"
#' @param validation Either an M-fold cross-valiadtion where M is an integr between 1-n (total number of samples)
#' @param folds Number of folds, used only if validation = "Mfold"
#' @param progressBar display progress bar (TRUE) or Not (FALSE)
#' @import mixOmics
#' @seealso NULL
#' @return NULL
#' @aliases NULL
#' @examples \dontrun{
#'
#'}
## Predict new data with splsda model
perf.splsda2 = function(object, method.predict = c("all", "max.dist", "centroids.dist",
  "mahalanobis.dist"), validation = c("Mfold", "loo"), folds = 10,
  progressBar = TRUE, near.zero.var = FALSE){
  X = object$X
  level.Y = object$names$Y
  Y = object$ind.mat
  Y = mixOmics::map(Y)
  Y = factor(Y, labels = level.Y)
  ncomp = object$ncomp
  n = nrow(X)
  keepX = object$keepX
  tol = object$tol
  max.iter = object$max.iter
  features <- list()
  for (k in 1:ncomp) {
    features[[k]] = NA
  }
  method.predict = match.arg(method.predict, choices = c("all",
    "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
  if (any(method.predict == "all"))
    nmthdd = 3
  else nmthdd = length(method.predict)
  nzv = nearZeroVar(X)
  if (length(nzv$Position > 0)) {
    warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
    X = X[, -nzv$Position, drop = TRUE]
    if (ncol(X) == 0) {
      stop("No more predictors after Near Zero Var has been applied!")
    }
  }
  error.fun = function(x, y) {
    error.vec = sweep(x, 1, y, FUN = "-")
    error.vec = (error.vec != 0)
    error.vec = apply(error.vec, 2, sum)/length(y)
    return(error.vec)
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
  error.mat = array(0, dim = c(ncomp, nmthdd, M))
  score <- list()
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
    spls.res = splsda(X.train, Y.train, ncomp, max.iter,
      tol, keepX = keepX, near.zero.var = near.zero.var)
    for (k in 1:ncomp) {
      features[[k]] = c(unlist(features[[k]]), selectVar(spls.res,
        comp = k)$name)
    }
    if (!is.null(spls.res$nzv$Position))
      X.test = X.test[, -spls.res$nzv$Position]
    Y.predict = predict(spls.res, X.test, method = method.predict)$class
    score[[i]] = predict(spls.res, X.test, method = method.predict)$predict
    error.mat[, , i] = sapply(Y.predict, error.fun, y = as.numeric(Y[omit]))
  }
  if (stop.user == TRUE & validation == "Mfold")
    stop("The folds value was set too high to perform cross validation. Choose validation = \"loo\" or set folds to a lower value")
  if (progressBar == TRUE)
    cat("\n")
  res = apply(error.mat, 1:2, mean)
  rownames(res) = paste("ncomp", 1:ncomp, sep = " ")
  colnames(res) = names(Y.predict)
  list.features = list()
  for (k in 1:ncomp) {
    remove.na = which(is.na(features[[k]]))
    list.features[[k]] = sort(table(as.factor(features[[k]][-remove.na]))/M,
      decreasing = TRUE)
  }
  names(list.features) = paste("comp", 1:ncomp)

  ## re-organize scores
  score <- lapply(1:ncomp, function(i){
    do.call(rbind, lapply(score, function(j){
      j[, , i]
    }))
  })

  result = list()
  result$error.rate = res
  result$features$stable = list.features
  result$nzvX = nzv$Position
  result$score = score
  method = "plsda.mthd"
  result$meth = "splsda.mthd"
  class(result) = c("perf", method)
  return(invisible(result))
}

