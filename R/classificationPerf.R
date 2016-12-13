#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
## determine AUC from predictions and true labels
tperformance = function(weights, trueLabels){
  ## Determine optimal cut-off values and associated performance measures
  df = data.frame(prob = as.numeric(weights),
    status = model.matrix(~factor(as.character(trueLabels), levels = levels(trueLabels)))[, 2])
  roc.score = roc(response = df$status, predictor = weights, plot = FALSE, percent = TRUE, na.rm = TRUE, direction = "<")

  optimal.cutpoint.Youden <- optimal.cutpoints(X = "prob", status = "status", tag.healthy = 0, methods = "Youden",
    data = df, control = control.cutpoints(), ci.fit = FALSE,
    conf.level = 0.95, trace = FALSE, pop.prev = 0.5)
  optimalValues <- round(c(summary(optimal.cutpoint.Youden)$p.table$Global$Youden[[1]][1:5, ], roc.score$auc/100), 3)
  names(optimalValues) <- c(names(optimalValues)[-length(names(optimalValues))], "AUC")
  optimalValues
}

#' createFolds function copied from the caret R-package
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param y Y vector of phenotypic lables
#' @param k number of iterations of the cross-valiation
#' @export
createFolds = function (y, k = 10, list = TRUE, returnTrain = FALSE)
{
  if (class(y)[1] == "Surv")
    y <- y[, "time"]
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2)
      cuts <- 2
    if (cuts > 5)
      cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0)
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k,
          size = numInClass[i])
      }
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
      sep = "")
    if (returnTrain)
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}
