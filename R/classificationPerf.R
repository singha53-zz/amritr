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
  roc.score = roc(response = trueLabels, predictor = weights, levels = levels(trueLabels), plot = TRUE, percent = TRUE, na.rm = TRUE, direction = "<")

  optimal.cutpoint.Youden <- optimal.cutpoints(X = "prob", status = "status", tag.healthy = 0, methods = "Youden",
    data = df, control = control.cutpoints(), ci.fit = FALSE,
    conf.level = 0.95, trace = FALSE, pop.prev = 0.5)
  optimalValues <- round(c(summary(optimal.cutpoint.Youden)$p.table$Global$Youden[[1]][1:5, ], roc.score$auc/100), 2)
  names(optimalValues) <- c(names(optimalValues)[-length(names(optimalValues))], "AUC")
  optimalValues
}
