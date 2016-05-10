library(pROC)
library(OptimalCutpoints)

## determine AUC from predictions and true labels
tperformance = function(weights, trueLabels, direction){
  ## Determine optimal cut-off values and associated performance measures
  df = data.frame(prob = weights,
    status = model.matrix(~factor(as.character(trueLabels), levels = levels(trueLabels)))[, 2])
  roc.score = roc(response = as.character(trueLabels), predictor = weights, plot = TRUE, percent = TRUE, na.rm =TRUE, direction = direction)

  optimal.cutpoint.Youden <- optimal.cutpoints(X = "prob", status = "status", tag.healthy = 0, methods = "Youden",
    data = df, control = control.cutpoints(), ci.fit = FALSE,
    conf.level = 0.95, trace = FALSE, pop.prev = 0.5)
  optimalValues <- round(c(summary(optimal.cutpoint.Youden)$p.table$Global$Youden[[1]][1:5, ], roc.score$auc), 2)
  names(optimalValues) <- c(names(optimalValues)[-length(names(optimalValues))], "AUC")
  optimalValues
}
