#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
compVar= function (demo, eset, variables, ncomp = 10)
{
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  pcaX <- prcomp(eset, scale. = TRUE, center = TRUE)
  pval <- do.call(rbind, lapply(variables, function(i) {
    apply(pcaX$x, 2, function(j) {
      predictor <- demo[, i]
      if (class(predictor) == "factor") {
        if (nlevels(predictor) == 2) {
          coef(summary(lm(as.numeric(j) ~ predictor)))[2,
            "Pr(>|t|)"]
        }
        else {
          anova(lm(as.numeric(j) ~ predictor))["predictor",
            "Pr(>F)"]
        }
      }
      else {
        coef(summary(lm(as.numeric(j) ~ predictor)))[2,
          "Pr(>|t|)"]
      }
    })
  }))
  rownames(pval) <- variables
  colnames(pval) <- paste(colnames(pval), paste0(round(100 *
      (pcaX$sdev^2/sum(pcaX$sdev^2)), 1), "%"), sep = "-")
  pval <- pval[, 1:ncomp]
  pvalheatmap <- pval
  pvalheatmap[pvalheatmap < 0.01] <- 0.01
  pvalheatmap[pvalheatmap > 0.1] <- 1
  pvalheatmap[pvalheatmap > 0.01 & pvalheatmap < 0.05] <- 0.05
  pvalheatmap[pvalheatmap > 0.05 & pvalheatmap < 0.1] <- 0.1
  pvalheatmap[pvalheatmap == "0.01"] <- "p < 0.01"
  pvalheatmap[pvalheatmap == "0.05"] <- "0.01 < p < 0.05"
  pvalheatmap[pvalheatmap == "0.1"] <- "0.05 < p < 0.10"
  pvalheatmap[pvalheatmap == "1"] <- "p > 0.10"
  p <- pvalheatmap %>% as.data.frame %>% mutate(Variable = rownames(.)) %>%
    gather(Threshold, Value, -Variable) %>% mutate(Threshold = factor(Threshold,
      levels = unique(Threshold))) %>% mutate(Variable = factor(as.character(Variable),
        levels = variables)) %>%
    mutate(Value = factor(Value, levels = c("p < 0.01", "0.01 < p < 0.05", "0.05 < p < 0.10", "p > 0.10"))) %>%
    ggplot(aes(Threshold, Variable)) +
    geom_tile(aes(fill = Value), colour = "white") + scale_fill_manual(values = rev(brewer.pal(n = 8,
      name = "Blues")[c(2, 4, 6, 8)])) + customTheme(sizeStripFont = 10,
        xAngle = 40, hjust = 1, vjust = 1, xSize = 10, ySize = 10,
        xAxisSize = 10, yAxisSize = 10) + xlab("") + ylab("")
  return(list(pval = pval, pvalheatmap = pvalheatmap, p = p))
}
