#' split demographics into two datasets (continuous and categorical variables)
#'
#' @param demo demographics dataset
#' @param group string stateting the column name of the grouping variable
#' @param trim remove variables with 20% missing data
#' @export
splitData = function(demo, group, trim = 0.8){
  if(nrow(demo) < 10)
    stop("too few samples")

  ## Class balance
  print("Number of samples in each group")
  print(table(demo[, group]))
  imbal <- scale(table(demo[, group]), center = FALSE, scale = nrow(demo)) %>% drop
  if(sum(imbal < 0.4) > 0)
    print("Class imbalance exists")

  ## Remove missing demo
  missing <- colSums(!is.na(demo))/nrow(demo)
  varKeep <- missing[missing >= trim]

  ## List missing variables
  if(length(setdiff(names(missing), names(varKeep))) > 0){
    print("Missing variables and percentage of avaiable data")
    print(missing[missing < trim][order(missing[missing < trim])])
  }

  ## remove missing variables
  data <- demo[, names(varKeep)]

  ## Split data into 2 datasets (categorical and continuous variables)
  dataList <- lapply(1 : ncol(data), function(i){
    if(nlevels(factor(as.character(data[, i]))) < 10){
      x <- factor(as.character(data[, i]))
    } else {
      x <- as.numeric(data[, i])
    }
    x
  })
  names(dataList) <- colnames(data)
  data.new <- dataList %>% do.call(cbind.data.frame, .)

  ## data with categorical variables
  data.cat <- data.new[, !sapply(data.new, is.numeric)]

  ## how many variables with one level?
  oneLevel <- which(dataList[colnames(data.cat)] %>% sapply(., function(i) nlevels(i)) == 1)
  if(length(oneLevel) > 1){
    print(paste0("The categorical variable ",
      paste(names(which(dataList[colnames(data.cat)] %>% sapply(., function(i) nlevels(i)) == 1)), collapse = "/"),
      " only has 1 level; User should remove this variable"))
  } else {
    print("All categorial variables have more than 1 level: Good!")
  }

  ## data with continuous variables
  data.cont <- data.new[, sapply(data.new, is.numeric), drop = FALSE]

  ## Summariize
  summary <- data.frame(Total_nVar = ncol(demo),
    missing_nVar = (ncol(demo) - length(varKeep)),
    nonMissing_nVar = length(varKeep),
    nonMissing_cat_nVar = ncol(data.cat),
    nonMissing_cont_nVar = ncol(data.cont))

  return(list(summary=summary, data.cat=data.cat, data.cont=data.cont))
}

#' perform hypothesis test based on variable type
#'
#' @param data input dataset
#' @param group string stateting the column name of the grouping variable
#' @export
hypothesisTests = function(data, group){
  library("lmtest")
  if (!is.factor(data[, group]))
    stop("group variable must be a factor!")
  classes <- droplevels(data[, group])
  if (nlevels(classes) < 2)
    stop("at least 2 levels in the group variable are required")
  data <- data[, setdiff(colnames(data), group), drop = FALSE]
  if (all(sapply(data, is.numeric))) {
    isCont = TRUE
  }
  else {
    isCont = FALSE
  }
  if (isTRUE(isCont)) {
    summary <- data %>% dplyr::mutate(group = classes) %>%
      tidyr::gather(Var, Value, -group) %>% dplyr::group_by(group,
        Var) %>% dplyr::summarize(Mean.SD = paste(c(round(mean(Value,
          na.rm = TRUE), 1), round(sd(Value, na.rm = TRUE),
            1)), collapse = "+/-")) %>% tidyr::spread(group,
              Mean.SD) %>% as.data.frame
    if (nlevels(classes) == 2) {
      result <- apply(data, 2, function(i) {
        fit <- lm(i ~ classes)
        sigTest <- data.frame(effectSize = fit$coefficients[2],
          lm.Pval = coef(summary(fit))[2, "Pr(>|t|)"],
          wilcoxon.Pval = wilcox.test(i[classes == levels(classes)[1]],
            i[classes == levels(classes)[2]])$p.value,
          Bartlett.Test = bartlett.test(i ~ classes)$p.value,
          Breusch.Pagan.Test = bptest(fit)$p.value, Shapiro.Test = shapiro.test(fit$residuals)$p.value)
        if (sigTest$Bartlett.Test < 0.05) {
          sigTest$Bartlett.Test_HO_ConstantVar = "Reject_Null"
        }
        else {
          sigTest$Bartlett.Test_HO_ConstantVar = "Dont_Reject_Null"
        }
        if (sigTest$Breusch.Pagan.Test < 0.05) {
          sigTest$Breusch.Pagan.Test_HO_ConstantVar = "Reject_Null"
        }
        else {
          sigTest$Breusch.Pagan.Test_HO_ConstantVar = "Dont_Reject_Null"
        }
        if (sigTest$Shapiro.Test < 0.05) {
          sigTest$Shapiro.Test_HO_normal = "Reject_Null"
        }
        else {
          sigTest$Shapiro.Test_HO_normal = "Dont_Reject_Null"
        }
        assumptions <- table(as.character(sigTest[, c("Bartlett.Test_HO_ConstantVar",
          "Breusch.Pagan.Test_HO_ConstantVar", "Shapiro.Test_HO_normal")]))
        if (sum(names(assumptions) %in% "Reject_Null")) {
          sigTest$WhichTest <- "wilcoxon"
        }
        else {
          sigTest$WhichTest <- "lm"
        }
        if (sigTest$WhichTest == "lm") {
          if (sigTest$lm.Pval < 0.05) {
            sigTest$Decision <- "Significant"
          }
          else {
            sigTest$Decision <- "Not_Significant"
          }
        }
        else {
          if (sigTest$wilcoxon.Pval < 0.05) {
            sigTest$Decision <- "Significant"
          }
          else {
            sigTest$Decision <- "Not_Significant"
          }
        }
        sigTest
      }) %>% do.call(rbind, .)
    }
    else {
      result <- apply(data, 2, function(i) {
        fit <- aov(i ~ classes)
        sigTest <- data.frame(effectSize = fit$coefficients[2],
          anova.Pval = summary(fit)[[1]][1, "Pr(>F)"], Kruskal.Pval = kruskal.test(i ~
              classes)$p.value, Bartlett.Test = bartlett.test(i ~
                  classes)$p.value, Breusch.Pagan.Test = bptest(fit)$p.value,
          Shapiro.Test = shapiro.test(fit$residuals)$p.value)
        if (sigTest$Bartlett.Test < 0.05) {
          sigTest$Bartlett.Test_HO_ConstantVar = "Reject_Null"
        }
        else {
          sigTest$Bartlett.Test_HO_ConstantVar = "Dont_Reject_Null"
        }
        if (sigTest$Breusch.Pagan.Test < 0.05) {
          sigTest$Breusch.Pagan.Test_HO_ConstantVar = "Reject_Null"
        }
        else {
          sigTest$Breusch.Pagan.Test_HO_ConstantVar = "Dont_Reject_Null"
        }
        if (sigTest$Shapiro.Test < 0.05) {
          sigTest$Shapiro.Test_HO_normal = "Reject_Null"
        }
        else {
          sigTest$Shapiro.Test_HO_normal = "Dont_Reject_Null"
        }
        assumptions <- table(as.character(sigTest[, c("Bartlett.Test_HO_ConstantVar",
          "Breusch.Pagan.Test_HO_ConstantVar", "Shapiro.Test_HO_normal")]))
        if (sum(names(assumptions) %in% "Reject_Null")) {
          sigTest$WhichTest <- "Kruskal-Wallis"
        }
        else {
          sigTest$WhichTest <- "ANOVA"
        }
        if (sigTest$WhichTest == "ANOVA") {
          if (sigTest$anova.Pval < 0.05) {
            sigTest$Decision <- "Significant"
          }
          else {
            sigTest$Decision <- "Not_Significant"
          }
        }
        else {
          if (sigTest$Kruskal.Pval < 0.05) {
            sigTest$Decision <- "Significant"
          }
          else {
            sigTest$Decision <- "Not_Significant"
          }
        }
        sigTest
      }) %>% do.call(rbind, .)
    }
    cbind(summary, result[summary$Var, ]) %>% arrange(effectSize)
  }
  else {
    result <- apply(data, 2, function(i) {
      sigTest <- data.frame(chisq.Pval = chisq.test(table(i,
        classes))$p.value)
      if (sigTest$chisq.Pval < 0.05) {
        sigTest$Decision <- "Significant"
      }
      else {
        sigTest$Decision <- "Not_Significant"
      }
      sigTest
    }) %>% do.call(rbind, .)
    summary <- data %>% dplyr::mutate(group = classes) %>%
      tidyr::gather(Var, lvl, -group) %>% dplyr::group_by(Var,
        group) %>% dplyr::summarise(lvl = paste(paste(names(table(lvl)),
          round(100*table(lvl)/sum(table(lvl)), 1), sep = "_"), collapse = "/")) %>%
      tidyr::spread(group, lvl) %>% as.data.frame
    cbind(summary, result[summary$Var, ])
  }
}

