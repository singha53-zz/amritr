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
#' @param data input dataset (all variables of a given dataset must be of the same type)
#' @param group string stating the column name of the grouping variable
#' @details
#' This function performs either a parametric (linear model) or non-parametric test (Wilcoxon or Kruskal-Wallis) for continuous data, or the Chi-Square test for categorical variables. The diagnostic assumptions of linear models are assessed using global statistic that assessing the four assumptions of linear models (skewness, kurtosis, link function, homoscadasticity (based on a test for heteroscadasticity)) based on the gvlma R-library, \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2820257/}).
#' @export
hypothesisTests = function (data, group, details = FALSE){
  library("lmtest")
  if (!is.factor(data[, group]))
    stop("group variable must be a factor!")
  classes <- droplevels(data[, group])
  if (nlevels(classes) < 2)
    stop("at least 2 levels in the group variable are required")
  data <- data[, setdiff(colnames(data), group), drop = FALSE]
  if (all(sapply(data, is.numeric))) {
    isCont = TRUE
  } else {
    isCont = FALSE
  }
  if (isTRUE(isCont)) {
    ## Check linear model assumptions
    diagnostics = data %>%
      dplyr::mutate(group = classes) %>%
      tidyr::gather(Var, Value, -group) %>%
      dplyr::group_by(Var) %>%
      nest() %>%
      dplyr::mutate(assumptions = purrr::map(data, ~{
        fit <- lm(Value ~ group, data = .)
        assumptions <- gvlma::gvlma(fit)
        assumptions <- assumptions$GlobalTest[2:6] %>%
          do.call(rbind, .) %>%
          as.data.frame %>%
          mutate(Test = c("Global Stat", "Skewness", "Kurtosis", "Link Function", "Heteroscedasticity"),
            lmDiagnostics = ifelse(pvalue < 0.05, "Assumptions NOT satisfied!", "Assumptions acceptable."))
        assumptions$Parametric <- paste("lm", signif(pf(summary(fit)$fstatistic[1],
          summary(fit)$fstatistic[2],
          summary(fit)$fstatistic[3],lower.tail=FALSE), 2), sep="_")
        if(nlevels(classes) == 2){
          assumptions$NonParametric <- paste("wilcoxon", signif(wilcox.test(Value ~ group, data = .)$p.value, 2), sep="_")
        } else {
          assumptions$NonParametric <- paste("Kruskal.Wallis", signif(kruskal.test(Value ~ group, data = .)$p.value, 2), sep="_")
        }
        assumptions
      })) %>%
      dplyr::mutate(summary = purrr::map(data, ~{
        group_by(., group) %>% dplyr::summarize(Mean.SD = paste(c(round(mean(Value,
          na.rm = TRUE), 1), round(sd(Value, na.rm = TRUE),
            1)), collapse = "+/-")) %>%
          tidyr::spread(group, Mean.SD)
      })) %>%
      unnest(summary) %>% unnest(assumptions) %>% unnest() %>%
      dplyr::select(-c(Value, pvalue, Decision)) %>%
      gather(TestType, Pval, Parametric:NonParametric) %>%
      separate(Pval, c("Method", "p.value"), "_") %>%
      mutate(Decision = ifelse(as.numeric(p.value) < 0.05, "Significant", "Non-Significant"))

    ## based on diagnostic run parameteric or non-parametric test
    if(details){
      split(diagnostics, diagnostics$Decision)
    } else {
      diagnostics <- diagnostics %>% filter(Test == "Global Stat") %>%
        group_by(Var) %>%
        filter(lmDiagnostics == "Assumptions acceptable." & TestType == "Parametric" | lmDiagnostics == "Assumptions NOT satisfied!" & TestType == "NonParametric")
      split(diagnostics, diagnostics$Decision)
    }

  } else {
    result <- apply(data, 2, function(i) {
      sigTest <- data.frame(chisq = chisq.test(table(i,
        classes))$p.value)
      if (sigTest$chisq < 0.05) {
        sigTest$Decision <- "Significant"
      } else {
        sigTest$Decision <- "Not_Significant"
      }
      sigTest
    }) %>% do.call(rbind, .)
    summary <- data %>% dplyr::mutate(group = classes) %>%
      tidyr::gather(Var, lvl, -group) %>% dplyr::group_by(Var,
        group) %>% dplyr::summarise(lvl = paste(paste(names(table(lvl)),
          round(100 * table(lvl)/sum(table(lvl)), 1), sep = "_"),
          collapse = "/")) %>% tidyr::spread(group, lvl) %>%
      as.data.frame
    diagnostics <- cbind(summary, result[summary$Var, ])
    split(diagnostics, diagnostics$Decision)
  }
}
