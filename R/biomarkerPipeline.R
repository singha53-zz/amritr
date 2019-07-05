#' biomarkerPipeline
#'
#' build various biomarker panels
#' @param X nxp matrix - training dataset
#' @param Y binary variable
#' @param topranked - top features ranked using p-value to build a classification panel
#' @param validation (Mfold/loocv)
#' @param M - # of folds
#' @param iter - Number of times to repeat cross-validation
#' @param threads - number of nodes (each CV runs on a separate node)
#' @param progressBar = (TRUE/FALSE)
#' @param pathways - list of data.frame containing pathway to genes mapping
#' @export
biomarkerPipeline = function(X = X, Y = Y, topranked = 50, validation = "Mfold", M = 5, iter = 1, threads = 1, progressBar = TRUE, pathways = pathways){
  ## Enet panel
  alphaSeq <- seq(0, 0.9, by = 0.1)
  enetAUC <- lapply(c("none", "p.value"), function(i){
    aucEnet <- lapply(alphaSeq, function(x){
      result = enet(X, Y, alpha=x, lambda = NULL, family = "binomial", X.test = NULL,
        Y.test = NULL, filter = i, topranked = topranked)
      cv <- perf.enet(object = result, validation = validation, M = M, iter = iter,
        threads = threads, progressBar = progressBar)
      filter(cv$perf, ErrName == "AUC") %>% dplyr::select(Mean, SD) %>%
        mutate(Panel = paste("Enet", x, sep = "_"), Genes = paste(result$enet.panel, collapse = "_"))
    })
    as.data.frame(do.call(rbind, aucEnet))
  })
  enetAUC <- do.call(rbind, enetAUC)
  enetAUC$Panel <- paste(enetAUC$Panel, rep(c("none", "p.value"), each = length(alphaSeq)), sep = "_")

  ## Random forest panel
  rfAUC <- as.data.frame(do.call(rbind, lapply(c("none", "p.value"), function(i){
    res.rf <- rforest(X, Y, X.test = NULL, Y.test = NULL, family = "binomial", filter = i, topranked = 10)
    cvRf <- perf.rf(object = res.rf, validation = validation, M = M, iter = iter,
      threads = threads, progressBar = progressBar)
    filter(cvRf$perf, ErrName == "AUC") %>% dplyr::select(Mean, SD) %>%
      mutate(Panel = paste("RF", i, sep = "_"), Genes = paste(res.rf$rf.panel, collapse = "_"))
  })))

  ## svm
  svmAUC <- as.data.frame(do.call(rbind, lapply(c("none", "p.value"), function(i){
    res.svm <- supportVectorMachine(X, Y, X.test = NULL, Y.test = NULL, family = "binomial",
      filter = i, topranked = topranked)
    cvSvm <- perf.svm(object = res.svm, validation = validation, M = M, iter = iter,
      threads = threads, progressBar = TRUE)
    filter(cvSvm$perf, ErrName == "AUC") %>% dplyr::select(Mean, SD) %>%
      mutate(Panel = paste("SVM", i, sep = "_"), Genes = paste(res.svm$svm.panel, collapse = "_"))
  })))

  ## single biomarkers
  singlePerf <- do.call(rbind, lapply(1 : ncol(X), function(i){
    res.single <- glmPanel(X[, i, drop = FALSE], Y, X.test = NULL, Y.test = NULL)
    cvSingle <- perf.glm(object = res.single, validation = validation, M = M, iter = iter,
      threads = threads, progressBar = TRUE)
    filter(cvSingle$perf, ErrName == "AUC")
  })) %>% mutate(Panel = colnames(X), Genes = colnames(X)) %>% dplyr::select(Mean, SD, Panel, Genes)

  if(!is.null(pathways)){
    ## change pathway datafrome to list
    pathways <- lapply(split(pathways, pathways$Database), function(i){
      lapply(split(i, i$Term), function(j){
        j$Genes
      })
    })
    ## pathway biomarkers
    convertPathwayToDataFrame = function(X, pathway){
      library(dplyr)
      pathwayList <- lapply(pathway, function(i){
        intersect(unlist(i), colnames(X))
      })
      pathwayList <- pathwayList[unlist(lapply(pathwayList, length)) != 0]
      pathwayList <- lapply(pathwayList, function(j){
        paste(as.character(j)[order(as.character(j))], collapse = "_")
      })

      pathways <- do.call(rbind, pathwayList) %>% as.data.frame %>% mutate(Pathway = rownames(.)) %>%
        group_by(V1) %>% dplyr::summarise_each(funs(paste(., collapse = ";")))
      pathways
    }
    Dat <- lapply(pathways, function(i){
      convertPathwayToDataFrame(X, i)
    })

    pathwayAUC <- lapply(Dat, function(i){
      pathwayDatabaseBiomarkers <- do.call(rbind, lapply(as.character(i$V1), function(x){
        res.single <- rforest(X[, unlist(strsplit(x, "_")), drop = FALSE], Y, X.test = NULL, Y.test = NULL, family = "binomial",
          filter = "none", topranked = 50)
        cvSingle <- perf.rf(object = res.single, validation = validation, M = M, iter = iter,
          threads = threads, progressBar = progressBar)
        filter(cvSingle$perf, ErrName == "AUC")
      }))
      cbind(i, pathwayDatabaseBiomarkers) %>% mutate(Genes = V1, Panel = Pathway) %>%
        dplyr::select(Mean, SD, Panel, Genes)
    })
    pathwayAUCList <- do.call(rbind, pathwayAUC)
    pathwayAUCList$Panel <- paste(rep(names(Dat), unlist(lapply(Dat, nrow))), pathwayAUCList$Panel, sep = "_")

    allPanels <- rbind(enetAUC, rfAUC, svmAUC, singlePerf, pathwayAUCList)
    allPanels$Type <- rep(c("Enet", "RF", "SVM", "GLM_single", paste(names(Dat), "RF", sep = "_")), c(nrow(enetAUC), nrow(rfAUC), nrow(svmAUC), nrow(singlePerf), unlist(lapply(Dat, nrow))))
  } else {
    allPanels <- rbind(enetAUC, rfAUC, svmAUC, singlePerf)
    allPanels$Type <- rep(c("Enet", "RF", "SVM", "GLM_single"), c(nrow(enetAUC), nrow(rfAUC), nrow(svmAUC), nrow(singlePerf)))
  }
  allPanels
}
