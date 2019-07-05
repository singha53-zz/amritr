#' integrativePanels()
#'
#' build integrative panels (Concatenation, Ensemble, DIABLO)
#' @param X.train list of nxp datasets of length k
#' @param Y.train - n-vector of class labels
#' @param concat_alpha - 0 < value < 1
#' @param concat_lambda - value that controls the strength of the penalization
#' @param single_alphaList - list of alpha values of length k
#' @param single_lambdaList - list of lambda values of length k
#' @export
integrativePanels = function(X.train, Y.train, X.test, Y.test, concat_alpha, concat_lambda, single_alphaList, single_lambdaList, M, iter, threads){

  ## Concatenation-Enet
  if(is.null(X.test)){
    result <- amritr::enet(X = do.call(cbind, X.train), Y = Y.train, alpha = concat_alpha,
      family="multinomial", lambda = concat_lambda, X.test = NULL,
      Y.test = NULL, filter = "none", topranked = 50)
    concat_enetPanel <- lapply(X.train, function(i){ intersect(colnames(i), result$enet.panel)})

    ## Estimate panel performance using cross-validation
    cv <- amritr::perf.enet(result, validation = "Mfold", M = M, iter = iter, threads = cpus, progressBar = FALSE)
    concat_enetErrTrain_tuneConcat <- filter(cv$perf, ErrName == "BER")[-1]
    concat_enetErrTest_tuneConcat <- c(NA, NA)
    names(concat_enetErrTest_tuneConcat) <- names(concat_enetErrTrain_tuneConcat)

    ## Concateantion panel error rate
    concatErr <- rbind(concat_enetErrTrain_tuneConcat, concat_enetErrTest_tuneConcat) %>%
      mutate(Set = c("Train", "Test"))
    concatErr$Method <- "Concatenation"
  } else {
    if(length(X.train) == length(X.test)){
      result <- amritr::enet(X = do.call(cbind, X.train), Y = Y.train, alpha = concat_alpha,
        family="multinomial", lambda = concat_lambda, X.test = do.call(cbind, X.test),
        Y.test = Y.test, filter = "none", topranked = 50)
      concat_enetPanel <- lapply(X.train, function(i){ intersect(colnames(i), result$enet.panel)})

      ## Estimate panel performance using cross-validation
      cv <- amritr::perf.enet(result, validation = "Mfold", M = M, iter = iter, threads = cpus, progressBar = FALSE)
      concat_enetErrTrain_tuneConcat <- filter(cv$perf, ErrName == "BER")[-1]
      concat_enetErrTest_tuneConcat <- c(result$perfTest["BER"], NA)
      names(concat_enetErrTest_tuneConcat) <- names(concat_enetErrTrain_tuneConcat)

      ## Concateantion panel error rate
      concatErr <- rbind(concat_enetErrTrain_tuneConcat, concat_enetErrTest_tuneConcat) %>%
        mutate(Set = c("Train", "Test"))
      concatErr$Method <- "Concatenation"
    } else {
      result <- amritr::enet(X = do.call(cbind, X.train), Y = Y.train, alpha = concat_alpha,
        family="multinomial", lambda = concat_lambda, X.test = NULL,
        Y.test = NULL, filter = "none", topranked = 50)
      concat_enetPanel <- lapply(X.train, function(i){ intersect(colnames(i), result$enet.panel)})

      ## Estimate panel performance using cross-validation
      cv <- amritr::perf.enet(result, validation = "Mfold", M = M, iter = iter, threads = cpus, progressBar = FALSE)
      concat_enetErrTrain_tuneConcat <- filter(cv$perf, ErrName == "BER")[-1]
      concat_enetErrTest_tuneConcat <- c(NA, NA)
      names(concat_enetErrTest_tuneConcat) <- names(concat_enetErrTrain_tuneConcat)

      ## Concateantion panel error rate
      concatErr <- rbind(concat_enetErrTrain_tuneConcat, concat_enetErrTest_tuneConcat) %>%
        mutate(Set = c("Train", "Test"))
      concatErr$Method <- "Concatenation"
    }
  }

  ## Ensemble-Enet
  ensembleMod <- amritr::ensembleEnet(X.train = X.train, Y.train = Y.train, alphaList = single_alphaList,
    lambdaList = single_lambdaList, X.test = X.test, Y.test = Y.test, filter = "none", topranked = 50)
  ensembleResult <- ensembleMod$result %>% zip_nPure()
  ensemble_enetPanel <- ensembleResult$enet.panel

  ## Estimate panel performance using cross-validation
  ensembleTrain <- perfEnsemble(object=ensembleMod, validation = "Mfold", M = M, iter = iter,
    threads = cpus, progressBar = TRUE)
  ensemble_enetLength <- lapply(ensemble_enetPanel, length)
  ensemble_enetErrTrain_tuneEnsemble <- filter(ensembleTrain$perf, ErrName == "BER")[-1]
  ensemble_enetErrTest_tuneEnsemble <- c(ensembleMod$perfTest["BER"], NA)
  names(ensemble_enetErrTest_tuneEnsemble) <- names(ensemble_enetErrTrain_tuneEnsemble)

  ## Ensemble panel error rate
  ensembleErr <- rbind(ensemble_enetErrTrain_tuneEnsemble, ensemble_enetErrTest_tuneEnsemble) %>%
    mutate(Set = c("Train", "Test"))
  ensembleErr$Method <- "Ensemble"

  # DIABLO
  design <- setDesign(X.train, corCutOff = 0.6, plotMat = FALSE)
  ncomp <- nlevels(Y.train) - 1
  list.keepX <- lapply(ensemble_enetLength, function(i){
    rep(round(i/ncomp, 0), ncomp)
  })
  TCGA.block.splsda = block.splsda(X = X.train, Y = Y.train,
    ncomp = ncomp, keepX = list.keepX, design = design,
    scheme = "centroid")
  diabloPanel <- lapply(TCGA.block.splsda$loadings[-(length(X.train)+1)], function(x)
    unique(as.character(as.matrix(apply(x, 2, function(i) names(i)[which(i != 0)])))))

  ## training error
  cv <- perf(TCGA.block.splsda, validation = "Mfold", folds = M, cpus = cpus, nrepeat = iter)
  err <- extractErr(cv)
  err$Comp <- factor(err$Comp, levels = paste("comp", 1:ncomp, sep = "."))
  diablo_enetErrTrain <- err %>% filter(Type == "centroids.dist", Class == "Overall.BER",
    EnsembleMode == "wMajVote", Dataset == "DIABLO", Comp == paste("comp", ncomp, sep = "."))
  diablo_enetErrTrain <- diablo_enetErrTrain[c("meanErr", "sdErr")]

  ## test error
  if(is.null(X.test)){
    diablo_enetErrTest <- c(NA, NA)
    names(diablo_enetErrTest) <- names(diablo_enetErrTrain)
    ## DIABLO panel error rate
    diabloErr <- rbind(diablo_enetErrTrain, diablo_enetErrTest)
    colnames(diabloErr) <- c("Mean", "SD")
    diabloErr$Set <- c("Train", "Test")
    diabloErr$Method <- "DIABLO"
  } else {
    diabloTest <- predict(TCGA.block.splsda, X.test, method = "all")
    diabloTestConsensus <- lapply(diabloTest$WeightedVote, function(i){
      predY <- apply(i, 2, function(z){
        temp <- table(factor(z, levels = levels(Y.test)), Y.test)
        diag(temp) <- 0
        error = c(colSums(temp)/summary(Y.test), sum(temp)/length(Y.test), mean(colSums(temp)/summary(Y.test)))
        names(error) <- c(names(error)[1:nlevels(Y.test)], "ER", "BER")
        error
      })
    })
    diablo_enetErrTest <- c(diabloTestConsensus$centroids.dist["BER", paste("comp", ncomp, sep = "")], NA)
    names(diablo_enetErrTest) <- names(diablo_enetErrTrain)
    ## DIABLO panel error rate
    diabloErr <- rbind(diablo_enetErrTrain, diablo_enetErrTest)
    colnames(diabloErr) <- c("Mean", "SD")
    diabloErr$Set <- c("Train", "Test")
    diabloErr$Method <- "DIABLO"
  }

  ## Error rates
  error <- rbind(concatErr, ensembleErr, diabloErr)

  ## Panels
  panels = list(Concatenation = concat_enetPanel, Ensemble = ensemble_enetPanel,
    DIABLO = diabloPanel)

  return(list(error = error, panels = panels))
}

#' networkStats()
#'
#' compuate network statistics given a adjacency matrix
#' @param adjMat p x p matrix
#' @param mode = "graph"
#' @export
networkStats = function(adjMat, mode = "graph"){
  return(1)
}

#' graphIndices()
#'
#' build integrative panels (Concatenation, Ensemble, DIABLO)
#' @param panels - list of panels (character vector corresponding to column names of X.train)
#' @param X.train - list of datasets (n x p)
#' @param cut-off - of person correlatio
#' @param concat_lambda - value that controls the strength of the penalization
#' @param single_alphaList - list of alpha values of length k
#' @param single_lambdaList - list of lambda values of length k
#' @export
graphIndices = function(panels = panels, X.train = X.train, cutoff = cutoff){
  library(dplyr)
  ## Determine adjacency matrices
  adjMat <- lapply(panels, function(i){
    adjMat = cor(do.call(cbind, mapply(function(x, y){
      y[, x]
    }, x = i, y = X.train)))
    adjMat[abs(adjMat) < cutoff] <- 0
    adjMat[abs(adjMat) > cutoff] <- 1
    adjMat
  })

  # Estimate graph statistics
  graphs <- lapply(adjMat, function(x) amritr::networkStats(adjMat = x))

  ## graphDensity, recipScore, transScore, cliques
  gIndices <- lapply(graphs, function(i){
    data.frame(Value = c(i$graphDensity, i$recipScore["edgewise.lrr"], i$transScore["weakcensus"], length(i$cliques)),
      Statistic = c("Graph Density", "Edgewise.lrr", "Transitivity", "NumOfCliques"))

  })
  gIndicesDat <- do.call(rbind, gIndices) %>% mutate(Method = rep(names(gIndices), each = 4))

  # Dyads and Triads
  dyads <- do.call(rbind, lapply(graphs, function(x){x$dyadCensus}))
  dyads <- as.data.frame(dyads)[,-2] %>% mutate(Method = names(graphs)) %>%
    tidyr::gather(Type, Number, -Method)
  triads <- do.call(rbind, lapply(graphs, function(x){x$triadCensus}))
  triads <- as.data.frame(triads) %>% mutate(Method = names(graphs)) %>%
    tidyr::gather(Type, Number, -Method)

  # Number of isolates
  isolatedFeat <- lapply(graphs, function(x){ x$isolatedVertices})
  isolates <- do.call(rbind, lapply(isolatedFeat, function(x){
    unlist(lapply(X.train, function(y){
      length(intersect(x, colnames(y)))
    }))
  })) %>% as.data.frame %>% mutate(Method = names(graphs)) %>%
    tidyr::gather(Dataset, NumOfIsolates, -Method)

  return(list(adjMat = adjMat, gIndicesDat = gIndicesDat, dyads = dyads, triads = triads, isolates = isolates))
}
