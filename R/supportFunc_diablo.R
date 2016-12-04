#' clean up diablo perf() output
#'
#' @param cv - perf(obj) - obj-block.splsda object
#' @export
extractErr = function(cv){
## Single Omics
er0 <- cv$error.rate.all %>% zip_nPure() %>% lapply(., function(x){
  do.call(rbind, x) %>% data.frame %>% mutate(Comp = rep(rownames(x[[1]]), each = length(x)), Iteration = rep(1:length(x), each = nrow(x[[1]]))) %>%
    gather(Type, ErrorRate, -c(Comp:Iteration))
})
er <- do.call(rbind, er0) %>% data.frame %>% mutate(Dataset = rep(names(cv$error.rate.all[[1]]), each = nrow(er0[[1]])))
er$EnsembleMode <- "SingleOmics"
er$Class <- "Overall.ER"
er$Comp <- gsub(" ", ".", er$Comp)

a <- cv$error.rate.per.class.all %>% zip_nPure()
classErr0 <- lapply(a, function(i){
  classErr0 <- i %>% zip_nPure() %>% lapply(., function(x){
    do.call(rbind, x) %>% data.frame %>% mutate(Class = rep(rownames(x[[1]]), each = length(x)), Iteration = rep(1:length(x), each = nrow(x[[1]]))) %>%
      gather(Comp, ErrorRate, -c(Class:Iteration))
  })
  classErr <- do.call(rbind, classErr0) %>% data.frame %>% mutate(Type = rep(names(i[[1]]), each = nrow(classErr0[[1]])))
  classErr$EnsembleMode <- "SingleOmics"
  classErr
})
classErr <- do.call(rbind, classErr0) %>% data.frame %>% mutate(Dataset = rep(names(classErr0), each = nrow(classErr0[[1]])))
ber <- classErr %>% group_by(Iteration, Comp, Type, EnsembleMode, Dataset) %>% summarise(ErrorRate = mean(ErrorRate))
ber$Class <- "Overall.BER"
singleOmics <- rbind(er, classErr[, colnames(er)], as.data.frame(ber[, colnames(er)]))

## Average prediction
avgPred <- do.call(rbind, cv$AveragedPredict.error.rate.all) %>% data.frame %>% mutate(Class = rep(rownames(cv$AveragedPredict.error.rate.all[[1]]), each = length(cv$AveragedPredict.error.rate.all)), Iteration = rep(1:length(cv$AveragedPredict.error.rate.all), each = nrow(cv$AveragedPredict.error.rate.all[[1]]))) %>%
  gather(Comp, ErrorRate, -c(Class:Iteration))
avgPred$Type  = rep("max.dist", nrow(avgPred))
avgPred$EnsembleMode <- "AvgPred"
## weighted average prediction
WeightedAvgPred <- do.call(rbind, cv$WeightedPredict.error.rate.all) %>% data.frame %>% mutate(Class = rep(rownames(cv$WeightedPredict.error.rate.all[[1]]), each = length(cv$WeightedPredict.error.rate.all)), Iteration = rep(1:length(cv$WeightedPredict.error.rate.all), each = nrow(cv$WeightedPredict.error.rate.all[[1]]))) %>%
  gather(Comp, ErrorRate, -c(Class:Iteration))
WeightedAvgPred$Type  = rep("max.dist", nrow(WeightedAvgPred))
WeightedAvgPred$EnsembleMode <- "wAvgPred"
## majority vote
majVote0 <- cv$MajorityVote.error.rate.all %>% zip_nPure() %>% lapply(., function(x){
  do.call(rbind, x) %>% data.frame %>% mutate(Class = rep(rownames(x[[1]]), each = length(x)), Iteration = rep(1:length(x), each = nrow(x[[1]]))) %>%
    gather(Comp, ErrorRate, -c(Class:Iteration))
})
majVote <- do.call(rbind, majVote0) %>% data.frame %>% mutate(Type = rep(names(cv$MajorityVote.error.rate.all[[1]]), each = nrow(majVote0[[1]])))
majVote$EnsembleMode <- "MajVote"
## weighted majority vote
wMajVote0 <- cv$WeightedVote.error.rate.all %>% zip_nPure() %>% lapply(., function(x){
  do.call(rbind, x) %>% data.frame %>% mutate(Class = rep(rownames(x[[1]]), each = length(x)), Iteration = rep(1:length(x), each = nrow(x[[1]]))) %>%
    gather(Comp, ErrorRate, -c(Class:Iteration))
})
wMajVote <- do.call(rbind, wMajVote0) %>% data.frame %>% mutate(Type = rep(names(cv$WeightedVote.error.rate.all[[1]]), each = nrow(wMajVote0[[1]])))
wMajVote$EnsembleMode <- "wMajVote"

multiOmics <- rbind(avgPred, WeightedAvgPred, majVote, wMajVote)
multiOmics$Dataset <- "DIABLO"
multiOmics <- multiOmics[, colnames(singleOmics)]

rbind(singleOmics, multiOmics) %>%
  dplyr::group_by(Type, Class, Comp, EnsembleMode, Dataset) %>% dplyr::summarise(meanErr = mean(ErrorRate), sdErr = sd(ErrorRate, na.rm = TRUE))
}

#' clean up diablo perf() output
#'
#' @param cv - perf(obj) - obj-block.splsda object
#' @export
setDesign = function(X, corCutOff, plotMat = FALSE){
  x.xList <- list()
  for(i in 1:length(X)){
    corDat <- rep(0, length(X))
    names(corDat) <- paste("cor", names(X)[i], names(X), sep = "_")
    for(j in 1:length(X)){
      result <- pls(X = X[[i]], Y = X[[j]], ncomp = 1)
      corDat[j] <- as.numeric(cor(result$variates$X, result$variates$Y))
    }
    x.xList[[i]] <- corDat
  }
  corMat <- do.call(rbind, x.xList)
  rownames(corMat) <- colnames(corMat) <- names(X)

  design <- matrix(0, nrow = nrow(corMat), ncol = ncol(corMat))
  design[corMat > corCutOff] <- 1
  diag(design) <- 0
  rownames(design) <- colnames(design) <- names(X)

  if(plotMat == TRUE){
    corrplot(corMat)
    corrplot(corMat,add=TRUE, type="upper", method="ell",order="original",
      diag=FALSE,tl.pos="d", cl.pos="n")
    corrplot(corMat,add=TRUE, type="lower", method="number",order="original",
      diag=TRUE,tl.pos="d", cl.pos="n")
  }
  design
}
