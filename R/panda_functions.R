#' Caculate differential regulatory network statistics
#'
#'
#' @param motif 3 column matrix, 1st column (regVar), 2nd column (gene), 3rd (weight)
#' @param geneExp nxp matrix (gene by sample)
#' @param regulatory expression (regVar by sample)
#' @param phenotype factor specifiying the group of each sample
#' @param nperms # of permutations
#' @export
pandaxBoot = function (motif, geneExp, regExp, phenotype, nperms, threads){
  ## Estimate signal
  diff <- pandax(motif, geneExp, regExp, phenotype)

  ## Estimate noise
  phenotype <- lapply(1:nperms, function(i) sample(phenotype))
  require(parallel)
  cl <- parallel::makeCluster(mc <- getOption("cl.cores",  threads))
  parallel::clusterExport(cl, varlist = c("pandax", "pandaModif", "nperms",
    "normalizeNetwork", "tanimoto", "update.diagonal", "motif", "geneExp",
    "regExp", "phenotype"), envir = environment())
  pandaxSD <- parallel::parLapply(cl, phenotype, function(phenotypei,
    motif, geneExp, regExp, phenotype) {
    library(dplyr); library(tidyr); library(matrixStats);
    pandax(motif=motif, geneExp=geneExp, regExp=regExp, phenotype = phenotypei)
  }, motif, geneExp, regExp, phenotype) %>% amritr::zip_nPure()
  parallel::stopCluster(cl)

  bootStd <- lapply(pandaxSD, function(i){
    bootStd0 <- do.call(rbind, i) %>% as.data.frame %>%
      mutate(regVar = rownames(.)) %>% group_by(regVar) %>%
      summarise_each(funs(sd))
    bootStd <- as.matrix(bootStd0[, -1])
    rownames(bootStd) <- bootStd0$regVar
    bootStd
  })

  ## Estimate singal to noise ratio
  snr <- rbind(cbind(diff$regNetDiff/bootStd$regNetDiff, diff$tfCoopNetworkDiff/bootStd$tfCoopNetworkDiff),
    cbind(diff$geneCoregDiff/bootStd$geneCoregDiff, t(diff$regNetDiff/bootStd$regNetDiff)))
  snr
}

#' Caculate difference between regulatory coefficients of the two groups
#'
#'
#' @param motif 3 column matrix, 1st column (regVar), 2nd column (gene), 3rd (weight)
#' @param geneExp nxp matrix (gene by sample)
#' @param regulatory expression (regVar by sample)
#' @param phenotype factor specifiying the group of each sample
#' @export
pandax = function (motif, geneExp, regExp, phenotype){
  lvls <- levels(phenotype)
  expr1 <- geneExp[, phenotype == lvls[1]]
  expr2 <- geneExp[, phenotype == lvls[2]]
  regexp1 <- cor(t(regExp[, phenotype == lvls[1]])) %>% as.data.frame %>%
    mutate(miRNA = rownames(.)) %>% gather(miRNA1, Exp, -miRNA)
  regexp2 <- cor(t(regExp[, phenotype == lvls[2]])) %>% as.data.frame %>%
    mutate(miRNA = rownames(.)) %>% gather(miRNA1, Exp, -miRNA)
  ## Run PANDA for each phenotype
  pandaGroup1 <- pandaModif(motif = motif, expr = expr1, ppi = regexp1,
    alpha = 0.1, hamming = 0.01, iter = NA)
  pandaGroup2 <- pandaModif(motif = motif, expr = expr2, ppi = regexp2,
    alpha = 0.1, hamming = 0.01, iter = NA)
  ## Gene Co-regulation
  geneCoregDiff <-  pandaGroup2$geneCoreg - pandaGroup1$geneCoreg
  regNetDiff <- pandaGroup2$regulatoryNetwork - pandaGroup1$regulatoryNetwork
  tfCoopNetworkDiff <- pandaGroup2$tfCoopNetwork - pandaGroup1$tfCoopNetwork

  return(list(geneCoregDiff = geneCoregDiff, regNetDiff = regNetDiff, tfCoopNetworkDiff = tfCoopNetworkDiff))
}


#' run PANDA (calculate geneCoreg, regNet and co-operative matrices)
#'
#'
#' @param motif 3 column matrix, 1st column (regVar), 2nd column (gene), 3rd (weight)
#' @param geneExp nxp matrix (gene by sample)
#' @param regulatory expression (regVar by sample)
#' @param phenotype factor specifiying the group of each sample
#' @export
pandaModif = function (motif, expr = NULL, ppi = NULL, alpha = 0.1, hamming = 0.01,
  iter = NA){
  expr <- expr[which(rownames(expr) %in% motif[, 2]), ]
  motif <- motif[which(motif[, 2] %in% rownames(expr)), ]
  expr <- expr[order(rownames(expr)), ]
  num.conditions <- ncol(expr)
  tf.names <- sort(unique(motif[, 1]))
  gene.names <- sort(unique(rownames(expr)))
  num.TFs <- length(tf.names)
  num.genes <- length(gene.names)
  geneCoreg <- cor(t(expr), method = "pearson", use = "pairwise.complete.obs")
  Idx1 = match(motif[, 1], tf.names)
  Idx2 = match(motif[, 2], gene.names)
  Idx = (Idx2 - 1) * num.TFs + Idx1
  regulatoryNetwork = matrix(data = 0, num.TFs, num.genes)
  regulatoryNetwork[Idx] = motif[, 3]
  colnames(regulatoryNetwork) <- gene.names
  rownames(regulatoryNetwork) <- tf.names
  tfCoopNetwork <- diag(num.TFs)
  ppi <- ppi[which(ppi[, 1] %in% tf.names & ppi[, 2] %in% tf.names),
    ]
  Idx1 <- match(ppi[, 1], tf.names)
  Idx2 <- match(ppi[, 2], tf.names)
  Idx <- (Idx2 - 1) * num.TFs + Idx1
  tfCoopNetwork[Idx] <- ppi[, 3]
  Idx <- (Idx1 - 1) * num.TFs + Idx2
  tfCoopNetwork[Idx] <- ppi[, 3]
  colnames(tfCoopNetwork) <- tf.names
  rownames(tfCoopNetwork) <- tf.names
  regulatoryNetwork = normalizeNetwork(regulatoryNetwork)
  tfCoopNetwork = normalizeNetwork(tfCoopNetwork)
  geneCoreg = normalizeNetwork(geneCoreg)
  minusAlpha = 1 - alpha
  step = 0
  hamming_cur = 1
  while (hamming_cur > hamming) {
    if ((!is.na(iter)) && step >= iter) {
      print(paste("Reached maximum iterations, iter =",
        iter), sep = "")
      break
    }
    Responsibility = tanimoto(tfCoopNetwork, regulatoryNetwork)
    Availability = tanimoto(regulatoryNetwork, geneCoreg)
    RA = 0.5 * (Responsibility + Availability)
    hamming_cur = sum(abs(regulatoryNetwork - RA))/(num.TFs *
        num.genes)
    regulatoryNetwork = minusAlpha * regulatoryNetwork +
      alpha * RA
    ppi = tanimoto(regulatoryNetwork, t(regulatoryNetwork))
    ppi = update.diagonal(ppi, num.TFs, alpha, step)
    tfCoopNetwork = minusAlpha * tfCoopNetwork + alpha *
      ppi
    CoReg2 = tanimoto(t(regulatoryNetwork), regulatoryNetwork)
    CoReg2 = update.diagonal(CoReg2, num.genes, alpha, step)
    geneCoreg = minusAlpha * geneCoreg + alpha * CoReg2
    step = step + 1
  }
  return(list(geneCoreg = geneCoreg, regulatoryNetwork = regulatoryNetwork,
    tfCoopNetwork = tfCoopNetwork))
}


#' normalizeNetwork
#'
#'
#' @export
normalizeNetwork = function (X)
{
  X <- as.matrix(X)
  nr = nrow(X)
  nc = ncol(X)
  dm = c(nr, nc)
  mu0 = mean(X)
  std0 = sd(X) * sqrt((nr * nc - 1)/(nr * nc))
  mu1 = rowMeans(X)
  std1 = rowSds(X) * sqrt((nc - 1)/nc)
  mu1 = rep(mu1, nc)
  dim(mu1) = dm
  std1 = rep(std1, nc)
  dim(std1) = dm
  Z1 = (X - mu1)/std1
  mu2 = colMeans(X)
  std2 = colSds(X) * sqrt((nr - 1)/nr)
  mu2 = rep(mu2, each = nr)
  dim(mu2) = dm
  std2 = rep(std2, each = nr)
  dim(std2) = dm
  Z2 = (X - mu2)/std2
  normMat = Z1/sqrt(2) + Z2/sqrt(2)
  Z0 = (X - mu0)/std0
  f1 = is.na(Z1)
  f2 = is.na(Z2)
  normMat[f1] = Z2[f1]/sqrt(2) + Z0[f1]/sqrt(2)
  normMat[f2] = Z1[f2]/sqrt(2) + Z0[f2]/sqrt(2)
  normMat[f1 & f2] = 2 * Z0[f1 & f2]/sqrt(2)
  normMat
}

#' rowSds
#'
#'
#' @export
rowSds = function (x, rows = NULL, cols = NULL, ...)
{
  x <- rowVars(x, rows = rows, cols = cols, ...)
  sqrt(x)
}

#' rowVars
#'
#'
#' @export
rowVars = function (x, rows = NULL, cols = NULL, na.rm = FALSE, center = NULL,
  dim. = dim(x), ...)
{
  dim. <- as.integer(dim.)
  if (is.null(center)) {
    na.rm <- as.logical(na.rm)
    hasNAs <- TRUE
    sigma2 <- .Call("rowVars", x, dim., rows, cols, na.rm,
      hasNAs, TRUE, PACKAGE = "matrixStats")
    return(sigma2)
  }
  if (is.vector(x))
    dim(x) <- dim.
  if (!is.null(rows) && !is.null(cols))
    x <- x[rows, cols, drop = FALSE]
  else if (!is.null(rows))
    x <- x[rows, , drop = FALSE]
  else if (!is.null(cols))
    x <- x[, cols, drop = FALSE]
  dim. <- dim(x)
  if (!is.null(rows))
    center <- center[rows]
  ncol <- ncol(x)
  if (ncol <= 1L) {
    x <- rep(NA_real_, times = nrow(x))
    return(x)
  }
  if (na.rm) {
    nNA <- rowCounts(x, value = NA_real_, na.rm = FALSE)
    n <- ncol - nNA
    hasNA <- any(nNA > 0L)
    if (hasNA) {
      n[n <= 1L] <- NA_integer_
    }
    else {
      na.rm <- FALSE
    }
  }
  else {
    n <- ncol
  }
  x <- x * x
  x <- rowMeans(x, na.rm = na.rm)
  x <- (x - center^2)
  x * (n/(n - 1))
}

#' colSds
#'
#'
#' @export
colSds = function (x, rows = NULL, cols = NULL, ...)
{
  x <- colVars(x, rows = rows, cols = cols, ...)
  sqrt(x)
}

#' colVars
#'
#'
#' @export
colVars = function (x, rows = NULL, cols = NULL, na.rm = FALSE, center = NULL,
  dim. = dim(x), ...)
{
  dim. <- as.integer(dim.)
  if (is.null(center)) {
    dim. <- as.integer(dim.)
    na.rm <- as.logical(na.rm)
    hasNAs <- TRUE
    sigma2 <- .Call("rowVars", x, dim., rows, cols, na.rm,
      hasNAs, FALSE, PACKAGE = "matrixStats")
    return(sigma2)
  }
  if (is.vector(x))
    dim(x) <- dim.
  if (!is.null(rows) && !is.null(cols))
    x <- x[rows, cols, drop = FALSE]
  else if (!is.null(rows))
    x <- x[rows, , drop = FALSE]
  else if (!is.null(cols))
    x <- x[, cols, drop = FALSE]
  dim. <- dim(x)
  if (!is.null(cols))
    center <- center[cols]
  nrow <- nrow(x)
  if (nrow <= 1L) {
    x <- rep(NA_real_, times = ncol(x))
    return(x)
  }
  if (na.rm) {
    nNA <- colCounts(x, value = NA_real_, na.rm = FALSE)
    n <- nrow - nNA
    hasNA <- any(nNA > 0L)
    if (hasNA) {
      n[n <= 1L] <- NA_integer_
    }
    else {
      na.rm <- FALSE
    }
  }
  else {
    n <- nrow
  }
  x <- x * x
  x <- colMeans(x, na.rm = na.rm)
  x <- (x - center^2)
  x * (n/(n - 1))
}

#' tanimoto
#'
#'
#' @export
tanimoto = function (X, Y)
{
  nc = ncol(Y)
  nr = nrow(X)
  dm = c(nr, nc)
  Amat = (X %*% Y)
  Bmat = colSums(Y * Y)
  Bmat = rep(Bmat, each = nr)
  dim(Bmat) = dm
  Cmat = rowSums(X * X)
  Cmat = rep(Cmat, nc)
  dim(Cmat) = dm
  den = (Bmat + Cmat - abs(Amat))
  Amat = Amat/sqrt(den)
  return(Amat)
}

#' update.diagonal
#'
#'
#' @export
update.diagonal = function (diagMat, num, alpha, step)
{
  seqs = seq(1, num * num, num + 1)
  diagMat[seqs] = NaN
  diagstd = rowSds(diagMat, na.rm = TRUE) * sqrt((num - 2)/(num -
      1))
  diagMat[seqs] = diagstd * num * exp(2 * alpha * step)
  return(diagMat)
}
