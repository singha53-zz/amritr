#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
extendedBIC = function (gamma, omegahat, S, n)
{
  p = nrow(omegahat)
  es = sum(omegahat[upper.tri(omegahat)] != 0)
  return(-log(det(omegahat)) + sum(diag(omegahat %*% S)) +
      es * (log(n)/n) + es * gamma * (4 * log(p)/n))
}

#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
fast.dingo = function (dat, x, rhoarray = NULL, diff.score = T, B = 30, verbose = T, threads = 4)
{
  n = nrow(dat)
  p = ncol(dat)
  II = diag(p)
  w.upper = which(upper.tri(II))
  w.mat = which(upper.tri(II), arr.ind = T)
  mdat = apply(dat, 2, mean)
  sdat = apply(dat, 2, sd)
  stddat = t((t(dat) - mdat)/sdat)
  S = cov(stddat)
  if (is.null(rhoarray)) {
    if (n > p) {
      rhoarray = exp(seq(log(0.001), log(1), length = 100))
    }
    else {
      rhoarray = exp(seq(log(0.1), log(3), length = 100))
    }
  }

  cl <- parallel::makeCluster(mc <- getOption("cl.cores", threads))
  parallel::clusterExport(cl, varlist=c("rhoarray", "S", "extendedBIC", "dat"), envir=environment())
  parallel::clusterEvalQ(cl, library(glasso))
  BIC <- unlist(parLapply(cl, rhoarray, function(rhoarrayi, S, dat){
    fit.gl1 = glasso(S, rho = rhoarrayi)
    fit.gl2 = glasso(S, rho = rhoarrayi, w.init = fit.gl1$w,
      wi.init = fit.gl1$wi)
    extendedBIC(gamma = 0, omegahat = fit.gl2$wi, S = S, n = nrow(dat))
  }, S, dat))
  parallel::stopCluster(cl)

  rho = rhoarray[which.min(BIC)]
  fit.gl1 = glasso(S, rho = rho)
  fit.gl2 = glasso(S, rho = rho, w.init = fit.gl1$w, wi.init = fit.gl1$wi)
  Omega = fit.gl2$wi
  diag.Omega = diag(Omega)
  P = -Omega/diag.Omega
  diag(P) = 0
  if (verbose)
    cat("Step 1 of DINGO is finished at", date(), "\n")
  tY.org = stddat %*% (II - t(P))
  mdat = apply(tY.org, 2, mean)
  sdat = apply(tY.org, 2, sd)
  std.tY = t((t(tY.org) - mdat)/sdat)
  fit.g = DINGO::Greg.em(std.tY ~ x)
  if (verbose)
    cat("Step 2 of DINGO is finished at", date(), "\n")
  if (diff.score) {
    if (verbose)
      cat("Bootstrap scoring is started at", date(), "\n")
    boot.fit = fast.scoring.boot(stddat = stddat, z = x, Omega = Omega,
      A = fit.g$A, B = fit.g$B, boot.B = B, verbose = verbose)
    if (verbose)
      cat("Bootstrap scoring is done at", date(), "\n")
    return(list(genepair = boot.fit$genepair, levels.x = boot.fit$levels.z,
      R1 = boot.fit$R1, R2 = boot.fit$R2, boot.diff = boot.fit$boot.diff,
      diff.score = boot.fit$diff.score, rho = rho, P = P,
      Q = fit.g$B, Psi = fit.g$A))
  }
  else {
    genepair = data.frame(gene1 = colnames(dat)[w.mat[, 1]],
      gene2 = colnames(stddat)[w.mat[, 2]])
    levels.x = unique(x)
    R1 = -scaledMat(solve(Sigmax(Q = fit.g$B, P = P, Psi = fit.g$A,
      x = c(1, levels.x[1]))))[w.upper]
    R2 = -scaledMat(solve(Sigmax(Q = fit.g$B, P = P, Psi = fit.g$A,
      x = c(1, levels.x[2]))))[w.upper]
    return(list(genepair = genepair, levels.x = levels.x,
      R1 = R1, R2 = R2, boot.diff = NULL, diff.score = NULL,
      rho = rho, P = P, Q = fit.g$B, Psi = fit.g$A))
  }
}

#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
fast.scoring.boot = function (stddat, z, Omega, A, B, boot.B = 100, verbose = T)
{
  p = ncol(stddat)
  n = nrow(stddat)
  w.upper = which(upper.tri(Omega))
  w.mat = which(upper.tri(Omega), arr.ind = T)
  II = diag(p)
  levels.z = unique(z)
  stopifnot(length(levels.z) == 2)
  diag.Omega = diag(Omega)
  P = -Omega/diag.Omega
  diag(P) = 0
  tY.org = stddat %*% (II - t(P))

  cl <- makeCluster(threads, type = "SOCK")
  clusterExport(cl, c("n", "tY.org", "z", "levels.z", "P", "w.upper", "calculateDiff"))
  clusterEvalQ(cl, library(DINGO))
  bootList <- parLapply(cl, 1:boot.B, function(x){
    tryCatch(calculateDiff(n, tY.org, z, levels.z, P, w.upper), error = function(e) NA)
  })
  stopCluster(cl)
  boot.diff <- t(na.omit(do.call(rbind, bootList)))

  genepair = data.frame(gene1 = colnames(stddat)[w.mat[, 1]],
    gene2 = colnames(stddat)[w.mat[, 2]])
  R1 = -scaledMat(solve(Sigmax(Q = B, P = P, Psi = A, x = c(1,
    levels.z[1]))))[w.upper]
  R2 = -scaledMat(solve(Sigmax(Q = B, P = P, Psi = A, x = c(1,
    levels.z[2]))))[w.upper]
  diff.score = (trans.Fisher(R1) - trans.Fisher(R2))/apply(boot.diff,
    1, sd)
  return(list(genepair = genepair, levels.z = levels.z, R1 = R1,
    R2 = R2, boot.diff = boot.diff, diff.score = diff.score))
}

#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
calculateDiff = function(n, tY.org, z, levels.z, P, w.upper){
  w.id = sample(1:n, replace = T)
  tY = tY.org[w.id, ]
  mdat = apply(tY, 2, mean)
  sdat = apply(tY, 2, sd)
  std.tY = t((t(tY) - mdat)/sdat)
  fit.g = DINGO::Greg.em(std.tY ~ z[w.id])
  smat = diag(sdat)
  sigmaX1 = smat %*% Sigmax(Q = fit.g$B, P = P, Psi = fit.g$A,
    x = c(1, levels.z[1])) %*% smat
  omegaX1 = solve(sigmaX1)
  boot.RX1 = DINGO::trans.Fisher(-scaledMat(omegaX1)[w.upper])
  sigmaX2 = smat %*% Sigmax(Q = fit.g$B, P = P, Psi = fit.g$A,
    x = c(1, levels.z[2])) %*% smat
  omegaX2 = solve(sigmaX2)
  boot.RX2 = DINGO::trans.Fisher(-scaledMat(omegaX2)[w.upper])
  return(DINGO::trans.Fisher(-scaledMat(omegaX1)[w.upper]) -
      trans.Fisher(-scaledMat(omegaX2)[w.upper]))
}
