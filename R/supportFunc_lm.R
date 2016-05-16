#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
lm_singlePredictor = function(x, y, xlab, ylab, main, lim = NULL) {
    fit <- lm(y ~ x)
    est <- round(coef(summary(fit))["x", "Estimate"], 3)
    pval <- round(coef(summary(fit))["x", "Pr(>|t|)"], 3)
    if (is.null(lim)) {
        plot(y ~ x, xlab = xlab, ylab = ylab, main = main)
        abline(fit)
    } else {
        plot(y ~ x, xlab = xlab, ylab = ylab, main = main, xlim = lim$xlim,
            ylim = lim$ylim)
        abline(fit)
        legend(x = lim$xlim[1], y = lim$ylim[2], c(paste("slope",
            est, sep = "="), paste("p-value", pval, sep = "=")),
            bty = "n")
        text(x = x, y = y, labels = names(x))
    }
}

#' table of classification performances
#'
#' takes in predited weights and true labels and determines performance characterisitcs
#' @param weights are the predicted scores/probablities of test data
#' @param trubeLabels are the true labels associated with the test data
#' @param direction = "auto", ">", "<"
#' @export
lme_interactionBinaryCont = function(x, y, binary, replicates,
    xlab, ylab, main, lim = NULL) {
    fit.lme <- summary(lme(y ~ x * binary, random = ~1 | replicates))$tTable
    est.lme <- round(fit.lme["x:binaryDEP", "Value"], 3)
    pval.lme <- round(fit.lme["x:binaryDEP", "p-value"], 3)

    if (is.null(lim)) {
        plot(y[binary == levels(binary)[1]] ~ x[binary == levels(binary)[1]],
            xlab = xlab, ylab = ylab, main = main)
        abline(lm(y[binary == levels(binary)[1]] ~ x[binary == levels(binary)[1]]))
        points(y[binary == levels(binary)[2]] ~ x[binary == levels(binary)[2]])
        abline(lm(y[binary == levels(binary)[2]] ~ x[binary == levels(binary)[2]]))
    } else {
        plot(y[binary == levels(binary)[1]] ~ x[binary == levels(binary)[1]],
            xlab = xlab, ylab = ylab, main = main, col = 1, pch = 19,
            xlim = lim$xlim, ylim = lim$ylim)
        abline(lm(y[binary == levels(binary)[1]] ~ x[binary == levels(binary)[1]]),
            col = 1)
        points(y[binary == levels(binary)[2]] ~ x[binary == levels(binary)[2]],
            col = 2, pch = 19)
        abline(lm(y[binary == levels(binary)[2]] ~ x[binary == levels(binary)[2]]),
            col = 2)
        legend(x = lim$xlim[1], y = lim$ylim[2], c(paste("slope",
            est.lme, sep = "="), paste("p-value", pval.lme, sep = "=")),
            bty = "n")
        legend("topright", levels(binary), col = 1:2, pch = 19)
    }
}

