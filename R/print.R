#' Print a summary of the fitted model.

#' @return The posterior mean and 95 percent credible intervals, n_eff, Rhat and WAIC.
#' @param x An nmadasfit object from \link{fit}.
#' @param digits An optional positive value to control the number of digits to print when printing numeric values. The default is 3.
#' @param ... other \link[rstan]{stan} options.
#' @examples
#'
#' \dontrun{
#'
#' data(demodata)
#'
#' fit1 <- fit(S.ID='study',
#'			   T.ID = 'Test',
#'			   tp = 'TP',
#'			   tn = 'TN',
#'			   fp = 'FP',
#'			   fn = 'FN',
#'             data = demodata,
#'             iter = 6000,
#'             warmup = 2000,
#'             thin = 5,
#'             seed = 3)
#'
#' print(fit1)
#'
#'}
#' @references {Watanabe S (2010). Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in Singular
#' Learning Theory. Journal of Machine Learning Research, 11, 3571-3594.}
#' @references {Vehtari A, Gelman A (2014). WAIC and Cross-validation in Stan. Unpublished, pp. 1-14.}
#' @export
#' @author Victoria N Nyaga
print.nmadasfit <- function(x, digits=3, ...){

#=======================Extract Model Parameters ===================================#
   sm <- summary.nmadasfit(x, ...)

# =============================  Print Summaries =====================================


    cat("Sampling Features\n\n")
    cat(paste("Sampling algorithm: ", attr(x@fit@sim$samples[[1]], "args")$sampler_t, "\n", sep=""))

    cat(paste("\n", x@fit@sim$chains, " chain(s)", "each with iter=", x@fit@sim$iter,"; ", "warm-up=",
              x@fit@sim$warmup, "; ", "thin=", x@fit@sim$thin, ".\n",  sep=""))
    cat(paste("post-warmup draws per chain=",
              (x@fit@sim$iter-x@fit@sim$warmup)/x@fit@sim$thin, ";", "total post-warmup draws=",
              ((x@fit@sim$iter-x@fit@sim$warmup)/x@fit@sim$thin)*x@fit@sim$chains, ".\n", sep=""))

    w <- waic(x@fit)

    cat("\nPredictive accuracy of the model\n\n")
    cat(paste("Log point-wise predictive density (LPPD): ", sprintf(paste("%.", digits, "f", sep=''),w$lppd), sep=''))
    cat("\n")
    cat(paste("Effective number of parameters: ",  sprintf(paste("%.", digits, "f", sep=''),w$p_waic), sep=''))
    cat("\n")
    cat(paste("Watanabe-Akaike information Criterion (WAIC): ", sprintf(paste("%.", digits, "f", sep=''),w$waic), sep=''))
    cat("\n\n")

    cat("Posterior marginal mean sensitivity and specificity\n\twith 95% credible intervals\n\n")
    print(sm$MU[,c(1:4, 6)], digits=digits, row.names=FALSE)
    cat("\n\n")

    if (is.data.frame(sm$RR)) {
      cat("Posterior marginal relative sensitivity and specificity\n\twith 95% credible intervals\n\n")
      print(sm$RR[,c(1:4, 6)], digits=digits, row.names=FALSE)
      cat("\n\n")
    }

    if (is.data.frame(sm$S)) {
      cat("Superiority index with 95% credible intervals\n\n")
      print(sm$S[order(-sm$S$Mean), c(1:3, 5)], digits=digits, row.names=FALSE)
      cat("\n\n")
    }
}

