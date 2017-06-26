#' summary
#' @description Generate a summary of the results.
#' @return The posterior mean and 95 percent credible intervals, n_eff, Rhat and WAIC.
#' @param object An object from \link{fit}.
#' @param digits An optional positive value to control the number of digits to print when printing numeric values.
#' @param ... other \link[rstan]{stan} options.
#' @examples
#'
#' \dontrun{
#'
#' fit1 <- fit(data=telomerase,
#'              SID = "ID",
#'              copula="fgm",
#'              iter = 400,
#'              warmup = 100,
#'              seed=1,
#'              cores=1)
#'
#' ss <- summary(fit1)
#'
#' }
#' @references {Watanabe S (2010). Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in Singular
#' Learning Theory. Journal of Machine Learning Research, 11, 3571-3594.}
#' @references {Vehtari A, Gelman A (2014). WAIC and Cross-validation in Stan. Unpublished, pp. 1-14.}
#' @export
#' @author Victoria N Nyaga

summary.nmadasfit <- function(object,
                            RR = TRUE,
                            SIndex = TRUE,
                            digits=3,
                            ...){


#=======================Extract Model Parameters ===================================#
   sm <- rstan::summary(object@fit, ...)

   #Obtain the summaries
   obtainsummary <- function(par) {
     x <- data.frame(summary(object@fit, pars=par)$summary[, c("mean", "2.5%", "50%", "97.5%", "n_eff", "Rhat")])
     names(x) <- c("Mean", "Lower", "Median", "Upper", "n_eff", "Rhat")

      if (par != "S"){

        if (RR){
          param <- c("RR.Sens", "RR.Spec")
        }
        if (par == "MU") {
          param <- c("Sensitivity", "Specificity")
        }

        x$Parameter <- rep(param, each=nrow(x)/2)
        x$Test <- rep(object@labels, 2)
        x <- x[, c("Test", "Parameter", "Mean", "Lower", "Median", "Upper", "n_eff", "Rhat")]
        x <- x[order(x$Test),]
      }
     else{
       x$Test <- object@labels
       x <- x[, c("Test", "Mean", "Lower", "Median", "Upper", "n_eff", "Rhat")]
     }

     row.names(x) <- NULL
     x
   }

   MU <- obtainsummary("MU")

   if (RR){
     RR <- obtainsummary("RR")
   }

   if (SIndex) {
     S <- obtainsummary("S")
   }

    w <- waic(object@fit)

    out <- list(MU=MU,
                RR = RR,
                S = S,
                WAIC=w,
                allsm=sm)

	out
}

