#' prepdata
#' @description Prepare the data
#' @param data A data-frame with no missing values containg TP, TN, FP, FN, study and test names.
#' @param S.ID A string indicating the name of the column with the study identifier.
#' @param T.ID A string indicating the name of the column with the test identifier.
#' @param tp A string indicating the name of the column with the true positives. Default is TP.
#' @param fn A string indicating the name of the column with the false negatives. Defautl is FN.
#' @param tn A string indicating the name of the column with the true negatives. Default is TN.
#' @param fp A string indicating the name of the column with the false positives. Default is FP.
#' @author Victoria N Nyaga
#' @keywords internal
prepdata.nmadas <- function(data,
    S.ID,
    T.ID,
    tp = NULL,
    fn = NULL,
    tn = NULL,
    fp = NULL){

    data$SID <- data[ , grepl(S.ID, names(data), ignore.case = TRUE)]
    data$SID <- as.numeric(factor(data$SID));

    data$TID <- data[ , grepl(T.ID, names(data), ignore.case = TRUE)]

    tests <- levels(factor(data$TID));
    data$TID <- as.numeric(factor(data$TID))

    if (!is.null(tp)) {
      if (tp != "TP") data$TP <- data[, grepl(tp, names(data), ignore.case = TRUE)]

      if (fn != "FN") data$FN <- data[, grepl(fn, names(data), ignore.case = TRUE)]

      if (tn != "TN") data$TN <- data[, grepl(tn, names(data), ignore.case = TRUE)]

      if (fp != "FP") data$FP <- data[, grepl(fp, names(data), ignore.case = TRUE)]

      data$Dis <- data$TP + data$FN
      data$NDis <- data$FP + data$TN
    }

    out <- new("nmadasdata",
               data = data,
               S.ID = S.ID,
               T.ID = T.ID,
               labels = tests)

    out
}


#'  Compute log pointwise predictive density, effective number of parameters and WAIC.
#'
#' @param model stanfit object
#'
#' @author Victoria N Nyaga
#' @keywords internal
#============================== WAIC =====================================================#

waic <- function (model){
    log_sum_exp <- function(x) {
        x_max <- base::max(x)
        x_max + log(base::sum(exp(x - x_max)))
    }
	#Calculate posterior variances from simulation
	colVars <- function (a){
		diff <- a - base::matrix (base::colMeans(a), nrow(a), ncol(a), byrow=TRUE)
		vars <- base::colMeans(diff^2)*base::nrow(a)/(base::nrow(a)-1)
		return (vars)
	}

    log_lik <- rstan::extract(model, "loglik")$loglik
    summands <- apply(log_lik, 2, log_sum_exp)
    correc <- - ncol(log_lik) * log(nrow(log_lik))
    lppd <- sum(summands) + correc
    p_waic_1 <- 2 * (sum(summands - base::colMeans(log_lik)) + correc)
    p_waic_2 <- sum (colVars(log_lik))
    waic_2 <- -2*lppd + 2*p_waic_2
    return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd,
                  p_waic_1=p_waic_1))
}

#'@examples
#'
#' data(demodata)
#'
#' df <- prepdata(data = demodata,
#'                 S.ID = "study",
#'                 T.ID = "Test")
#'
#' class(df)=="nmadasdata"


