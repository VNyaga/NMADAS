#' fit
#' @description Fit a NMA model to the data.
#' @param nma.model A model written in the stan format from \link{nmamodel}. If the model is not specified,
#' a hierachical beta-binomial model with frank copula is fitted.
#' @param data A data-frame with no missing values containg TP, TN, FP, FN, SID and TID.
#' @param S.ID A string indicating the name of the column with the study identifier.
#' @param T.ID A string indicating the name of the column with the test identifier.
#' @param tp A string indicating the name of the columnt with the true positives.
#' @param fn A string indicating the name of the columnt with the false negatives.
#' @param Comparator The name of the comparator test when relative sensitivity and specificity are required. By default the first test as
#' arranged alphabetically is the comparator.
#' @param tn A string indicating the name of the columnt with the true negatives.
#' @param fp A string indicating the name of the columnt with the false positives.
#' @param chains A positive numeric value specifying the number of chains, default is 3.
#' @param iter A positive numeric value specifying the number of iterations per chain. The default is 6000.
#' @param warmup A positive numeric value (<iter) specifying the number of iterations to be discarded(burn-in/warm-up). The default is 1000.
#' @param thin A positive numeric value specifying the interval in which the samples are stored. The default is 10.
#' @param cores A positive numeric values specifying the number of cores to use to execute parallel sampling. When the hardware has more at least 4 cores,
#' the default is 3 cores and otherwise 1 core.
#' @param ... Other optional parameters as specified in \link[rstan]{stan}.
#' @return An object of nmadasfit class.
#'
#' @examples
#' \dontrun{
#' data(demodata)
#'
#' modelcode <- nmadasmodel()
#'
#' fit1 <- fit(nma.model = modelcode
#'         S.ID='study',
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
#' modelcode <- nmadasmodel(copula = "fgm", marginals = "beta")
#'
#' fit2 <- fit(nma.model = modelcode,
#'			   S.ID='study',
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
#' }
#'
#'@references {Agresti A (2002). Categorical Data Analysis. John Wiley & Sons, Inc.}
#'@references {Clayton DG (1978). A model for Association in Bivariate Life Tables and its Application in
#'Epidemiological Studies of Familial Tendency in Chronic Disease Incidence. Biometrika,65(1), 141-151.}
#'@references {Frank MJ (1979). On The Simultaneous Associativity of F(x, y) and x + y - F(x, y). Aequationes Mathematicae, pp. 194-226.}
#'@references {Farlie DGJ (1960). The Performance of Some Correlation Coefficients for a General Bivariate
#'Distribution. Biometrika, 47, 307-323.}
#'@references {Gumbel EJ (1960). Bivariate Exponential Distributions. Journal of the American Statistical Association, 55, 698-707.}
#'@references {Meyer C (2013). The Bivariate Normal Copula. Communications in Statistics - Theory and Methods, 42(13), 2402-2422.}
#'@references {Morgenstern D (1956). Einfache Beispiele Zweidimensionaler Verteilungen. Mitteilungsblatt furMathematische Statistik, 8, 23 - 235.}
#'@references {Sklar A (1959). Fonctions de Repartition a n Dimensions et Leurs Marges. Publications de l'Institut de Statistique de L'Universite de Paris, 8, 229-231.}
#'@export
#'@importFrom rstan sampling
#'@importFrom rstan stan_model
#' @author Victoria N Nyaga <victoria.nyaga@outlook.com>
fit.nmadasmodel <- function(
  nma.model,
  data,
  S.ID,
  T.ID,
  Comparator = 'NA',
  tp = NULL,
  fn = NULL,
  tn = NULL,
  fp = NULL,
  cores = 3,
  chains = 3,
  iter = 6000,
  warmup = 1000,
  thin = 10,
  ...){

  #=================================================================================
  #====================              prepare the data           ====================
  df <- prepdata.nmadas(
      data = data,
      S.ID = S.ID,
      T.ID = T.ID,
      tp = tp,
      fn = fn,
      tn = tn,
      fp = fp)

  if (Comparator != 'NA'){
    tdf <- data.frame(TID = 1:length(data$labels), Test=data$labels)
    CIndex <- df$TID[df$Test==Comparator]
  }
  else  {
    CIndex <- 1
  }

    N <- nrow(df@data)
    Ns <- max(df@data$SID)
    Nt <- max(df@data$TID)

    datalist <- list(
      N = N,
      Ns = Ns,
      Nt = Nt,
      TP = data$TP,
      Dis = data$Dis,
      TN = data$TN,
      NDis = data$NDis,
      Test = data$TID,
      Study = data$SID,
      CIndex = CIndex)

	stanmodel <- rstan::stan_model(model_code = nma.model@model)

  mod <- rstan::sampling(object = stanmodel,
           data=datalist,
           warmup=warmup,
           thin=thin,
           chains=chains,
           cores=cores,
           iter=iter,
           ...)

  out <- new("nmadasfit",
             data = df@data,
             S.ID = S.ID,
             T.ID = T.ID,
             labels = df$labels,
             comparator = Comparator,
             fit = mod)

  out
}

