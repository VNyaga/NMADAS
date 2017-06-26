#' Traceplot
#' @description Trace plot using ggplot2.
#' @param x An nmadas object from \link{fit}.
#' @param ... additional options. See \link[base]{array} for more details.
#' @return A ggplot trace plot of the parameters of the models mean structure.
#' @examples
#' \dontrun{
#' data(demodata)
#'
#' frank <- nmadasmodel()
#'
#' fit1 <- fit(
#'         nma.model = frank,
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
#' traceplot(fit1)
#' }
#'@export

#' @author Victoria N Nyaga

traceplot.nmadasfit <- function(x, pars=c('MU'), colourset = "mix-red-brightblue", ...){

  draws <- as.array(x@fit, pars = pars, ...)
  bayesplot::color_scheme_set(colourset)
  g <- bayesplot::mcmc_trace(draws)
  if (grDevices::dev.interactive()) grDevices::dev.new()
  print(g)
  }


