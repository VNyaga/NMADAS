#' @import methods


#' @rdname fit
#' @export
#' @description Fit a NMA model to the data.
#' @param object A model written in the stan format from \link{nmadasmodel}. If the model is not specified,
#' a hierachical beta-binomial model with frank copula is fitted.
#' @param data A data-frame with no missing values containg TP, TN, FP, FN, SID and TID.
#' @param S.ID A string indicating the name of the column with the study identifier.
#' @param T.ID A string indicating the name of the column with the test identifier.
#' @param tp A string indicating the name of the columnt with the true positives.
#' @param fn A string indicating the name of the columnt with the false negatives.
#' @param tn A string indicating the name of the columnt with the true negatives.
#' @param fp A string indicating the name of the columnt with the false positives.
#' @param Comparator The name of the comparator test when relative sensitivity and specificity are required. By default the first test as
#' arranged alphabetically is the comparator.
#' @param chains A positive numeric value specifying the number of chains, default is 3.
#' @param iter A positive numeric value specifying the number of iterations per chain. The default is 6000.
#' @param warmup A positive numeric value (<iter) specifying the number of iterations to be discarded(burn-in/warm-up). The default is 1000.
#' @param thin A positive numeric value specifying the interval in which the samples are stored. The default is 10.
#' @param cores A positive numeric values specifying the number of cores to use to execute parallel sampling. When the hardware has more at least 4 cores,
#' the default is 3 cores and otherwise 1 core.
#' @param ... Other optional parameters as specified in \link[rstan]{stan}.
#' @return An object of nmadasfit class.

setGeneric(name="fit",
          function(object, ...){ standardGeneric("fit") })


#' A function to fit the model.
#' @rdname fit
#' @method fit nmadasmodel
#' @export
setMethod("fit",
          signature = "nmadasmodel",
          function(object,
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
              fit.nmadasmodel(nma.model = object,
                            data = data,
                            S.ID = S.ID,
                            T.ID = T.ID,
                            Comparator = Comparator,
                            tp = tp,
                            fn = fn,
                            tn = tn,
                            fp = fp,
                            cores = cores,
                            chains = chains,
                            iter =iter,
                            warmup = warmup,
                            thin = thin,
                            ...)
          })

#===============================================================================================================
#' A function to print the model.
#' @rdname show
#' @param object A nmadasmodel object returned by \link{nmadasmodel} function.
#' @method show nmadasmodel
#' @keywords internal
setMethod("show", signature = "nmadasmodel",
          function(object){
              cat(text = object@model)
          })



#' A function to print the results.
#' @rdname show
#' @method show nmadasfit
#' @keywords internal
setMethod("show", signature = "nmadasfit",
          function(object){
              print.nmadasfit(object)
})

#' A function to produce traceplots.
#' @param object A cdtafit object from \link{fit}
#' @param ... Extra optional arguments as defined in \link[rstan]{stan_trace}.
#' @rdname traceplot
#' @export

#===============================================================================================================
# setGeneric(name="traceplot",
#            def=function(object,
#                         pars=c('MU'),
#                         colourset = "mix-red-brightblue",
#                         ...){
#              standardGeneric("traceplot")
#              })

#' A function to produce traceplots.
#' @rdname  traceplot
#' @method traceplot cdtafit
#' @export
#'
setMethod("traceplot", signature = "nmadasfit",
          function(object,
                   pars=c('MU'),
                   colourset = "mix-red-brightblue",
                   ...){
              traceplot.nmadasfit(x=object,
                                  pars=pars,
                                  colourset = colourset,
                                  ...)

})


#' @rdname  traceplot
#' @param object A nmadasfit object from \link{fit}.
#' @export

#===============================================================================================================

#' @rdname Forestplot
#' @description  Produce forest plots for categorical covariates.
#' @param object A nmadasfit object from \link{fit}.
#' @param pointcolour A text indicating the colour of the study specific points. Default is "grey70".
#' @param pointsize Size of the study specific points. Default is 2.
#' @param RR Logical which is by default TRUE to draw a forest plot of the relative sensitivity and relative specificity.
#' @param textsize Size of the texts.
#' @param dodgewidth An optional numeric value to adjust the dogding position. The default is 1. See \link[ggplot2]{position_dodge}.
#' @param dp An optional positive value to control the number of digits to print when printing numeric values. The default is 2.
#' @param vlinecolour colour of the line at RR = 1. By default it is "blue".
#' @param textlabel The text that appear below the plots. By default it is "Mean [95\% CI]".
#' @param ... other \link[rstan]{stan} options.
#' @return A forestplot by ggplot2.

# setGeneric(name="forestplot", function(object,
#                                        vlinecolour = "blue",
#                                        textsize = 4,
#                                        pointcolour = "grey70",
#                                        pointsize = 2,
#                                        dp = 2,
#                                        textlabel ="Mean [95% CI]",
#                                        dodgewidth = 1,
#                                        RR = TRUE,
#                                        ...){
#   standardGeneric("forestplot")
#   })

#' A function to produce forest plots.
#' @rdname  forestplot
#' @method forestplot nmadasfit
#' @export
#'
setMethod("plot", signature = "nmadasfit",
          function(object,
                   vlinecolour = "blue",
                   textsize = 4,
                   pointcolour = "grey70",
                   pointsize = 2,
                   dp = 2,
                   textlabel ="Mean [95% CI]",
                   dodgewidth = 1,
                   RR = TRUE,
                   ...){
              forestplot.nmadasfit(object = object,
                                   vlinecolour = vlinecolour,
                                   textsize = textsize,
                                   pointcolour = pointcolour,
                                   pointsize = pointsize,
                                   dp = dp ,
                                   textlabel = textlabel,
                                   dodgewidth = dodgewidth,
                                   RR = RR ,
                                   ...)

})

#===============================================================================================================
#' #' @description  Draw a network plot.
#' #' @param data A data.frame with atleast two variables i.e study and test identifier.
#' #' @param SID A text indicating the name of the variable that identifies the different studies.
#' #' @param TID A text indicating the name of the variable that identifies the different tests.
#' #' @param ... additional options. See \link[pcnetmeta]{nma.networkplot} for more details.
#'
#' setGeneric(name="networkplot", function(object,
#'                                         S.ID,
#'                                         T.ID,
#'                                         alphabetic = FALSE,
#'                                         edge.col = "grey",
#'                                         node.col = "orange",
#'                                         ...){
#'   standardGeneric("networkplot")
#' })
#'
#' #' A function for plotting a network.
#' #' @rdname  networkplot
#' #' @method networkplot nmadasdata
#' #' @export

