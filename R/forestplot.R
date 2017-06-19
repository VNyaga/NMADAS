#' Forest plot
#' @description  Produce forest plots for categorical covariates.
#' @param object A nmadasfit object from \link{fit}.
#' @param pointcolour A text indicating the colour of the study specific points. Default is "grey70".
#' @param pointsize Size of the study specific points. Default is 2.
#' @param RR Logical which is by default TRUE to draw a forest plot of the relative sensitivity and relative specificity.
#' @param vlinecolor A text indication the colour of the line at RR = 1. Default is "blue".
#' @param textsize Size of the texts.
#' @param dodgewidth An optional numeric value to adjust the dogding position. The default is 1. See \link[ggplot2]{position_dodge}.
#' @param dp An optional positive value to control the number of digits to print when printing numeric values. The default is 2.
#' @param ... other \link[rstan]{stan} options.
#' @return A forestplot by ggplot2.

#' @examples
#' \dontrun{
#' data(demodata)
#'
#' fit1 <- fit(S.ID='study',
#'			   T.ID = 'Test',
#'			   tp = 'TP',
#'			   tn = 'TN',
#'			   fp = 'FP',
#'			   fn = 'FN',
#'         data = demodata,
#'         iter = 6000,
#'         warmup = 2000,
#'         thin = 5,
#'         seed = 3)
#'
#'
#' forestplot(fit1)
#' }
#' @references {Watanabe S (2010). Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in Singular
#' Learning Theory. Journal of Machine Learning Research, 11, 3571-3594.}
#' @references {Vehtari A, Gelman A (2014). WAIC and Cross-validation in Stan. Unpublished, pp. 1-14.}
#' @export
#' @author Victoria N Nyaga <victoria.nyaga@outlook.com>
forestplot.nmadas <- function(object,
                              vlinecolour = "blue",
                              textsize = 4,
                              pointcolour = "grey70",
                              pointsize = 2,
                              dp = 2,
                              textlabel ="Mean [95% CI]",
                              dodgewidth = 1,
                              RR = TRUE,
                              ...
                              ){
  ##########################################################################################
  object@data$Sensitivity <- object@data$TP/object@data$Dis
  object@data$Specificity <- object@data$TN/object@data$NDis

  longdata <- reshape2::melt(object@data[, c("SID", "TID", "Sensitivity", "Specificity")],
                             id.vars = c("SID", "TID"))
  names(longdata)[2:4] <- c("Test", "Parameter", "Mean")
  longdata$Test <- factor(longdata$Test)

  sm <- summary.nmadasfit(object, ...)

  dodge <- position_dodge(width = dodgewidth)
  #==================================================================================================

  sm$MU$MU.Mean <- formattable::formattable(sm$MU$Mean, digits = dp, format = "f")
  sm$MU$MU.Lower <- formattable::formattable(sm$MU$Lower, digits = dp, format = "f")
  sm$MU$MU.Upper <- formattable::formattable(sm$MU$Upper, digits = dp, format = "f")
  sm$MU$Text <- paste(sm$MU$MU.Mean, '[', sm$MU$MU.Lower, ', ', sm$MU$MU.Upper, ']', sep='')

  sm$MU$Test <- as.numeric(factor(sm$MU$Test))

  MU.plot <- ggplot(data=sm$MU,
                aes(x = as.factor(Test),
                    y = Mean)) +
    geom_point(data=longdata,
               aes(x = Test,
                   y = Mean),
               size = pointsize,
               shape = 21,
               fill = 'white',
               colour = pointcolour) +
      coord_flip(expand=TRUE) +
    scale_x_discrete(name="",
                     labels=object@labels) +

      facet_grid(~ Parameter) +
      geom_point(position=dodge,
                 size=pointsize,
                 shape=5) +
      geom_errorbar(aes(ymin=Lower,
                        ymax=Upper),
                    width=0,
                    size=0.5,
                    position=dodge) +
    geom_text(aes(x = Test,
                  y = 1.4,
                  label = Text),
              size = textsize,
              colour="black") +
      theme(axis.text.x=element_text(size=13,
                                     colour='black'),
            axis.text.y=element_text(size=13,
                                     colour='black'),
            axis.title.x =element_text(size=13,
                                       colour='black'),
            axis.ticks.y = element_line(size=rel(0),
                                        color='white'),
            strip.text.y = element_text(size = 13,
                                        colour='black',
                                        angle=180),
            strip.text.x = element_text(size = 13,
                                        colour='black'),
            panel.grid.major = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(color = 'black'),
            axis.line.y = element_blank(),
            strip.background = element_blank(),
            legend.position ='right',
            legend.text=element_text(size=13,
                                     colour='black'),
            legend.key = element_rect(fill="white"),
            plot.background = element_rect(fill = "white",
                                           colour='white'),
            panel.spacing = unit(2, "lines")) +
      scale_y_continuous(name="",
                         limits=c(-0.04, 1.6),
                         breaks=c(0, 0.5, 1, 1.4),
                         expand = c(0.005, 0.005),
                         labels = c("0", "0.5", "1", textlabel))


  if (grDevices::dev.interactive()) {
    grDevices::dev.new()
    print(MU.plot)
  }
  #==================================================================================================
  if(RR){
    sm$RR$RR.Mean <- formattable(sm$RR$Mean, digits = dp, format = "f")
    sm$RR$RR.Lower <- formattable(sm$RR$Lower, digits = dp, format = "f")
    sm$RR$RR.Upper <- formattable(sm$RR$Upper, digits = dp, format = "f")
    sm$RR$Text <- paste(sm$RR$RR.Mean, '[', sm$RR$RR.Lower, ', ', sm$RR$RR.Upper, ']', sep='')
    ytext <- max(sm$RR$Upper) + 0.65
    yulim <- max(sm$RR$Upper) + 1 #max RR
    yllim <- min(sm$RR$Lower) - 0.5 #Min RR


    sm$RR$Test <- as.numeric(factor(sm$RR$Test))
    sm$RR$Parameter <- rep(c("Relative sensitivity", "Relative specificity"), length.out=nrow(sm$RR))

    brks <- c(seq(yllim, max(sm$RR$Upper), max(sm$RR$Upper)/3), max(sm$RR$Upper) + 0.65) #Define where axis ticks appear
    labels <- c(as.character(formattable(seq(yllim, max(sm$RR$Upper), max(sm$RR$Upper)/3), digits = dp, format = "f")), textlabel) #Define labels for the axis ticks appear

    RR.plot <- ggplot(data=sm$RR,
                      aes(x = as.factor(Test),
                          y = Mean)) +
        coord_flip(expand=TRUE) +
      scale_x_discrete(name="",
                       labels=object@labels) +

      facet_grid(~ Parameter) +
      geom_point(position=dodge,
                 size=2,
                 shape=5) +
      geom_errorbar(aes(ymin=Lower,
                        ymax=Upper),
                    width=0,
                    size=0.5,
                    position=dodge) +
      geom_text(aes(x=Test,
                    y= ytext,
                    label=Text),
                size=textsize,
                colour="black")+
      geom_text(aes(x = 0.75,
                    y = 0.5,
                    label= "Worse"),
                size=textsize,
                colour="black")+
      geom_text(aes(x = 0.75,
                    y = 1.5,
                    label= "Better"),
                size=textsize,
                colour="black")+
      geom_hline(yintercept = 1,
                 color = vlinecolour,
                 linetype=2) +
      theme(axis.text.x=element_text(size=13,
                                     colour='black'),
            axis.text.y=element_text(size=13,
                                     colour='black'),
            axis.title.x =element_text(size=13,
                                       colour='black'),
            axis.ticks.y = element_line(size=rel(0),
                                        color='white'),
            strip.text.y = element_text(size = 13,
                                        colour='black',
                                        angle=180),
            strip.text.x = element_text(size = 13,
                                        colour='black'),
            panel.grid.major = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(color = 'black'),
            axis.line.y = element_blank(),
            strip.background = element_blank(),
            legend.position ='right',
            legend.text=element_text(size=13, colour='black'),
            legend.key = element_rect(fill="white"),
            plot.background = element_rect(fill = "white",
                                           colour='white'),
            panel.spacing = unit(2, "lines")) +
      scale_y_continuous(name="",
                         limits=c(yllim, yulim),
                         breaks = brks,
                         expand = c(0.005, 0.005),
                         labels = labels)

    if (grDevices::dev.interactive()) {
      grDevices::dev.new()
      print(RR.plot)
    }
  }
}
