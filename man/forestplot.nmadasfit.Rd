% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forestplot.R
\name{forestplot.nmadasfit}
\alias{forestplot.nmadasfit}
\title{Forest plot}
\usage{
forestplot.nmadasfit(object, vlinecolour = "blue", textsize = 4,
  pointcolour = "grey70", pointsize = 2, dp = 2,
  textlabel = "Mean [95\% CI]", dodgewidth = 1, RR = TRUE, ...)
}
\arguments{
\item{object}{A nmadasfit object from \link{fit}.}

\item{textsize}{Size of the texts.}

\item{pointcolour}{A text indicating the colour of the study specific points. Default is "grey70".}

\item{pointsize}{Size of the study specific points. Default is 2.}

\item{dp}{An optional positive value to control the number of digits to print when printing numeric values. The default is 2.}

\item{textlabel}{The text that appear below the plots. By default it is "Mean [95\% CI]".}

\item{dodgewidth}{An optional numeric value to adjust the dogding position. The default is 1. See \link[ggplot2]{position_dodge}.}

\item{RR}{Logical which is by default TRUE to draw a forest plot of the relative sensitivity and relative specificity.}

\item{...}{other \link[rstan]{stan} options.}

\item{vlinecolor}{A text indication the colour of the line at RR = 1. Default is "blue".}

\item{vline}{colour of the line at RR = 1. By default it is "blue".}
}
\value{
A forestplot by ggplot2.
}
\description{
Produce forest plots for categorical covariates.
}
\examples{
\dontrun{
data(demodata)

frank <- nmadasmodel()

fit1 <- fit(
        nma.model = frank,
        S.ID='study',
		   T.ID = 'Test',
		   tp = 'TP',
		   tn = 'TN',
		   fp = 'FP',
		   fn = 'FN',
            data = demodata,
            iter = 6000,
            warmup = 2000,
            thin = 5,
            seed = 3)


forestplot(fit1)
}
}
\references{
{Watanabe S (2010). Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in Singular
Learning Theory. Journal of Machine Learning Research, 11, 3571-3594.}

{Vehtari A, Gelman A (2014). WAIC and Cross-validation in Stan. Unpublished, pp. 1-14.}
}
\author{
Victoria N Nyaga <victoria.nyaga@outlook.com>
}
