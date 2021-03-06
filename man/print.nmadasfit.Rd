% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.nmadasfit}
\alias{print.nmadasfit}
\title{Print a summary of the fitted model.}
\usage{
\method{print}{nmadasfit}(x, digits = 3, ...)
}
\arguments{
\item{x}{An nmadasfit object from \link{fit}.}

\item{digits}{An optional positive value to control the number of digits to print when printing numeric values. The default is 3.}

\item{...}{other \link[rstan]{stan} options.}
}
\value{
The posterior mean and 95 percent credible intervals, n_eff, Rhat and WAIC.
}
\description{
Print a summary of the fitted model.
}
\examples{

\dontrun{

data(demodata)

fit1 <- fit(S.ID='study',
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

print(fit1)

}
}
\references{
{Watanabe S (2010). Asymptotic Equivalence of Bayes Cross Validation and Widely Applicable Information Criterion in Singular
Learning Theory. Journal of Machine Learning Research, 11, 3571-3594.}

{Vehtari A, Gelman A (2014). WAIC and Cross-validation in Stan. Unpublished, pp. 1-14.}
}
\author{
Victoria N Nyaga
}
