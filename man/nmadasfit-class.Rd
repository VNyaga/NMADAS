% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{class}
\name{nmadasfit-class}
\alias{nmadasfit-class}
\title{Class nmadasfit}
\description{
A nmadasfit class in the NMADAS package.
}
\section{Slots}{

\describe{
\item{\code{data}}{data A data-frame with no missing values containg TP, TN, FP, FN, study and test names.}

\item{\code{comparator}}{Name of comparator test.}

\item{\code{fit}}{an object of class stanfit returned by the function sampling.}

\item{\code{labels}}{A character vector witht the names of the tests evaluated.}

\item{\code{S.ID}}{Study identifier}

\item{\code{T.ID}}{Test identifier}
}}

\seealso{
\link{fit}

Other nmadas: \code{\link{nmadasdata-class}},
  \code{\link{nmadasmodel-class}}
}
\author{
Victoria N Nyaga \email{victoria.nyaga@outlook.com}
}
