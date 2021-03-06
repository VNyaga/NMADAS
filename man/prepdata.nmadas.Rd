% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{prepdata.nmadas}
\alias{prepdata.nmadas}
\title{prepdata}
\usage{
prepdata.nmadas(data, S.ID, T.ID, tp = NULL, fn = NULL, tn = NULL,
  fp = NULL)
}
\arguments{
\item{data}{A data-frame with no missing values containg TP, TN, FP, FN, study and test names.}

\item{S.ID}{A string indicating the name of the column with the study identifier.}

\item{T.ID}{A string indicating the name of the column with the test identifier.}

\item{tp}{A string indicating the name of the column with the true positives. Default is TP.}

\item{fn}{A string indicating the name of the column with the false negatives. Defautl is FN.}

\item{tn}{A string indicating the name of the column with the true negatives. Default is TN.}

\item{fp}{A string indicating the name of the column with the false positives. Default is FP.}
}
\description{
Prepare the data
}
\author{
Victoria N Nyaga
}
\keyword{internal}
