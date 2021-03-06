% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class.R, R/convert.R
\name{as.icd_long_data}
\alias{as.icd_long_data}
\alias{as.icd_wide_data}
\alias{icd_long_data}
\alias{icd_wide_data}
\alias{is.icd_long_data}
\alias{is.icd_wide_data}
\alias{icd_long_data}
\title{Convert between and identify 'long' and 'wide' patient data formats}
\usage{
as.icd_long_data(x)

as.icd_wide_data(x)

icd_long_data(...)

icd_wide_data(...)

is.icd_long_data(x)

is.icd_wide_data(x)
}
\arguments{
\item{x}{\code{data.frame} or \code{matrix} to set class, or convert.}

\item{...}{arguments passed on to create a \code{data.frame}}
}
\description{
Long and Wide Formats: As is common with many datasets, key
  variables can be concentrated in one column or spread over several. Tools
  format of clinical and administrative hospital data, we can perform the
  conversion efficiently and accurately, while keeping some metadata about
  the codes intact, e.g. whether they are ICD-9 or ICD-10.

\code{icd_long_data} and \code{icd_wide_data} create
  \code{data.frame}s using all the arguments, and sets the class, whereas
  \code{as.icd_long_data} and \code{as.icd_wide_data} just set the class of
  an existing \code{data.frame}.
}
\details{
Long or wide format ICD data is expected to be in a data frame. It
  does not carry any ICD classes at the top level, even if it only contains
  one type of code, but its constituent columns may have a class specified,
  e.g. 'icd9'.
}
\section{Functions}{
\itemize{
\item \code{as.icd_long_data}: Set class on a matrix or data.frame to
\code{icd_long_data}. To convert wide to long data, use
\code{\link{wide_to_long}}.

\item \code{as.icd_wide_data}: Set class on a matrix or data.frame to
\code{icd_wide_data}. To convert long to wide data, use
\code{\link{long_to_wide}}.

\item \code{icd_long_data}: Construct a \code{data.frame}, adding the \code{icd_long_data} class.

\item \code{icd_wide_data}: Construct a \code{data.frame}, adding the \code{icd_wide_data} class.

\item \code{is.icd_long_data}: Return \code{TRUE} if \code{x} has the \code{icd_long_data} class.

\item \code{is.icd_wide_data}: Return \code{TRUE} if \code{x} has the \code{icd_wide_data} class.
}}

\seealso{
Other ICD code conversion: \code{\link{convert}},
  \code{\link{long_to_wide}}, \code{\link{wide_to_long}}

Other ICD code conversion: \code{\link{convert}},
  \code{\link{long_to_wide}}, \code{\link{wide_to_long}}
}
