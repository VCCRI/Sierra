% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_util.R
\name{ReadPeakCounts}
\alias{ReadPeakCounts}
\title{Read in peak data saved in MEX format}
\usage{
ReadPeakCounts(data.dir = NULL, mm.file = NULL, barcodes.file = NULL,
  sites.file = NULL)
}
\arguments{
\item{data.dir}{directory where output from CountPeaks is stored}

\item{mm.file}{count matrix in MEX format}

\item{barcodes.file}{file containing cell barcodes corresponding to columns in the matrix}

\item{sites.file}{file containing peak coordinate names corresponding to rows in the matrix}
}
\value{
a sparseMatrix
}
\description{
Read in peak data saved in MEX format
}
\examples{
count.mat = ReadPeakCounts(data.dir = data.dir)

}