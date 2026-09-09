\docType{data}
\name{dataLiver}
\alias{dataLiver}
\title{scRNA-Seq data from Yang et al}
\format{Normalized and log-transformed scRNA-Seq data representing a
  differentiation timecourse of hepatoblasts into hepatocytes and
  cholangiocytes. Genes have been annotated to unique human Entrez gene
  IDs. There are 447 cells in total and 7 different timepoints
  (embryonic stages).
}

\usage{
load(file.path(Sys.getenv("DIFFR_TEST_DATA_DIR"), "dataLiver.rda"))
}
\description{
The data has been normalized for library-size with cell-dependent
scale-factors, and has been log-transformed with a pseudocount of 1.1.
It is distributed in \file{DiffR_test_data.tar.gz}, not inside the
installed DiffR package.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{scLiverSparse.m}: the sparse scRNA-Seq data matrix.
  \item \code{phenoLiver.v}: a vector of timepoints matching the
  columns of \code{scLiverSparse.m}.
  \item \code{colorLiver.v}: a color vector for the timepoints.
}
}
\references{
Yang L et al. A single-cell transcriptomic analysis reveals precise
pathways and regulatory mechanisms underlying hepatoblast
differentiation. Hepatology 2017.
}
\keyword{datasets}
