\docType{data}
\name{ppi_PC_2012}
\alias{ppi_PC_2012}
\alias{ppi_PC_2012.m}
\title{PPI network from Pathway Commons v1}
\format{Adjacency matrix for the PPI network  
}

\usage{
data(ppi_PC_2012)
}
\description{
Adjacency matrix describing a curated PPI network derived from Pathway
Commons on the 13th Jun 2012, containing over 8000 proteins. Proteins
have been annotated to the genes (Entrez gene ID) they encode. Diagonal
entries are 0. Other entries are zero or 1 if an interaction has been
reported. Please see references for detailed construction of the PPI network.
}
\details{
Upon loading the data object, we find the following:
\itemize{
  \item \code{ppi_PC_2012.m}: the adjacency matrix
}
}
\references{
Teschendorff AE, Enver T. Single-cell entropy for accurate estimation of
differentiation potency from a cell’s transcriptome. Nat.Commun.2017
}
\keyword{datasets}
