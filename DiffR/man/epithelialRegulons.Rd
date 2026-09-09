\name{epithelialRegulons}
\alias{epithelialRegulons}
\alias{netBREAST.m}
\alias{netCOL.m}
\alias{netESOPH.m}
\alias{netGAST2.m}
\alias{netKID.m}
\alias{netLIV22.m}
\alias{netLUNGm.m}
\alias{netPANC.m}
\alias{netSKIN.m}
\docType{data}
\title{Tissue-Specific Epithelial Regulon Matrices}
\description{
Nine human tissue-specific epithelial regulon matrices used by
\code{SciraEstRegAct}. Rows are Entrez Gene IDs, columns are transcription
factor Entrez Gene IDs, and entries encode signed regulatory associations.
}
\usage{
data(netBREAST.m)
data(netCOL.m)
data(netESOPH.m)
data(netGAST2.m)
data(netKID.m)
data(netLIV22.m)
data(netLUNGm.m)
data(netPANC.m)
data(netSKIN.m)
}
\format{
Each object is a numeric matrix with 18,165 genes. The numbers of
transcription factors are 9 (breast), 17 (colon), 30 (esophagus), 32
(stomach), 14 (kidney), 14 (liver), 8 (lung), 24 (pancreas), and 20 (skin).
}
\details{
The objects are selected automatically from the \code{tissue} argument of
\code{SciraEstRegAct}; users generally do not need to load them directly.
}
\keyword{datasets}
