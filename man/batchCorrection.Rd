% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiBaC.R
\name{batchCorrection}
\alias{batchCorrection}
\title{batchCorrection}
\usage{
batchCorrection(mbac, multiBatchDesign, Interaction = FALSE,
  Variability = 0.9)
}
\arguments{
\item{mbac}{mbac object generated by createMbac. PLS models slot must be present.}

\item{multiBatchDesign}{A list containing the original and predicted omic for each batch. All omics must be present in every batch. Output object of genMissingOmics function}

\item{Interaction}{Logical. Whether to model the interaction between experimental factors and bacth factor in ARSyN models. By default, FALSE.}

\item{Variability}{From 0 to 1. Minimum percent of data variability that must be explained for each ARSyN model. By default, 0.90.}
}
\value{
Custom mbac object. Elements in a mbac object:
\enumerate{
    \item ListOfBatches: A list of MultiAssayExperiment objects (one per batch).
    \item commonOmic: Name of the common omic between the batches.
    \item CorrectedData: Same structure than ListOfBatches but with the corrected data instead of the original.
    \item PLSmodels: PLS models created during MultiBaC method performance (one model per non-common omic data type).
    \item ARSyNmodels: ARSyN models created during MultiBaC performance (one per omic data type).
    \item InnerRelation: Table of class data.frame containing the inner correlation (i.e. correlation between the scores of X (t) and Y (u) matrices) for each PLS model across all components.
}
}
\description{
This function performs the ARSyNbac correction [1] for each omic contained in mulBatchDesign input object.
}
\examples{
data('multiyeast')

my_mbac <- createMbac (inputOmics = list(A.rna, A.gro, B.rna, B.ribo, C.rna, C.par),
                       batchFactor = c("A", "A", "B", "B", "C", "C"),
                       experimentalDesign = list("A" =  c("Glu+", "Glu+", "Glu+",
                       "Glu-", "Glu-", "Glu-"),
                       "B" = c("Glu+", "Glu+", "Glu-", "Glu-"),
                       "C" = c("Glu+", "Glu+", "Glu-", "Glu-")),
                       omicNames = c("RNA", "GRO", "RNA", "RIBO", "RNA", "PAR"),
                       commonOmic = "RNA")

my_mbac_2 <- genModelList (my_mbac, test.comp = NULL,
                           scale = FALSE, center = TRUE,
                           crossval = NULL,
                           showinfo = TRUE)
multiBatchDesign <- genMissingOmics(my_mbac_2)
my_finalwise_mbac <- batchCorrection(my_mbac_2,
                                     multiBatchDesign = multiBatchDesign,
                                     Interaction = FALSE,
                                     Variability = 0.9)

}
\references{
[1] Nueda MJ, Ferrer A, Conesa A. ARSyN: A method for the identification and removal of systematic noise in multifactorial time course microarray experiments. Biostatistics. 2012;13(3):553-566. doi:10.1093/biostatistics/kxr042
}
