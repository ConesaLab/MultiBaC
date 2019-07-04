#' A distributed yeast multi-omic dataset.
#'
#' The whole dataset is composed by six different matrices containing four different omic technologies collected from three different studies.
#' Original data were reduced to 200 features.
#'
#' @format Six matrices (fetarues x samples):
#' \describe{
#'   \item{A.rna}{Gene expression data from lab A}
#'   \item{A.gro}{Genomic-Run-On measurements (transcription rates) from lab A}
#'   \item{B.rna}{Gene expression data from lab B}
#'   \item{B.ribo}{RIBO-seq data (translation rates) from lab B}
#'   \item{C.rna}{Gene expression data from lab C}
#'   \item{C.par}{Global PAR-CLIP data from lab C}
#' }
#' @source Datasets were downloaded from GEO database \url{"https://www.ncbi.nlm.nih.gov/geo/"}
"multiyeast"
