#' summary.mbac
#'
#' Displays the structure and the content of the object of class mbac.
#'
#' @param object An object of class mbac.
#' @param ... additional arguments affecting the summary produced.
#'
#' @return Custom mbac object structure.
#' @export
#'
#' @rdname summary
#' @examples
#' data('multiyeast')
#'
#' my_mbac <- createMbac (inputOmics = list(A.rna, A.gro, B.rna, B.ribo, C.rna, C.par),
#'                        batchFactor = c("A", "A", "B", "B", "C", "C"),
#'                        experimentalDesign = list("A" =  c("Glu+", "Glu+",
#'                        "Glu+", "Glu-", "Glu-", "Glu-"),
#'                        "B" = c("Glu+", "Glu+", "Glu-", "Glu-"),
#'                        "C" = c("Glu+", "Glu+", "Glu-", "Glu-")),
#'                        omicNames = c("RNA", "GRO", "RNA", "RIBO", "RNA", "PAR"),
#'                        commonOmic = "RNA")
#' summary(my_mbac)
#'
summary.mbac <- function(object, ...) {
  batches <- names(object$ListOfBatches)
  omics <- lapply(object$ListOfBatches, names)
  print(paste0("Object of class mbac: It contains ", length(batches), " different bacthes and ",
               length(unique(unlist(omics))), " omic type(s)."))
  toplot <-  matrix(NA, length(batches), length(unique(unlist(omics))))
  rownames(toplot) <- batches; colnames(toplot) <- unique(unlist(omics))
  for ( i in batches) {
    toplot[i,names(object$ListOfBatches[[i]])] <- TRUE
  }
  knitr::kable(toplot)
}
