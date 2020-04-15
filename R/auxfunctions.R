#' @import ggplot2
#' @import graphics
NULL

#' Class MultiBaC
#'
#' @param list List of objects. At least 1 of them.
#' \enumerate{
#'     \item ListOfBatches: A list of MultiAssayExperiment objects (one per batch).
#'     \item commonOmic: Name of the common omic between the batches. It must be one of the names in omicNames argument. If NULL (default), the omic names that appears more times is selected as commonOmic.
#'     \item CorrectedData: Same structure than ListOfBatches but with the corrected data instead of the original.
#'     \item PLSmodels: PLS models created during MultiBaC method performance (one model per non-common omic data type).
#'     \item ARSyNmodels: ARSyN models created during MultiBaC performance (one per omic data type).
#'     \item InnerRelation: Table of class data.frame containing the inner correlation (i.e. correlation between the scores of X (t) and Y (u) matrices) for each PLS model across all components.
#' }
#' @param ... Other attributes
#'
#' @return Instance of class "mbac"
#'
mbacClass = function(list, ...) {
  model = structure(list, class = "mbac",
                    my_attribute = "MultiBaC_output")
  return(model)
}


#' createPLSmodel
#'
#' This function creates a PLS model between two omics data matrices. It also performs an optimization of the
#' number of component of the model that maximize the Q^2 as a measure of prediction performance.
#'
#' @param test.comp Maximum number of components allowed in PLS models. If NULL (default), the minimal effective rank of the matrices is used as the maximum number of components.
#' @param scale Logical. Whether X and Y matrices must be scaled. By default, FALSE.
#' @param center Logical. Whether X and Y matrices must be centered. By default, TRUE.
#' @param omicslist A list of length 2. Each slot must contain a data matrix (features x samples).
#' @param messages If TRUE, messages about the algorithm steps will be displayed in console window.
#' @param crossval Integer indicating the number of cross-validation segments. The number of samples (rows of 'x') must be at least >= crossvalI. If NULL (default), a leave-one-out crossvalidation is conducted.
#' @param regressor Integer (1 or 2): Which of the matrices is the X matrix for PLS model.
#' @param showinfo Logical. Whether to show function process in prompt.
#'
#' @return Instance of class "opls". See package ropls. The final model has the optimal number of components used for prediction.
#'
#' @examples
#' \dontrun{
#' ## Using example data provided by MultiBaC package
#' data('multiyeast')
#'
#' omicsList <- list("RNA" = A.rna, "GRO" = A.gro)
#'
#' ## Create model
#' plsModel <- createPLSmodel(omicsList, crossval = NULL, test.comp = 3,
#' scale = FALSE, center = TRUE, regressor = 1)
#' }
createPLSmodel <- function(omicslist, test.comp, messages = TRUE,
                           scale = FALSE, center = TRUE,
                           crossval, regressor, showinfo = TRUE) {
  # Set preprocessing -----------------------------------------------------------
  if ( center == FALSE & scale == FALSE) {
    pret <- "none"
  } else if (center == TRUE & scale == FALSE) {
    pret <- "center"
  } else {
    pret <- "standard"
  }

  # Set crossval ----------------------------------------------------------------
  if (is.null(crossval)) {
    crossval <- dim(t(omicslist[[1]]))[1]
  }

  # Set test.comp ---------------------------------------------------------------
  max.comp <- min(dim(omicslist[[1]])) - 1
  if ( test.comp > max.comp ) {
    if (showinfo) {
      message(paste0("Warning: (Input test.comp exceeds the minimum dimension of data matrix. test.comp set to ", max.comp, ")"))
    }
    test.comp <- max.comp
  }
  # Create models ---------------------------------------------------------------
  models <- list()
  if (methods::is(regressor, "character")) {
    regressor <- which(names(omicslist) == regressor)
  }
  for ( i in seq_along(names(omicslist))[-regressor]) {

    # COmpute Q2 ------------------------------------------------------------------
    plsModel <- ropls::opls(t(omicslist[[regressor]]), t(omicslist[[i]]),
                            predI = test.comp, fig.pdfC = NULL, info.txtC = NULL,
                            crossvalI = dim(t(omicslist[[1]]))[1], scaleC = pret)
    q2v <- plsModel@modelDF$`Q2(cum)`

    # Built final model -----------------------------------------------------------
    plsModel <- ropls::opls(t(omicslist[[1]]), t(omicslist[[i]]),
                            predI = which(q2v == max(q2v))[1], fig.pdfC = NULL, info.txtC = NULL,
                            crossvalI = dim(t(omicslist[[1]]))[1], scaleC = pret)
    models[[(i-1)]] <- plsModel
  }
  names(models) <- names(omicslist)[-1]
  return(c(models))
}

#' getData
#'
#' @param ListOfBatches A list of MultiAssayExperiment elements. Object returned by inputData function.
#'
#' @return A list of matrices.
#'
#' @examples
#' \dontrun{
#' ## Using example data provided by MultiBaC package
#' data(multiyeast)
#' my_mbac <- createMbac (inputOmics = list(A.rna, A.gro, B.rna, B.ribo, C.rna, C.par),
#'                        batchFactor = c("A", "A", "B", "B", "C", "C"),
#'                        experimentalDesign = list("A" =  c("Glu+", "Glu+",
#'                        "Glu+", "Glu-", "Glu-", "Glu-"),
#'                        "B" = c("Glu+", "Glu+", "Glu-", "Glu-"),
#'                        "C" = c("Glu+", "Glu+", "Glu-", "Glu-")),
#'                        omicNames = c("RNA", "GRO", "RNA", "RIBO", "RNA", "PAR"))
#'
#' omicData <- getData (my_mbac$ListOfBatches)}
#'
getData <- function(ListOfBatches) {
  if ( inherits(ListOfBatches[[1]], "MultiAssayExperiment") ) {

    inputList <- list()

    for ( i in seq_along(names(ListOfBatches)) ) {
      inputList[[names(ListOfBatches)[i]]] <- ListOfBatches[[names(ListOfBatches)[i]]]@ExperimentList
    }


  } else {
    inputList <- ListOfBatches
  }
  return (inputList)
}

