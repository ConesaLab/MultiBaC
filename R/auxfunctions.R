#' @import ggplot2
#' @import graphics
NULL

mbacClass = function(list, ...) {
  model = structure(list, class = "mbac",
                    my_attribute = "MultiBaC_output")
  return(model)
}

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
                            predI = test.comp, fig.pdfC = "none", info.txtC = "none",
                            crossvalI = dim(t(omicslist[[1]]))[1], scaleC = pret)
    q2v <- plsModel@modelDF$`Q2(cum)`

    # Built final model -----------------------------------------------------------
    plsModel <- ropls::opls(t(omicslist[[1]]), t(omicslist[[i]]),
                            predI = which(q2v == max(q2v))[1], fig.pdfC = "none", info.txtC = "none",
                            crossvalI = dim(t(omicslist[[1]]))[1], scaleC = pret)
    models[[(i-1)]] <- plsModel
  }
  names(models) <- names(omicslist)[-1]
  return(c(models))
}

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

