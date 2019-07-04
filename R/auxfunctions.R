#' @import stats
#' @import graphics
#' @import ggplot2
#' @include ASCA3f.R
NULL

#' MultiBaC_plot
#'
#' Display two plots summarizing MultiBaC steps.
#'
#' @param ascamodels A list. Each slot contains a vector of cumulative explained variances of the batch estimation.
#' @param q2values A list. Each slot contains a vector of Q^2 values of a PLS model.
#'
#' @return Display two plots.
#'
#' @examples
#' \dontrun{
#' ## Simulating inputs
#' q2values <- list("model1" = c(0.51, 0.76),
#' "model2" = c(0.30, 0.69))
#'
#' ascamodels <- list("A" = c(50, 79),
#' "B" = c(55, 87), "C" = c(70, 90))
#'
#' MultiBaC_plots (q2values, ascamodels)
#' }
MultiBaC_plots <- function(q2values, ascamodels) {
  initpar <- par(c("mfrow"))
  on.exit(par(mfrow = initpar))
  par(mfrow = c(1,2), xpd = TRUE)

  test.comp <- max(unlist(lapply(q2values, length)))

  # Q2 plot -------------------------------------------------------------------

  pallete <- colors()[c(11,17,51,56,29,512,97,653,136,24)]

  # Make plot
  plot(1:3,1:3, type = "n", pch = 19,
       ylim = c(min(unlist(q2values)),1), xaxt = "n",
       xlim = c(1,test.comp), xlab="Number of Components", ylab = "Squared Q value",
       main = "Squared Q plot", bty = "n",
       cex.lab = 1.25, cex.axis = 1.25, font.lab = 2, cex.main=1.5)
  for ( i in seq_along(q2values)) {
    lines(1:length(q2values[[i]]), q2values[[i]], type = "b", pch = 19, col = pallete[i])
  }
  axis(1, seq_len(test.comp), c(seq_len(test.comp)), cex.axis = 1.25)
  legend(test.comp-2, 0.65,
         bty = "n", title = expression(underline(Batches)),
         legend = c(names(q2values)),
         col = c(pallete),
         cex = 1.5, lty = c(1,1), pch = c(19,19))

  # Explained batch-related variability plot ----------------------------------
  plot(c(0,0), type = "n", pch = 19,
       ylim = c(0,100), xaxt = "n",
       xlim = c(0,test.comp), xlab="Number of Components",
       ylab = "Explained batch-related variability (%)",
       main = "ARSyN nº of components", bty = "n",
       cex.lab = 1.25, cex.axis = 1.25, font.lab = 2, cex.main=1.5)
  pallete <- colors()[c(11,17,51,56,29,512,97,653,136,24)]

  for ( i in seq_along(ascamodels)) {
    lines(0:(length(ascamodels[[1]])-1),
          c(ascamodels[[i]]*100),
          type = "b", pch = 19, col = pallete[i])
  }

  axis(1, 0:(length(ascamodels[[1]])-1), 0:(length(ascamodels[[1]])-1), cex.axis = 1.25)
  legend(0.8, 75,
         bty = "n", title = expression(underline(Omics)),
         legend = c("common","non-common: ", names(ascamodels)[-1]),
         col = c(pallete[1], "white", pallete[2:4]),
         cex = 1.5, lty = c(1,0,rep(1,length(ascamodels)-1)), pch = rep(19,length(ascamodels)+1))

  # Advertising about superposition
  if (sum(ascamodels[[1]] == ascamodels[[2]])) {
    message("Caution: Explained variance could be similar for more than two omics. Thus lines and dots could be superpossed")
  }
}

#' createPLSmodel
#'
#' This function creates a PLS model between two omics data matrices. It also performs an optimization of the
#' number of component of the model that maximize the Q^2 as a measure of prediction performance.
#'
#' @param test.comp maximum number of latent variables (or components) to test
#' @param scale TRUE or FALSE. Whether X and Y matrices have to be scaled
#' @param center TRUE or FALSE. Whether X and Y matrices have to be centered
#' @param omicslist A list of length 2. Each slot must contain a data matrix (features x samples).
#' @param messages If TRUE, messages about the algorithm steps will be displayed in console window.
#' @param crossval Integer: number of cross-validation segments. The number of samples (rows of 'x') must be at least >= crossvalI.
#' If NULL (default) leave-one-out crossvalidation is performed
#' @param regressor Integer (1 or 2): Which of the matrices is the X matrix for PLS model.
#'
#' @return Instance of class "opls". See package ropls. The final model has the optimal number of components used for prediction.
#'
#' @examples
#' \dontrun{
#' ## Using example data provided by MultiBaC package
#' data(example)
#'
#' omicsList <- list("RNA" = A.rna, "GRO" = A.gro)
#'
#' ## Create model
#' plsModel <- createPLSmodel(omicsList, crossval = NULL, test.comp = 3,
#' scale = FALSE, center = TRUE, regressor = 1)
#' }
createPLSmodel <- function(omicslist, test.comp, messages = TRUE,
                           scale = FALSE, center = TRUE,
                           crossval, regressor) {
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
    message(paste0("Input test.comp exceeds the minimum dimension of data matrix. test.comp set to ", max.comp))
    test.comp <- max.comp
  }

  # Create models ---------------------------------------------------------------
  models <- list()
  for ( i in seq_along(names(omicslist))[-regressor]) {

    # COmpute Q2 ------------------------------------------------------------------
    plsModel <- ropls::opls(t(omicslist[[regressor]]), t(omicslist[[i]]),
                            predI = test.comp, printL=FALSE, plotL = FALSE,
                            crossvalI = dim(t(omicslist[[1]]))[1], scaleC = pret)
    q2v <- plsModel@modelDF$`Q2(cum)`

    # Built final model -----------------------------------------------------------
    plsModel <- ropls::opls(t(omicslist[[1]]), t(omicslist[[i]]),
                            predI = which(q2v == max(q2v)), printL=FALSE, plotL = FALSE,
                            crossvalI = dim(t(omicslist[[1]]))[1], scaleC = pret)
    models[[(i-1)]] <- plsModel
  }
  names(models) <- names(omicslist)[-1]
  return(c(models, data.frame(q2v)))
}

#' getData
#'
#' @param ListOfBatches A list. Object returned by inputData function.
#'
#' @return A list of lists.
#'
#' @examples
getData <- function(ListOfBatches) {
  if ( inherits(ListOfBatches[[1]], "MultiAssayExperiment") ) {

    inputList <- list()

    for ( i in seq_along(names(ListOfBatches)) ) {
      inputList[[names(ListOfBatches)[i]]] <- ListOfBatches[[names(ListOfBatches)[i]]]@ExperimentList
    }


  } else {
    if ( is.null(cond.factor) ) {
      stop("cond.factor = NULL but a value is needed")
    }
    inputList <- ListOfBatches
  }
  return (inputList)
}

