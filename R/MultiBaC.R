#' @import stats
#' @import graphics
#' @import abind
#' @import ggplot2
#' @include plsClass.R
#' @include auxfunctions.R
#' @include ASCA2f.R
NULL

#' MultiBaC
#'
#' Description
#'
#'
#' @docType package
#' @name MultiBaC
NULL
#> NULL


#' MultiBaC
#'
#' @param X list of common omics
#' @param YZ list of non-common omics
#' @param test.comp maximum number of components allowed
#' @param it maximum number of iterations allowed for nipals algorithm
#' @param tol threshold value to validate nipals algorithm convergence
#' @param scale TRUE or FALSE. Whether X and Y matrices have to be scaled
#' @param center TRUE or FALSE. Whether X and Y matrices have to be centered
#' @param cond.factor list of experimental design conditions for each lab
#' @param estBatchmagnitude TRUE or FALSE. Whether to show an estimation of batch effect magnitude or not
#' @param Fac Vector with numbers of components to extract in each submodel, by order:
#' \enumerate{
#'     \item Model a (time)
#'     \item Model b (second factor)
#'     \item Model c (third factor)
#'     \item Model ab (interaction) to not consider interations set to 0
#'     \item Model ac (interaction) to not consider interations set to 0
#'     \item Model bc (interaction) to not consider interations set to 0
#'     \item Model abc (interaction) to not consider interations set to 0
#'     \item Number components of residues
#' }
#' @param type Vector indicating whether the analyses of the model of the factors that interact with time (b.ab and c.ac) are studied jointly(1) or separately(2).
#' @param showvar Show variability associated to each model or not.
#' @param showscree Show a screenplot.
#' @param nm TRUE or FALSE. Whether to compute a final step to improve batch correction.
#'
#' @return A list with two slots
#' #' \describe{
#'     \item{X.star}{List of common omics matrices corrected}
#'     \item{YZ.star}{List of non-common omics matrices corrected}
#' }
#' @export
#'
#' @examples
#' Corrected_Matrices <- MultiBaC ( X = list("A" = rnaseq_a, "B" = rnaseq_b),
#'                                  YZ = list("A" = chipseq_a, "B" = metabolomics_b),
#'                                  test.comp = 3, cond.factor = design)
MultiBaC <- function(ListOfBatches, test.comp, cond.factor,
                     scale = FALSE, center = TRUE,
                     showplot = TRUE,
                     Fac = NULL,
                     nm = FALSE) {
  # Input data structure ------------------------------------------------------
  # do.call(inputeval, args = list(ListOfBatches = ListOfBatches,
  #                                test.comp = test.comp, cond.factor = cond.factor,
  #                                scale = scale, center = center,
  #                                showplot = showplot,
  #                                Fac = Fac))

  if ( inherits(ListOfBatches[[1]], "MultiAssayExperiment") ) {

    inputList <- list()

    for ( i in seq_along(names(ListOfBatches)) ) {
      inputList[[names(ListOfBatches)[i]]] <- ListOfBatches[[names(ListOfBatches)[i]]]@ExperimentList
    }

    returnMultiAssay = TRUE

  } else {
    if ( is.null(cond.factor) ) {
      stop("cond.factor = NULL but a value is needed")
    }
    inputList <- ListOfBatches
  }
  # browser()
  # Create PLS models ---------------------------------------------------------
  message("1: Create PLS models")
  modelList <- lapply(seq_along(inputList), function (index) {
    message(paste0("\t - Model for batch ",index))
    createPLSmodel(inputList[[index]], test.comp = test.comp, messages = FALSE,
                   it = it, tol = tol, scale = scale, center = center, plot = Q2plot)
  })
  names(modelList) <- names(inputList)

  # Generate missing omics ----------------------------------------------------
  message("2: Generating missing omics")
  missingOmics <- inputList
  for ( i in seq_along(inputList) ) {
    aux <- names(inputList)[-i]
    predictedOmics <- list()
    for ( n in aux ) {
      for ( j in seq_along(modelList[[n]])) {
        Ohat <- predict(modelList[[n]][[j]], newdata = t(inputList[[index]][[1]]))
        predictedOmics[[names(modelList[[n]])[j]]] <- Ohat
      }
    }
    missingOmics [[names(inputList)[i]]] <-  predictedOmics
  }

  # Batch effect correction using ARSyN ---------------------------------------
  message("3: Batch effect correction using ARSyN")
  omics_to_correct <- list()

  # Making model matrices
  cond_factor <- c()
  for ( i in seq_along(cond.factor) ) {
    cond_factor <- c(cond_factor, data.frame(cond.factor)[[i]])
  }
  cond_matrix <- model.matrix(~ 0 + factor(cond_factor))

  batch_factor <- c()
  for ( i in seq_along(inputList) ) {
    batch_factor <- c(batch_factor, rep(i, dim(inputList[[i]][[1]])[2]))
  }
  batch_matrix <- model.matrix(~ 0 + factor(batch_factor))

  # Creating omic blocks to correct
  multiBatchDesign <- inputList
  for ( i in seq_along(inputList) ) {
    multiBatchDesign[[i]] <- c(inputList[[i]], missingOmics[[i]])
  }

  #### insert function to correct


  # Applying ARSyN
  arsynmodels <- list()
  if (is.null(Fac)) {
    Fac <- c(ncol(cond_matrix)-1, ncol(batch_matrix)-1,
            (ncol(cond_matrix)*ncol(batch_matrix))-1,
            1)
  }

  ## Common omic
  Xunfold <- NULL
  for ( i in seq_along(X) ) {
    Xunfold <- rbind(Xunfold, X[[i]])
  }

  arsyn <- ASCA.2f(Xunfold, Designa = cond_matrix, Designb = batch_matrix,
                   Fac = Fac, type = 2, showvar = FALSE, showscree = FALSE)
  arsynmodels[[1]] <- arsyn$Model.b$var.exp[,2]

  corrected <- Xunfold - arsyn$Model.b$TP - arsyn$Model.ab$TP
  X.star <- list()
  roff <- 0
  for ( i in seq_along(X) ) {
    X.star[[i]] <- corrected[(roff+1):(roff + dim(X[[i]])[1]),]
    roff <- roff + dim(X[[i]])[1]
  }
  names(X.star) <- names(X)

  ## Non-common omics

  #### new approach
  if (nm) {
    corrected_noncommon <- new_module(YZ, missingOmics, cond.factor)
  } else {
    corrected_noncommon <- list()
    for ( i in seq_along(omics_to_correct) ) {
      arsyn <- ASCA.2f(omics_to_correct[[i]], Designa = cond_matrix, Designb = batch_matrix,
                       Fac = Fac, type = 2, showvar = FALSE, showscree = FALSE)
      arsynmodels[[i+1]] <- arsyn$Model.b$var.exp[,2]

      corrected <- omics_to_correct[[i]] - arsyn$Model.b$TP - arsyn$Model.ab$TP
      corrected_noncommon[[i]] <- corrected[seq_len(dim(YZ[[i]])[1]),]
    }
    names(corrected_noncommon) <- names(YZ)
  }
  names(arsynmodels) <- c("common", names(omics_to_correct))

  # Q2 and explained batch-related variability plots --------------------------
  if(showplot) {
    MultiBaC_plots (modelList, arsynmodels, test.comp)
  }


  return(
    multibacClass(
      list("X.star" = X.star,
           "YZ.star" = corrected_noncommon,
           "batchMagnitudePlot" = batchEstimation(Xunfold, batch_factor))
      )
    )
}


#' createPLSmodel
#'
#' @param X matrix of predictor variables
#' @param Y matrix of response variables
#' @param test.comp maximum number of latent variables (or components) to test
#' @param it maximum number of iterations allowed for nipals algorithm
#' @param tol threshold value to validate nipals algorithm convergence
#' @param scale TRUE or FALSE. Whether X and Y matrices have to be scaled
#' @param center TRUE or FALSE. Whether X and Y matrices have to be centered
#' @param plot TRUE or FALSE. Whether to plot Q^2 across components.
#'
#' @return Instance of class "plsClass"
#' @export
#'
#' @examples
#' \dontrun{
#' }
createPLSmodel <- function(omicslist, test.comp, it = 50, messages = TRUE,
                           tol = 1e-08, scale = FALSE, center = TRUE, plot = TRUE) {
  # Set preprocessing -----------------------------------------------------------
  if ( center == FALSE & scale == FALSE) {
    pret <- "none"
  } else if (center == TRUE & scale == FALSE) {
    pret <- "center"
  } else {
    pret <- "standard"
  }

  # Create models ---------------------------------------------------------------
  models <- list()
  for ( i in seq_along(names(omicslist))[-1]) {

    # COmpute Q2 ------------------------------------------------------------------
    plsModel <- ropls::opls(t(omicslist[[1]]), t(omicslist[[i]]),
                            predI = test.comp, printL=FALSE, plotL = FALSE,
                            crossvalI = , scaleC = pret)


    # Optimize number of components -----------------------------------------------
    new.comp <- optimizeComponents(plsModel)

    # Built final model -----------------------------------------------------------
    plsModel <- ropls::opls(t(omicslist[[1]]), t(omicslist[[i]]),
                            predI = new.comp, printL=FALSE, plotL = FALSE,
                            crossvalI = , scaleC = pret)
    models[[(i-1)]] <- plsModel
  }
  names(models) <- names(omicslist)[-1]
  models
}


#' MultiBaC_plot
#'
#' @param modelList
#' @param arsynmodels
#'
#' @return
#' @export
#'
#' @examples
MultiBaC_plots <- function(modelList, arsynmodels, test.comp) {
  initpar <- par(c("mfrow"))
  on.exit(par(mfrow = initpar))
  par(mfrow = c(1,2), xpd = TRUE)

  # Q2 plot -------------------------------------------------------------------

  q2values <- vapply(modelList, function(model) {
    names(model)
    model$Q2
  }, rep(0,length(modelList[[1]]$Q2)))

  pallete <- colors()[c(11,17,51,56,29,512,97,653,136,24)]

  plot(q2values[,1], type = "n", pch = 19,
       ylim = c(min(q2values), max(q2values)), xaxt = "n",
       xlim = c(1,test.comp+1), xlab="Number of Components", ylab = "Squared Q value",
       main = "Squared Q plot", bty = "n",
       cex.lab = 1.25, cex.axis = 1.25, font.lab = 2, cex.main=1.5)
  for ( i in seq_along(modelList)) {
    lines(q2values[,i], type = "b", pch = 19, col = pallete[i])
  }
  axis(1, seq_len(test.comp), c(seq_len(test.comp)), cex.axis = 1.25)
  legend(test.comp, mean(c(min(q2values), max(q2values))),
         bty = "n", title = expression(underline(Batches)),
         legend = c(names(modelList)),
         col = c(pallete),
         cex = 1.5, lty = c(1,1), pch = c(19,19))

  # Explained batch-related variability plot ----------------------------------
  plot(c(0,0), type = "n", pch = 19,
       ylim = c(0,100), xaxt = "n",
       xlim = c(0,3), xlab="Number of Components",
       ylab = "Explained batch-related variability (%)",
       main = "ARSyN nº of components", bty = "n",
       cex.lab = 1.25, cex.axis = 1.25, font.lab = 2, cex.main=1.5)
  pallete <- colors()[c(11,17,51,56,29,512,97,653,136,24)]

  for ( i in seq_along(arsynmodels)) {
    lines(c(0,seq_along(arsynmodels[[i]])),
          c(0, arsynmodels[[i]]*100),
          type = "b", pch = 19, col = pallete[i])
  }

  axis(1, 0:2, 0:2, cex.axis = 1.25)
  legend(0.8, 50,
         bty = "n", title = expression(underline(Omics)),
         legend = c("common","non-common: ", names(arsynmodels)[-1]),
         col = c(pallete[1], "white", pallete[2:4]),
         cex = 1.5, lty = c(1,0,rep(1,length(arsynmodels)-1)), pch = rep(19,length(arsynmodels)+1))

  # Advertising about superposition
  if (sum(arsynmodels[[1]] == arsynmodels[[2]])) {
    message("Caution: Explained variance could be similar for more than two omics. Thus lines and dots could be superpossed")
  }
}






