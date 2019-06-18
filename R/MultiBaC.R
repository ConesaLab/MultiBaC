#' @import stats
#' @import graphics
#' @import abind
#' @import ggplot2
#' @import MultiAssayExperiment
#' @import ropls
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
                     commonOmic = 1,
                     scale = FALSE, center = TRUE,
                     showplot = TRUE, crossval = NULL,
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
    returnMultiAssay = FALSE
  }
  # browser()
  # Create PLS models ---------------------------------------------------------
  message("1: Create PLS models")
  modelList <- list()
  q2ofmodels <- list()
  for ( i in names(inputList)) {
    message(paste0("\t - Model for batch ",i))
    aux <- createPLSmodel(inputList[[i]], test.comp = test.comp, messages = FALSE,
                   it = it, tol = tol, scale = scale, center = center,
                   crossval, commonOmic)
    modelList[[i]] <- aux[1]
    q2ofmodels[[i]] <- aux[[2]]
  }

  # Generate missing omics ----------------------------------------------------
  message("2: Generating missing omics")
  missingOmics <- inputList
  for ( i in seq_along(inputList) ) {
    aux <- names(inputList)[-i]
    predictedOmics <- list()
    for ( n in aux ) {
      for ( j in seq_along(modelList[[n]])) {
        Ohat <- predict(modelList[[n]][[j]], newdata = t(inputList[[i]][[commonOmic]]))
        predictedOmics[[names(modelList[[n]])[j]]] <- t(Ohat)
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
    cond_factor <- c(cond_factor, cond.factor[[i]])
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

  # Applying ARSyN
  if (is.null(Fac)) {
    Fac <- c(ncol(cond_matrix)-1, ncol(batch_matrix)-1,
             (ncol(cond_matrix)*ncol(batch_matrix))-1,
             1)
  }

  # Correction step
  correctedOmics <- batchCorrection (multiBatchDesign, cond_matrix, batch_matrix, Fac)

  ## Common omic
  Xunfold <- NULL
  for ( i in seq_along(inputList) ) {
    Xunfold <- rbind(Xunfold, t(inputList[[i]][[commonOmic]]))
  }
  #
  # arsyn <- ASCA.2f(Xunfold, Designa = cond_matrix, Designb = batch_matrix,
  #                  Fac = Fac, type = 2, showvar = FALSE, showscree = FALSE)
  # arsynmodels[[1]] <- arsyn$Model.b$var.exp[,2]
  #
  # corrected <- Xunfold - arsyn$Model.b$TP - arsyn$Model.ab$TP
  # X.star <- list()
  # roff <- 0
  # for ( i in seq_along(X) ) {
  #   X.star[[i]] <- corrected[(roff+1):(roff + dim(X[[i]])[1]),]
  #   roff <- roff + dim(X[[i]])[1]
  # }
  # names(X.star) <- names(X)
  #
  # ## Non-common omics
  #
  # #### new approach
  # if (nm) {
  #   corrected_noncommon <- new_module(YZ, missingOmics, cond.factor)
  # } else {
  #  corrected_noncommon <- list()
  arsynmodels <- list()
    for ( i in (names(multiBatchDesign[[1]])) ) {
      omicB <- NULL
      for ( j in seq_along(multiBatchDesign)) {
        omicB <- rbind(omicB, t(multiBatchDesign[[j]][[i]]))
      }
      arsyn <- ASCA.2f(omicB, Designa = cond_matrix, Designb = batch_matrix,
                       Fac = Fac, type = 2, showvar = FALSE, showscree = FALSE)
      arsynmodels[[i]] <- c(0,arsyn$Model.b$var.exp[,2])

  #    corrected <- omics_to_correct[[i]] - arsyn$Model.b$TP - arsyn$Model.ab$TP
  #    corrected_noncommon[[i]] <- corrected[seq_len(dim(YZ[[i]])[1]),]
    }
  #  names(corrected_noncommon) <- names(YZ)
  # }

  # Q2 and explained batch-related variability plots --------------------------
  if(showplot) {
    MultiBaC_plots (q2ofmodels, arsynmodels, test.comp)
  }


  # Preparing results
  if (returnMultiAssay) {
    output <- ListOfBatches
    for (lab in names(ListOfBatches)) {
      for ( omic in names(ListOfBatches[[lab]]@ExperimentList)) {
        output[[lab]]@ExperimentList[[omic]] <- correctedOmics[[lab]][[omic]]
      }
    }
  } else {
    output <- correctedOmics
  }


  return(
    multibacClass(
      list("ListOfCorrectedBatches" = output,
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
#' @return Instance of class "opls"
#' @export
#'
#' @examples
#' \dontrun{
#' }
createPLSmodel <- function(omicslist, test.comp, it = 50, messages = TRUE,
                           tol = 1e-08, scale = FALSE, center = TRUE,
                           crossval, commonOmic) {
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

  # Create models ---------------------------------------------------------------
  models <- list()
  for ( i in seq_along(names(omicslist))[-commonOmic]) {

    # COmpute Q2 ------------------------------------------------------------------
    plsModel <- ropls::opls(t(omicslist[[commonOmic]]), t(omicslist[[i]]),
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


#' MultiBaC_plot
#'
#' @param modelList
#' @param arsynmodels
#'
#' @return
#' @export
#'
#' @examples
MultiBaC_plots <- function(q2values, arsynmodels, test.comp) {
  initpar <- par(c("mfrow"))
  on.exit(par(mfrow = initpar))
  par(mfrow = c(1,2), xpd = TRUE)

  # Q2 plot -------------------------------------------------------------------

  pallete <- colors()[c(11,17,51,56,29,512,97,653,136,24)]

  plot(1:3,1:3, type = "n", pch = 19,
       ylim = c(min(unlist(q2values)),1), xaxt = "n",
       xlim = c(1,test.comp+1), xlab="Number of Components", ylab = "Squared Q value",
       main = "Squared Q plot", bty = "n",
       cex.lab = 1.25, cex.axis = 1.25, font.lab = 2, cex.main=1.5)
  for ( i in seq_along(q2values)) {
    lines(1:test.comp, q2values[[i]], type = "b", pch = 19, col = pallete[i])
  }
  axis(1, seq_len(test.comp), c(seq_len(test.comp)), cex.axis = 1.25)
  legend(test.comp, 0.5,
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

  for ( i in seq_along(arsynmodels)) {
    lines(0:(length(arsynmodels[[1]])-1),
          c(arsynmodels[[i]]*100),
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


#' Title
#'
#' @param cond_matrix
#' @param batch_matrix
#' @param multiBatchDesign
#'
#' @return
#' @export
#'
#' @examples
batchCorrection <- function( multiBatchDesign, cond_matrix, batch_matrix, Fac) {

  # Max predI
  ncomp <- min(unlist(lapply(multiBatchDesign,function(x) {
    min(unlist(lapply(x, dim)))
  })))

  # Create block to apply ARSyN
  boff <- c()
  ooff <- c()
  lsmodels <- list()
  for ( i in seq_along(multiBatchDesign)) {
    lsmodels [[names(multiBatchDesign)[i]]] <- list()
    for ( j in seq_along(multiBatchDesign[[i]])) {
      pc <- opls(t(multiBatchDesign[[i]][[j]]), predI = ncomp, scaleC = "center",
                 printL = FALSE, plotL = FALSE, crossvalI = dim(t(multiBatchDesign[[i]][[j]]))[1])
      lsmodels[[names(multiBatchDesign)[i]]][[names(multiBatchDesign[[i]])[j]]] <- pc
    }
  }

  #####

  blocksomics <- NULL
  cond.matrix <- NULL
  batch.matrix <- NULL

  for ( i in (names(lsmodels[[1]]))) {
    cond.matrix <- rbind(cond.matrix, cond_matrix)
    batch.matrix <- rbind(batch.matrix, batch_matrix)
    roff <- 0
    for ( j in seq_along(lsmodels)) {
      t <- table(sign(lsmodels[[j]][[i]]@scoreMN[which(cond_matrix[(roff+1):(roff+dim(multiBatchDesign[[j]][[i]])[2]),1]==1),1]))
      s <- which(t==max(t))
      lsmodels[[j]][[i]]@scoreMN <- lsmodels[[j]][[i]]@scoreMN * s
      lsmodels[[j]][[i]]@loadingMN <- lsmodels[[j]][[i]]@loadingMN * s
      blocksomics <- rbind(blocksomics, lsmodels[[j]][[i]]@scoreMN)
      roff <- roff +  dim(multiBatchDesign[[j]][[i]])[2]
    }
  }
  # browser()
  asca <- ASCA.2f(blocksomics, Designa = cond.matrix, Designb = batch.matrix,
                             Fac = Fac, type = 2, showvar = FALSE, showscree = FALSE)
  blocks.corr <- blocksomics - asca$Model.b$TP

  # Back to original structure
  correctedOmics <- list()
  roff <- 0
  for ( i in (names(lsmodels[[1]]))) {
    for ( j in seq_along(lsmodels)) {
      s <- blocks.corr[(roff+1):(roff+dim(lsmodels[[j]][[i]]@scoreMN)[1]),]
      l <- lsmodels[[j]][[i]]@loadingMN
      r <- apply(s%*%t(l),1,function(x) {x+lsmodels[[j]][[i]]@xMeanVn})
      correctedOmics[[names(lsmodels)[j]]][[i]] <- r
      roff <- roff + dim(lsmodels[[j]][[i]]@scoreMN)[1]
    }
  }

  return (correctedOmics)

}



