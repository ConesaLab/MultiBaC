#' @import graphics
#' @import stats
#' @import ggplot2
#' @include auxfunctions.R
#' @include ASCA2f.R
NULL

#' MultiBaC
#'
#' MultiBaC performs a multi-omic, multi-batch correction
#'
#'
#' @docType package
#' @name MultiBaC
NULL
#> NULL


#' MultiBaC
#'
#' MultiBaC is a strategy to correct batch effects from multiomic datasets distributed across different labs or data acquisition events.
#' MultiBaC is the first Batch effect correction algorithm that dealing with batch effect correction in multiomics datasets.
#' MultiBaC is able to remove batch effects across different omics generated within separate batches provided that at least one common
#' omic data type is included in all the batches considered.
#'
#' @param ListOfBatches A list. Object returned by inputData function.
#' @param test.comp maximum number of components allowed in PLS models.
#' @param scale TRUE or FALSE. Whether X and Y matrices have to be scaled
#' @param center TRUE or FALSE. Whether X and Y matrices have to be centered
#' @param cond.factor A list with one slot for each batch containing the common experimental setting.
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
#' If NULL (default) the function optimize each value.
#' @param commonOmic Name or index of the commonOmic between the batches. If the index is given, the common data matrix must be in the same possition in each batch.
#' @param showplot If TRUE, the Q^2 and the Fac optimization will be plotted.
#' @param crossval Integer: number of cross-validation segments. The number of samples (rows of 'x') must be at least >= crossvalI.
#' If NULL (default) leave-one-out crossvalidation is performed
#'
#' @return A list with as many slots as the given number of batches.
#' Returns the same structure than the input data "ListOfBatches" but with batch-corrected data.
#' @export
#'
#' @examples
#' \dontrun{
#' ## Using example data provided by MultiBaC package
#' data(multiyeast)
#' inputData <- inputData(A.rna, A.gro, B.rna, B.ribo, C.rna, C.par, batches = c(1,1,2,2,3,3),
#' omicNames = c("RNA", "GRO", "RNA", "RIBO", "RNA", "PAR"),
#' batchesNames = c("A", "B", "C"))
#'
#' ## Creating cond.factor
#' cond.factor = list("A" = c("Glu+", "Glu+", "Glu+", "Glu-", "Glu-", "Glu+"),
#' "B" = c("Glu+", "Glu+", "Glu-", "Glu-"),
#' "C" = c("Glu+", "Glu+", "Glu-", "Glu-"))
#'
#' ## Applying MultiBaC
#' Corrected_Matrices <- MultiBaC (inputData, test.comp = 5,
#' cond.factor = cond.factor, commonOmic = "RNA")}
#'
MultiBaC <- function(ListOfBatches, test.comp, cond.factor,
                     commonOmic = 1,
                     scale = FALSE, center = TRUE,
                     showplot = TRUE, crossval = NULL,
                     Fac = NULL) {
  # Input data structure ------------------------------------------------------
  if (is.null(names(ListOfBatches))) {
    stop ("Stop at input evaluation: Elements in ListOfBatches must be named")
  }
  #inputList <- getData(ListOfBatches)

  # Create PLS models ---------------------------------------------------------
  message("1: Create PLS models")

  aux.1 <- genModelList(ListOfBatches, test.comp = test.comp,
                            scale = scale, center = center,
                            crossval = crossval, commonOmic = commonOmic)
  modelList <- aux.1$modelList
  q2ofmodels <- aux.1$q2ofmodels

  # Generate missing omics ----------------------------------------------------
  message("2: Generating missing omics")
  multiBatchDesign <- genMissingOmics(ListOfBatches, modelList, commonOmic)

  # Batch effect correction using ARSyN ---------------------------------------
  message("3: Batch effect correction using ARSyN")

  # Correction step
  aux.2 <- batchCorrection (multiBatchDesign, cond.factor, Fac)
  correctedOmics <- aux.2$correctedOmics
  ascamodels <- aux.2$ascamodels

  # Q2 and explained batch-related variability plots ---------------------------
  if(showplot) {
    MultiBaC_plots (q2ofmodels, ascamodels)
  }


  # Preparing results ----------------------------------------------------------
  output <- ListOfBatches
  for (lab in names(ListOfBatches)) {
    for ( omic in names(ListOfBatches[[lab]]@ExperimentList)) {
      output[[lab]]@ExperimentList[[omic]] <- correctedOmics[[lab]][[omic]]
    }
  }


  return(c(lapply(output, function(x) x)))
}


#' inputData
#'
#' This function creates a list object to be used by MultiBaC function from a set of matrix R objects.
#'
#' @param ... Matrices (features x samples) sepated by comma.
#' @param batches A vector or type factor. Indicates which batch each input matrix belongs to.
#' @param omicNames The names of input omics.
#' @param batchesNames The names of the different batches,
#'
#' @return A list containing one slot for each batch. Each slot contains a MultiAssayExperiment class object with the corresponding omics included.
#' @export
#'
#' @examples
#' \dontrun{
#' #' ## Using example data provided by MultiBaC package
#' data(multiyeast)
#' inputData <- inputData(A.rna, A.gro, B.rna, B.ribo, C.rna, C.par, batches = c(1,1,2,2,3,3),
#' omicNames = c("RNA", "GRO", "RNA", "RIBO", "RNA", "PAR"),
#' batchesNames = c("A", "B", "C"))
#' }
inputData <- function(..., batches, omicNames, batchesNames) {
  inputOmics <- list(...)
  omicList <- sapply(unique(batches), function(x) {
    aux.list <- inputOmics[which(batches==x)]
    names(aux.list) <- omicNames[which(batches==x)]
    MultiAssayExperiment::MultiAssayExperiment(aux.list)
  })

  names(omicList) <- batchesNames
  return(omicList)
}

#' batchEstPlot
#'
#' This function uses linear models to estimate the batch effect magnitude using the common data across batches. It compares
#' the result with theoretical distribution of diferrent levels of batch magnitude.
#'
#' @param ListOfBatches A list. Object returned by inputData function.
#' @param commonOmic Name or index of the commonOmic between the batches. If the index is given, the common data matrix must be in the same possition in each batch.
#'
#' @return An object of class ggplot.
#' @export
#'
#' @examples
#' \dontrun{
#' #' ## Using example data provided by MultiBaC package
#' data(multiyeast)
#' inputData <- inputData(A.rna, A.gro, B.rna, B.ribo, C.rna, C.par, batches = c(1,1,2,2,3,3),
#' omicNames = c("RNA", "GRO", "RNA", "RIBO", "RNA", "PAR"),
#' batchesNames = c("A", "B", "C"))
#'
#' batchEstPlot(inputData, commonOmic = 1)
#' }
batchEstPlot <- function(ListOfBatches, commonOmic) {

  initpar <- par(c("mai", "pin", "xpd"))
  mai <- c(initpar$mai[1:3], initpar$pin[2]-0.15)
  on.exit(par(mai = initpar$mai, xpd = initpar$xpd))
  par(mai = mai)

  # Input Data Structure -------------------------------------------------------
  inputList <- getData(ListOfBatches)

  # Extract common omic --------------------------------------------------------
  Xunfold <- NULL
  for ( lab in seq_along(inputList)) {
    Xunfold <- rbind(Xunfold, t(inputList[[lab]][[commonOmic]]))
  }

  # Create batch factor --------------------------------------------------------
  batch_factor <- c()
  for ( i in seq_along(inputList) ) {
    batch_factor <- c(batch_factor, rep(i, dim(inputList[[i]][[1]])[2]))
  }

  # Calculate Overall Mean -----------------------------------------------------

  offset<-apply(Xunfold,2,mean)
  Xoff<-Xunfold-(cbind(matrix(1,nrow=nrow(Xunfold),ncol=1))%*%rbind(offset))

  # Compute coefficients -------------------------------------------------------

  design.matrix <- model.matrix(~ factor(batch_factor))
  #
  design.matrix[which(design.matrix==0)] <- -1
  #
  beta.hat <- solve(t(design.matrix)%*%design.matrix) %*% t(design.matrix) %*% Xoff



  # Theoretical distributions --------------------------------------------------

  beta.hat_v <- as.vector(t(beta.hat))

  reg_v <- c()
  for ( i in 1:(length(beta.hat_v)/dim(beta.hat)[2])) {
    reg_v <- c(reg_v, rep(i,dim(beta.hat)[2]))
  }

  toplot <- data.frame("Coef" = (reg_v[(dim(beta.hat)[2]+1):length(beta.hat_v)]),
                       "Values" = beta.hat_v[(dim(beta.hat)[2]+1):length(beta.hat_v)])
  for ( i in seq_along(ListOfBatches)[-1]) {
    toplot$Coef[which(toplot$Coef==i)] <- names(ListOfBatches)[i]
  }



  # The label
  label = paste0("* Batch ", names(ListOfBatches)[1], " as\nreference")

  # The  plot
  p <- ggplot(toplot, aes(x=Coef, y=Values)) +
    geom_violin() +
    labs(title="Batch Magnitude Plot", x = "Batch", y = "Coefficients per Gene", size=2,
         col = "Theoretical Batch \nMagnitudes",
         caption = label, font_face = "italic") +
    annotate("text", x = Inf, y = -3, label = label, hjust = -0.2, size = 5) +
    scale_fill_brewer(palette="Blues") +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
    theme_classic(base_family = "") +
    theme(axis.text=element_text(size=18),
          title=element_text(size=16,face="bold"),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size=20),
          legend.title = element_text(size = 20, face = "plain"),
          legend.spacing = unit(.25,"cm")) +
    geom_hline(aes(yintercept=6, col = "High"), lty = 5, show.legend = TRUE) +
    geom_hline(aes(yintercept=-6, col = "High"), lty = 5,show.legend = TRUE) +
    geom_hline(aes(yintercept=4, col = "Medium"), lty = 5,show.legend = TRUE) +
    geom_hline(aes(yintercept=-4, col = "Medium"), lty = 5,show.legend = TRUE) +
    geom_hline(aes(yintercept=2, col = "Low"), lty = 5,show.legend = TRUE) +
    geom_hline(aes(yintercept=-2, col = "Low"), lty = 5,show.legend = TRUE)

  p
}

#' batchCorrection
#'
#' This function performs the ARSyN correction [1] for each omic contained in mulBatchDesign input object.
#'
#' @param multiBatchDesign A list containing the original and predicted omic for each batch. All omics must be present in every batch. Output object of genMissingOmics function
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
#' If NULL (default) the function optimize each value.
#' @param cond.factor A list with one slot for each batch containing the common experimental setting.
#'
#' @return A list with tow slots. The first one, "correctedOmics", has the same structure than multiBatchDesign input object, but containing corrected information.
#' The second slot contains a list with the batch effect estimation of each omic.
#' @export
#' @references [1] Nueda MJ, Ferrer A, Conesa A. ARSyN: A method for the identification and removal of systematic noise in multifactorial time course microarray experiments. Biostatistics. 2012;13(3):553-566. doi:10.1093/biostatistics/kxr042
#' @examples
#' \dontrun{
#' ## Using example data provided by MultiBaC package
#' data(multiyeast)
#' inputData <- inputData(A.rna, A.gro, B.rna, B.ribo, C.rna, C.par, batches = c(1,1,2,2,3,3),
#' omicNames = c("RNA", "GRO", "RNA", "RIBO", "RNA", "PAR"),
#' batchesNames = c("A", "B", "C"))
#'
#' ## Creating cond.factor
#' cond.factor = list("A" = c("Glu+", "Glu+", "Glu+", "Glu-", "Glu-", "Glu+"),
#' "B" = c("Glu+", "Glu+", "Glu-", "Glu-"),
#' "C" = c("Glu+", "Glu+", "Glu-", "Glu-"))
#'
#' ## Create PLS models
#' modelList <- genModelList(input, test.comp = 5)$modelList
#'
#' ## Generate missing omic data
#' missingOmics <- genMissingOmics(input, modelList, commonOmic = 1)
#'
#' ## Batch correction across omics
#'
#' correctedOmics <- batchCorrection(missingOmics, cond.factor)$correctedOmics}
batchCorrection <- function(multiBatchDesign, cond.factor, Fac = NULL) {

  # Making model matrices ------------------------------------------------------
  cond_factor <- c()
  for ( i in seq_along(cond.factor) ) {
    cond_factor <- c(cond_factor, cond.factor[[i]])
  }
  cond_matrix <- model.matrix(~ 0 + factor(cond_factor))

  batch_factor <- c()
  for ( i in seq_along(multiBatchDesign) ) {
    batch_factor <- c(batch_factor, rep(i, dim(multiBatchDesign[[i]][[1]])[2]))
  }
  batch_matrix <- model.matrix(~ 0 + factor(batch_factor))

  # Defining Fac ---------------------------------------------------------------
  if (is.null(Fac)) {
    Fac <- c(ncol(cond_matrix)-1, ncol(batch_matrix)-1,
             ((ncol(cond_matrix)*ncol(batch_matrix)))-1,
             1)
  }
  # Omic correction ------------------------------------------------------------
  correctedOmics <- multiBatchDesign
  ascamodels <- list()
  for ( omic in names(multiBatchDesign[[1]])) {
    Xunfold <- NULL
    set_lims <- 0
    for ( lab in names(multiBatchDesign)) {
      Xunfold <- rbind(Xunfold, t(multiBatchDesign[[lab]][[omic]]))
      set_lims <- c(set_lims, (set_lims[length(set_lims)])+dim(t(multiBatchDesign[[lab]][[omic]]))[1])
    }
    asca <- MultiBaC:::ASCA.2f(Xunfold, Designa = cond_matrix, Designb = batch_matrix,
                              Fac = Fac, type = 2, showvar = FALSE, showscree = FALSE)
    Xunfold.r <- Xunfold - asca$Model.b$TP - asca$Model.ab$TP #- asca$Model.res$TP
    ascamodels[[omic]] <- c(0,asca$Model.b$var.exp[,2])

    for( i in seq_along(multiBatchDesign)) {
      correctedOmics[[names(multiBatchDesign)[i]]][[omic]] <- t(Xunfold.r[(set_lims[i]+1):(set_lims[i+1]),])
    }
  }

  return (list("correctedOmics" = correctedOmics,
               "ascamodels" = ascamodels))
}

#' genModelList
#'
#' This function performs PLS models for every batch. A PLS model is generated for each non-common omic in each batch.
#'
#' @param ListOfBatches A list. Object returned by inputData function.
#' @param test.comp Maximum number of components allowed in PLS models.
#' @param scale TRUE or FALSE. Whether X and Y matrices have to be scaled
#' @param center TRUE or FALSE. Whether X and Y matrices have to be centered
#' @param crossval Integer: number of cross-validation segments. The number of samples (rows of 'x') must be at least >= crossvalI.
#' If NULL (default) leave-one-out crossvalidation is performed
#' @param commonOmic Name or index of the commonOmic between the batches. If the index is given, the common data matrix must be in the same possition in each batch.
#'
#' @return A list with two slots:
#' \enumerate{
#'  \item{modelList: A list of PLS models.}
#'  \item{q2ofmodels: A list. Each slot contains a vector of Q^2 values of a PLS model.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' ## Using example data provided by MultiBaC package
#' data(multiyeast)
#' inputData <- inputData(A.rna, A.gro, B.rna, B.ribo, C.rna, C.par, batches = c(1,1,2,2,3,3),
#' omicNames = c("RNA", "GRO", "RNA", "RIBO", "RNA", "PAR"),
#' batchesNames = c("A", "B", "C"))
#'
#' ## Creating cond.factor
#' cond.factor = list("A" = c("Glu+", "Glu+", "Glu+", "Glu-", "Glu-", "Glu+"),
#' "B" = c("Glu+", "Glu+", "Glu-", "Glu-"),
#' "C" = c("Glu+", "Glu+", "Glu-", "Glu-"))
#'
#' ## Create PLS models
#' modelList <- genModelList(input, test.comp = 5)$modelList
#'
#' ## Generate missing omic data
#' missingOmics <- genMissingOmics(input, modelList, commonOmic = 1)
#'
#' ## Batch correction across omics
#'
#' correctedOmics <- batchCorrection(missingOmics, cond.factor)$correctedOmics
#' }
genModelList <- function(ListOfBatches, test.comp, scale = FALSE, center = TRUE,
                         crossval = NULL, commonOmic = 1) {
  # Get matrices ---------------------------------------------------------------
  inputList <- getData(ListOfBatches)

  # Create models --------------------------------------------------------------
  modelList <- list()
  q2ofmodels <- list()
  for ( i in names(inputList)) {
    message(paste0("\t - Model for batch ",i))
    aux <- createPLSmodel(inputList[[i]], test.comp = test.comp, messages = FALSE,
                          scale = scale, center = center,
                          crossval, commonOmic)
    modelList[[i]] <- aux[1]
    q2ofmodels[[i]] <- aux[[2]]
  }

  return (list("modelList" = modelList, "q2ofmodels" = q2ofmodels))
}

#' genMissingOmics
#'
#' This function generates for all the batches the omic data they had not originally. This is the previous step to apply ARSyN [1] correction.
#'
#' @param ListOfBatches A list. Object returned by inputData function
#' @param modelList A list of PLS models. Oject returned by genModelList function
#' @param commonOmic Name or index of the commonOmic between the batches. If the index is given, the common data matrix must be in the same possition in each batch.
#'
#' @return A type list object.
#' @export
#' @references [1] Nueda MJ, Ferrer A, Conesa A. ARSyN: A method for the identification and removal of systematic noise in multifactorial time course microarray experiments. Biostatistics. 2012;13(3):553-566. doi:10.1093/biostatistics/kxr042
#' @examples
#' \dontrun{
#' ## Using example data provided by MultiBaC package
#' data(multiyeast)
#' inputData <- inputData(A.rna, A.gro, B.rna, B.ribo, C.rna, C.par, batches = c(1,1,2,2,3,3),
#' omicNames = c("RNA", "GRO", "RNA", "RIBO", "RNA", "PAR"),
#' batchesNames = c("A", "B", "C"))
#'
#' ## Creating cond.factor
#' cond.factor = list("A" = c("Glu+", "Glu+", "Glu+", "Glu-", "Glu-", "Glu+"),
#' "B" = c("Glu+", "Glu+", "Glu-", "Glu-"),
#' "C" = c("Glu+", "Glu+", "Glu-", "Glu-"))
#'
#' ## Create PLS models
#' modelList <- genModelList(input, test.comp = 5)$modelList
#'
#' ## Generate missing omic data
#' missingOmics <- genMissingOmics(input, modelList, commonOmic = 1)
#'
#' ## Batch correction across omics
#'
#' correctedOmics <- batchCorrection(missingOmics, cond.factor)$correctedOmics
#' }
genMissingOmics <- function(ListOfBatches, modelList, commonOmic) {

  # Get matrices ---------------------------------------------------------------
  inputList <- getData(ListOfBatches)

  # Generate missing omics -----------------------------------------------------
  missingOmics <- inputList
  for ( i in seq_along(inputList) ) {
    aux <- names(inputList)[-i]
    predictedOmics <- list()
    for ( n in aux ) {
      for ( j in seq_along(modelList[[n]])) {
        Ohat <- ropls::predict(modelList[[n]][[j]], newdata = t(inputList[[i]][[commonOmic]]))
        predictedOmics[[names(modelList[[n]])[j]]] <- t(Ohat)
      }
    }
    missingOmics [[names(inputList)[i]]] <-  predictedOmics
  }

  # Creating omic blocks to correct --------------------------------------------
  multiBatchDesign <- inputList
  for ( i in seq_along(inputList) ) {
    multiBatchDesign[[i]] <- c(inputList[[i]], missingOmics[[i]])
  }
  return(multiBatchDesign)
}

