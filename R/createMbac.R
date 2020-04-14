#' createMbac
#'
#' This function creates a list object to be used by MultiBaC function from a set of matrix R objects.
#'
#' @param inputOmics A list containing all the matrices or data.frame objects to be analysed. MultiAssayExperiment objects can alternatively be provided.
#' @param batchFactor Either a vector or a factor indicating the batch were each input matrix belongs to (i.e. study, lab, time point, etc.). If NULL (default) no batch is considered and just ARSyNbac noise reduction mode could be applied.
#' @param experimentalDesign A list with as many elements as batches. Each element can be a factor, a character vector or a data.frame indicating the experimental conditions for each sample in that batch. When being a data.frame with more than one column (multi-factorial experimental designs), the different columns will be combined into a single one to be used by MultiBaC or ARSyNbac. In any case, the experimental setting must be the same for all batches. In addition, the names of the elements in this list must be the same as declared in batches argument. If not (or if NULL), names are forced to be the same in as in batches argument and in the same order.
#' @param omicNames Vector of names for each input matrix. The common omic is required to have the same name across batches.
#' @param commonOmic Name of the common omic between the batches. It must be one of the names in omicNames argument. If NULL (default), the omic name which is common to all batches is selected as commonOmic.
#'
#' @return Custom mbac object. Elements in a mbac object:
#' \enumerate{
#'     \item ListOfBatches: A list of MultiAssayExperiment objects (one per batch).
#'     \item commonOmic Name of the common omic between batches.
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' my_mbac <- createMbac (inputOmics = list(A.rna, A.gro, B.rna, B.ribo, C.rna, C.par),
#'                        batchFactor = c("A", "A", "B", "B", "C", "C"),
#'                        experimentalDesign = list("A" =  c("Glu+", "Glu+",
#'                        "Glu+", "Glu-", "Glu-", "Glu-"),
#'                        "B" = c("Glu+", "Glu+", "Glu-", "Glu-"),
#'                        "C" = c("Glu+", "Glu+", "Glu-", "Glu-")),
#'                        omicNames = c("RNA", "GRO", "RNA", "RIBO", "RNA", "PAR"),
#'                        commonOmic = "RNA")
#' }
createMbac <- function(inputOmics, batchFactor = NULL,
                       experimentalDesign, omicNames,
                       commonOmic = NULL) {

  # Checking inputs ------------------------------------------------------------
  if ( !is.list(inputOmics)) {
    inputOmics <- list(inputOmics)
  }
  if (is.null(batchFactor)) {
    batchFactor <- c("Batch1")
  }
  batches <- as.character(unique(batchFactor))
  if (class(experimentalDesign) == "list") {
    namexp <- names(experimentalDesign)
    if (!is.null(namexp)) {
      if (!(length(intersect(batches, namexp)) == length(batches) && length(batches) == length(namexp))) {
        stop("Error: Batches in batchFactor and names in experimentalDesign do not match." )
      }
    } else {
      if (length(experimentalDesign) > 1) {
        stop("Error: experimentlDesign must be named." )
      } else {
        experimentalDesign <- list("Batch1" = experimentalDesign[[1]])
      }
    }
  } else {
    experimentalDesign <- list("Batch1" = experimentalDesign)
  }

  if ((length(levels(factor(batchFactor))) == length(inputOmics))) {
    commonOmic <- omicNames[1]
    if ( length(omicNames) == 1) {
      omicNames <- rep(omicNames, length(inputOmics))
    }
  } else {
    if (is.null(commonOmic)) {
      commonOmic <- names(which(table(omicNames) == max(table(omicNames)))[1])
    } else {
      if(!is.element(commonOmic, omicNames)) {
        stop("Error: commonOmic is not contained in omicNames.")
      }
    }
    if ( max(table(omicNames)) < length(batches)) {
      stop("Error: There is no at least two omic with the same identification. The common omic is requiered to be named equally" )
    }
  }
  #
  aux.list <- inputOmics[which(omicNames==commonOmic)]
  if (is.list(aux.list)) {
    if ( length( aux.list) > 1) {
      varSpace <- lapply(aux.list, function(l) {
        rownames(l)
      })
      len <- unlist(lapply(varSpace, length))
      init.var <- varSpace[[1]]
      common.var <- init.var
      for ( i in seq_along(varSpace)) {
        aux1 <- varSpace[[i]][is.element(varSpace[[i]], common.var)]
        common.var <- common.var[is.element(common.var, varSpace[[i]])]
      }
      #
      if ( length(common.var) == 0 || is.null(common.var)) {
        stop("Error: Common omics do not share the variable space. Note that rows must be named.")
      }
      #
      if ( length(common.var) != length(init.var)) {
        warning("Variable spaces of common omics have been modified in order to ensure the number of variables (and order) is the same across matrices.")
      }
      #
      aux.list <- lapply(aux.list, function(l) {
        l[common.var,]
      })
    }
  }

  inputOmics[which(omicNames==commonOmic)] <- aux.list

  # Set batch factor at a single vector if needed
  if (dim(as.data.frame(batchFactor))[2] > 1) {
    batchFactor <- apply(batchFactor, 1, paste, collapse = "/")
  }

  # cond.factor input
  if (class(experimentalDesign) != "list") {
    experimentalDesign <- sapply(levels(batchFactor), function(x) {
      data.frame(experimentalDesign)[which(batchFactor == x),]
    }, simplify = F)
  }
  omicList <- sapply(levels(factor(batchFactor)), function(x) {
    aux.list <- inputOmics[which(batchFactor==x)]
    names(aux.list) <- omicNames[which(batchFactor==x)]
    maplist <- lapply(aux.list, function(l) {
      data.frame(primary = paste0(x, data.frame(experimentalDesign[[x]])[,1],
                                  1:length(data.frame(experimentalDesign[[x]])[,1])),
                 colname = colnames(l),
                 stringsAsFactors = FALSE)
    })
    colDat <- data.frame("tfactor" = apply(as.data.frame(experimentalDesign[[x]]),
                                           1, paste, collapse = "/"),
                         row.names = paste0(x, data.frame(experimentalDesign[[x]])[,1],
                                            1:length(data.frame(experimentalDesign[[x]])[,1])))
    sampMap <- MultiAssayExperiment::listToMap(maplist)
    MultiAssayExperiment::MultiAssayExperiment(experiments = aux.list,
                                               colData = colDat,
                                               sampleMap = sampMap)
  })

  retobj <- list("ListOfBatches" = omicList,
                 "commonOmic" = commonOmic)

  # Test colData compatibility -------------------------------------------------
  gcoldata <- unlist(lapply(omicList, function(x) {
    factors <- names(table(as.data.frame(MultiAssayExperiment::colData(x))))
    lapply(omicList, function(l) {
      aux <- names(table(as.data.frame(MultiAssayExperiment::colData(l))))
      if (length(intersect(factors, aux)) != length(factors)) {
        val <- FALSE
      } else {
        val <- TRUE
      }
      val
    })
  }
  ))
  if ( is.element(FALSE, gcoldata)) {
    stop ("Experimental design is not exactly the same across batches. Such condition is mandatory")
  }


  return(mbacClass(retobj))
}
