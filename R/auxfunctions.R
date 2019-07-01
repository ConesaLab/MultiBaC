#' @import stats
#' @import ggplot2
#' @include ASCA3f.R
NULL

#' batchEstimation
#'
#' @param Xunfold Samples x genes matrix containing the common omic information for all batches (labs)
#' @param batch_factor Class factor indicating the batch of each sample in Xunfold
#'
#' @return
#' @export
#'
#' @examples
batchEstimation <- function(ListOfBatches, batch_factor, commonOmic) {

  initpar <- par(c("mai", "pin", "xpd"))
  mai <- c(initpar$mai[1:3], initpar$pin[2]-0.15)
  on.exit(par(mai = initpar$mai, xpd = initpar$xpd))
  par(mai = mai)

  # Extract common omic ---------------------------------------------------------------
  Xunfold <- NULL
  for ( lab in seq_along(ListOfBatches)) {
    Xunfold <- rbind(Xunfold, t(ListOfBatches[[lab]][[commonOmic]]))
  }

  #----------------------- Calculate Overall Mean -------------------------------------

  offset<-apply(Xunfold,2,mean)
  Xoff<-Xunfold-(cbind(matrix(1,nrow=nrow(Xunfold),ncol=1))%*%rbind(offset))

  #----------------------- Compute coefficients --------------------------------------

  design.matrix <- model.matrix(~ factor(batch_factor))
  #
  design.matrix[which(design.matrix==0)] <- -1
  #
  beta.hat <- solve(t(design.matrix)%*%design.matrix) %*% t(design.matrix) %*% Xoff



  #----------------------- Theoretical distributions ----------------------------------

  beta.hat_v <- as.vector(t(beta.hat))

  reg_v <- c()
  for ( i in 1:(length(beta.hat_v)/dim(beta.hat)[2])) {
    reg_v <- c(reg_v, rep(i,dim(beta.hat)[2]))
  }

  toplot <- data.frame("Coef" = factor(reg_v[(dim(beta.hat)[2]+1):length(beta.hat_v)]),
                       "Values" = beta.hat_v[(dim(beta.hat)[2]+1):length(beta.hat_v)])


  p <- ggplot(toplot, aes(x=Coef, y=Values)) +
    geom_violin() +
    labs(title="Batch Magnitude Plot", x = "Batch: 1 as reference", y = "Coefficients per Gene", size=2) +
    scale_fill_brewer(palette="Blues") +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
    theme_classic() +
    theme(axis.text=element_text(size=18),
          title=element_text(size=16,face="bold"),
          legend.text = element_text(size=20),
          legend.spacing = unit(.25,"cm")) +
    geom_hline(aes(yintercept=6, col = "High"), lty = 5, show.legend = TRUE) +
    geom_hline(aes(yintercept=-6, col = "High"), lty = 5,show.legend = TRUE) +
    geom_hline(aes(yintercept=4, col = "Medium"), lty = 5,show.legend = TRUE) +
    geom_hline(aes(yintercept=-4, col = "Medium"), lty = 5,show.legend = TRUE) +
    geom_hline(aes(yintercept=2, col = "Low"), lty = 5,show.legend = TRUE) +
    geom_hline(aes(yintercept=-2, col = "Low"), lty = 5,show.legend = TRUE) +
    labs(col = "Theorical
Batch
Magnitudes")

  return(p)

}


#' new_module
#'
#' @param YZ List of non-common omics
#' @param cond.factor List of experimental condition factor type objects for each matrix in YZ
#' @param missingOmics List of predicted omics for each batch (lab)
#'
#' @return
#' @export
#'
#' @examples
#'
new_module <- function (YZ, missingOmics, cond.factor) {

  omic_Sfactor <- c()
  cond_Sfactor <- c()
  batch_Sfactor <- c()
  #
  # omic_Pfactor <- c()
  # batch_Pfactor <- c()
  number_of_models <- seq_along(YZ)
  for ( i in seq_along(YZ) ) {
    omic_Sfactor <- c(omic_Sfactor, rep(i, nrow(YZ[[i]])))
    cond_Sfactor <- c(cond_Sfactor, cond.factor[[i]])
    batch_Sfactor <- c(batch_Sfactor, rep(i, nrow(YZ[[i]])))
    #
    # omic_Pfactor <- c(omic_Pfactor, rep(i, ncol(YZ[[i]])))
    # batch_Pfactor <- c(batch_Pfactor, rep(i, ncol(YZ[[i]])))
    #
    sel <- number_of_models[!is.element(number_of_models, i)]
    for ( j in seq_along(sel) ) {
      batch_Sfactor <- c(batch_Sfactor, rep(sel[j], nrow(YZ[[sel[j]]])))
      cond_Sfactor <- c(cond_Sfactor, cond.factor[[sel[j]]])
      omic_Sfactor <- c(omic_Sfactor, rep(i, nrow(YZ[[sel[j]]])))
      #
      # omic_Pfactor <- c(omic_Pfactor, rep(i, ncol(YZ[[sel[j]]])))
      # batch_Pfactor <- c(batch_Pfactor, rep(sel[j], ncol(YZ[[sel[j]]])))
    }

  }
  # Model matrices
  cond_Smatrix <- model.matrix(~ 0 + factor(cond_Sfactor))
  omic_Smatrix <- model.matrix(~ 0 + factor(omic_Sfactor))
  batch_Smatrix <- model.matrix(~ 0 + factor(batch_Sfactor))
  #
  # omic_Pmatrix <- model.matrix(~ 0 + factor(omic_Pfactor))
  # batch_Pmatrix <- model.matrix(~ 0 + factor(batch_Pfactor))

  # Select number of components
  ncomp <- min(unlist(lapply(YZ, function(x) min(dim(x)))))

  # Create PCA models
  PCAlist <- list()
  number_of_models <- seq_along(YZ)
  index <- 1
  for ( m in seq_along(YZ) ) {
    sel <- number_of_models[!is.element(number_of_models, m)]
    PCAlist[[index]] <- nipals(YZ[[m]], a = min(ncomp,5))
    aux <- list()
    for ( i in seq_along(sel) ) {
      aux[[i]] <- nipals(missingOmics[[sel[i]]][[m]], a = min(ncomp,5))
    }
    PCAlist <- c(PCAlist, aux)
    index <- index + length(aux) + 1
  }
  ## set names
  pcaNames <- c()
  for ( i in seq_along(YZ) ) {
    pcaNames <- c(pcaNames, names(YZ)[i],
                  unlist(lapply(1:(length(YZ)-1),
                         function(x) {paste0(names(YZ)[i], "fake", x)})))
  }
  names(PCAlist) <- pcaNames

  # create big matrix
  big_Smatrix <- PCAlist[[1]][[1]]
  for ( i in 2:length(PCAlist) ) {
    # check rotation
    ref <- mean(PCAlist[[1]][[1]][cond_Sfactor[1:dim(YZ[[1]])[1]],2])
    test <- mean(PCAlist[[i]][[1]][cond_Sfactor[1:dim(PCAlist[[i]][[1]])[1]],2])
    if ( sign(ref) == sign(test) ) c = 1 else c = -1
    big_Smatrix <- rbind(big_Smatrix, c*PCAlist[[i]][[1]])
    PCAlist[[i]][[2]] <- c*PCAlist[[i]][[2]]
  }

  # big_Pmatrix <- PCAlist[[1]][[2]]
  # for ( i in 2:length(PCAlist) ) {
  #   # check rotation
  #   ref <- mean(PCAlist[[1]][[2]][cond_Sfactor[1:dim(YZ[[1]])[2]],2])
  #   big_Pmatrix <- rbind(big_Pmatrix, PCAlist[[i]][[2]])
  # }

  arsyn <- ASCA.3f(big_Smatrix, Designa = cond_Smatrix, Designb = omic_Smatrix, Designc = batch_Smatrix,
                   Fac = c(1,1,1,1,1,1,1,1) , type = c(2,2,2), showvar = FALSE, showscree = FALSE)

  bigM.star <- big_Smatrix - arsyn$Model.c$TP - arsyn$Model.ac$TP

  # arsyn <- ASCA.2f(big_Pmatrix, Designa = omic_Pmatrix, Designb = batch_Pmatrix,
  #                  Fac = c(1,1,0,0) , type = 2, showvar = FALSE, showscree = FALSE)
  #
  # bigM.Pstar <- big_Pmatrix - arsyn$Model.b$TP #- arsyn$Model.bc$TP

  YZ.star <- list()
  roff <- 0; index <- 1
  for ( i in seq_along(PCAlist) ) {
    if ( is.element(names(PCAlist)[i], names(YZ)) ) {
      omic.star <- bigM.star[(roff+1):(roff+nrow(PCAlist[[i]][[1]])),]%*%t(PCAlist[[i]][[2]])
      roff <- roff + nrow(PCAlist[[i]][[1]])
      YZ.star[[index]] <- omic.star
      index <- index + 1
    } else {
      roff <- roff + nrow(PCAlist[[i]][[1]])
    }
  }
  names(YZ.star) <- names(YZ)

  return(YZ.star)

}


#' nipals
#'
#' @param X Samples x genes numerical matrix.
#' @param a Number of latent variables (or components)
#' @param it Maximum number of iterations allowed to nipals algorithm
#' @param tol threshold value to validate nipals algorithm convergence
#' @param scale TRUE or FALSE. Whether X and Y matrices have to be scaled
#' @param center TRUE or FALSE. Whether X and Y matrices have to be centered
#'
#' @return A list containing two slots:
#' \describe{
#'     \item{T}{scores}
#'     \item{P}{loadings}
#' }
#' @export
#'
#' @examples
#' pca <- nipals (X = t(assayData(my_counts)$exprs), a = 2, center = TRUE )
nipals <- function (X, a, it = 50, tol = 1e-10, center = FALSE, scale = FALSE)
{
  Xh <- scale(X, center = center, scale = scale)
  nr <- 0
  T <- NULL
  P <- NULL
  for (h in seq_len(a)) {
    th <- Xh[, 1]
    ende <- FALSE
    while (!ende) {
      nr <- nr + 1
      ph <- t((t(th) %*% Xh) * as.vector(1/(t(th) %*% th)))
      ph <- ph * as.vector(1/sqrt(t(ph) %*% ph))
      thnew <- t(t(ph) %*% t(Xh) * as.vector(1/(t(ph) %*% ph)))
      prec <- t(th - thnew) %*% (th - thnew)
      th <- thnew
      if (prec <= (tol)) {
        ende <- TRUE
      }
    }
    Xh <- Xh - (th %*% t(ph))
    T <- cbind(T, th)
    P <- cbind(P, ph)
    nr <- 0
  }
  list(T = T, P = P)
}


#' Title
#'
#' @param list
#'
#' @return
#' @export
#'
#' @examples
inputeval <- function(arguments) {
  # ListOfBatches ------------------------------------------------------
  #if (!all(unlist(lapply(list$ListOfBatches, class))=="MultiAssayExperiment")) {
    if (is.null(names(arguments$ListOfBatches))) {
      stop ("Stop at input evaluation: Elements in ListOfBatches must be named")
    }
  #}

}


#' optimizeComponents
#'
#' @param modelObject
#'
#' @return Instance of class "pls" with a new slot "optimal.ncomp"
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' }
optimizeComponents = function (modelObject) {

  # SYMPCA-P rule
  opt <- 0
  lim <- 0

  for ( i in seq_along(modelObject@modelDF$`Q2(cum)`) ) {
    if ( (modelObject@modelDF$`Q2(cum)`[i] - lim) > 0.0975 ) {
      opt <- i
      lim <- modelObject@modelDF$`Q2(cum)`[i]
    }
  }

  if ( modelObject@modelDF$`Q2(cum)`[opt] < 0.5 ) {
    warning("Predictian capacity is not good, it may affect correction performance")
  }

  opt
}
