#' @include ASCAfun1.R
#' @include ASCAfun12.R
#' @include ASCAfunres.R
NULL

#' #' ASCA for 2 experimental factors
#'
#' @param X Data matrix (p x n) that contains expression values of n genes (in columns) and p conditions (in rows).
#' @param Designa Matrix (p x I) that contains 0's and 1's for the time-points in the experiment.
#' @param Designb Matrix (p x J) experimental group.
#' @param Designc Unused argument.
#' @param Fac Vector with numbers of components to extract in each submodel, by order:
#' \enumerate{
#'     \item Model a (time)
#'     \item Model b (second factor)
#'     \item Model ab (interaction) to not consider interations set to 0
#'     \item Number components of residues
#' }
#' @param Interaction Logical. Whether to model the interaction between factors or not.
#' @param Variability From 0 to 1. Minimum percent of data variability that must be explained by each model. By default, 0.75.
#' @param beta Numeric. Components that represent more than beta times the average variability are identified as systematic noise in residuals. Used in noise reduction mode. By default, 2.
#' @return A list containing as many slots as computed models, each one with the following values:
#' \describe{
#'    \item{data}{data matrix used for PCA on that factor, corresponds to the estimated effects for each level of the factor}
#'    \item{scores}{corresponding to the PCA of that factor}
#'    \item{loadings}{corresponding to the PCA of that factor}
#'    \item{var.exp}{explained variability for each component in the PCA o that factor}
#'    \item{X}{the data matrix associated to each model}
#'    \item{TP}{the reconstructed data matrix corresponding that factor}
#'    \item{E}{residuals for the PCA model of that factor}
#'    \item{leverage}{gene leverages for the computed PCA model}
#'    \item{SPE}{gene SPEs for the computed PCA model}
#'  }
#'
#' @examples
#'
#' \dontrun{
#' my.asca <- ASCA.2f(X = t(data.example), Designa = Designa,
#'     Designb = Designb, Designc = NULL, Fac = c(1,2,2,2), type= 1)
#'}
ASCA.2f<-function(X = X, Designa = Designa, Designb = Designb,
                  Designc=NULL, Fac=c(1,2,2),
                  Interaction=TRUE, Variability, beta=2)
{
  #--------------------------------------------------------------------------------------
  #  Dimensions of the matrices:
  #  X (p x n) contains expression values of n genes (in columns) and p conditions (in rows)
  #  Designa (p x I) contains 0's and 1's for the TIME-POINTS in the experiment
  #  Designb (p x J) EXPERIMENTAL GROUP
  #  Designres (p x H) INDIVIDUALS
  #  Interaction = TRUE to consider interaction "ab" in the separated model

  n<-ncol(X)
  p<-nrow(X)
  I<-ncol(Designa)
  J<-ncol(Designb)

  Faca=Fac[1] # number components Model a (time)
  Facbab=Fac[2] # number components Model bab (second factor or second factor plus interaction)
  Facres=Fac[3] # number components Residues

  #----------------------- Calculate Overall Mean --------------------------------------

  offset<-apply(X,2,mean)
  Xoff<-X-(cbind(matrix(1,nrow=p,ncol=1))%*%rbind(offset))


  #-----------------------  PART I: Submodel a (TIME) -----------------------------------

  Model.a<-ASCAfun1(Xoff,Designa,Faca,Variability)
  Xres<-Xoff-Model.a$X

  #-------------------------- PART II: Submodel b and ab -------------------------------
  if (Interaction){
    Model.bab<-ASCAfun12(Xoff,Designa,Designb,Facbab,Variability)
  } else{
    Model.b<-ASCAfun1(Xoff,Designb,Facbab,Variability)
  }
  #browser()
  # ------------------------Collecting models ------------------------------------------

  models <- ls(pattern="Model")
  output <- vector(mode="list")
  Xres <- Xoff
  for (i in 1: length(models)) {
    mymodel <- get(models[i], envir=environment())
    output <- c(output, list(mymodel))
    Xres <- Xres - mymodel$X
    rm(mymodel)
    gc()

  }
  names(output) <- models

  #------------------------- PART III: Submodel res -----------------------------------
  if(Facres!=0){
    Model.res<-ASCAfun.res(Xres,Facres,beta)
  }
  Model.res<-list(Model.res)
  names(Model.res)<-c("Model.res")
  output<-c(output,Model.res)


  #------------------------- Add Input data to the Output ----------------------------

  Input<-list(X, Designa, Designb, Designc, Fac, Interaction,Xoff)
  names(Input)<-c("X", "Designa", "Designb", "Designc", "Fac", "Interaction","Xoff")
  Input<-list(Input)
  names(Input)<-"Input"
  output<-c(output,Input)

  output
}


