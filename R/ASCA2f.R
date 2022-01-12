#' @include ASCAfun1.R
#' @include ASCAfun12.R
#' @include ASCAfunres.R
NULL

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
  
  if (Interaction) {
    Fac <- c(ncol(Designa)-1,ncol(Designb)*ncol(Designa)-1,2)
    names(Fac) <- c("Model.a", "Model.bab", "Model.res")
  } else {
    Fac <- c(ncol(Designa)-1,ncol(Designb)-1,2)
    names(Fac) <- c("Model.a", "Model.b", "Model.res")
  }


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
    Model.res<-list(Model.res)
    names(Model.res)<-c("Model.res")
    output<-c(output,Model.res)
  }


  #------------------------- Add Input data to the Output ----------------------------

  Input<-list(X, Designa, Designb, Designc, Fac, Interaction,Xoff)
  names(Input)<-c("X", "Designa", "Designb", "Designc", "Fac", "Interaction","Xoff")
  Input<-list(Input)
  names(Input)<-"Input"
  output<-c(output,Input)

  output
}


