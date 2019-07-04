#' @import graphics
#' @import stats
#' @include ASCAfun1.R
#' @include ASCAfun2.R
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
#' @param type Number (1 or 2) to indicate if the analyses of the model b.ab is studied jointly(1) or separately(2).
#' @param showvar Show variability associated to each model or not.
#' @param showscree Show a screenplot.
#'
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
#' @export
#'
#' @examples
#'
#' \dontrun{my.asca <- ASCA.2f(X = t(data.example), Designa = Designa,
#'     Designb = Designb, Designc = NULL, Fac = c(1,2,2,2), type= 1)}
ASCA.2f<-function(X = X, Designa = Designa, Designb = Designb, Designc = NULL,
                  Fac = c(1,2,2,2), type = 1, showvar = TRUE, showscree = TRUE) {
#--------------------------------------------------------------------------------------
#  Dimensions of the matrices:
#  X (p x n) contains expression values of n genes (in columns) and p conditions (in rows)
#  Designa (p x I) contains 0's and 1's for the TIME-POINTS in the experiment
#  Designb (p x J) EXPERIMENTAL GROUP
#  Designres (p x H) INDIVIDUALS

#  type is a number (1 or 2) to indicate if the analyses of the model b.ab is studied jointly(1) or separately(2)

n<-ncol(X)
p<-nrow(X)
I<-ncol(Designa)
J<-ncol(Designb)

Faca=Fac[1] # number components Model a (time)
Facb=Fac[2] # number components Model b  (second factor)
Facab=Fac[3] # number components Model ab (interaction) to not consider interations set to 0
Facres=Fac[4] # number components Residues

#----------------------- Calculate Overall Mean --------------------------------------

offset<-apply(X,2,mean)
Xoff<-X-(cbind(matrix(1,nrow=p,ncol=1))%*%rbind(offset))


#-----------------------  PART I: Submodel a (TIME) -----------------------------------

Model.a<-ASCAfun1(Xoff,Designa,Faca)
Xres<-Xoff-Model.a$X

#print("OK PartI")

#-------------------------- PART II: Submodel b and ab -------------------------------

if(type==2)
{
	Model.b<-ASCAfun1(Xoff,Designb,Facb)
	if (Facab != 0 ) {
		Model.ab<-ASCAfun2(Xoff,Designa,Designb,Facab)
	}
	else {
		# Model.ab = 0
	}
	# output<-list(Model.a,Model.b,Model.ab)
	# names(output)<-c("Model.a","Model.b","Model.ab")

	# Xres<-Xres-Model.b$X-Model.ab$X
}

if(type==1)
{
	Model.bab<-ASCAfun12(Xoff,Designa,Designb,Facab)

	# output<-list(Model.a,Model.bab)
	# names(output)<-c("Model.a","Model.bab")

	#Xres<-Xres-Model.bab$X
}

#print("OK PartII")


# ------------------------Collecting models ------------------------------------------

models <- ls(pattern="Model")
output <- vector(mode="list")
Xres <- Xoff
for (i in seq_along(models)) {
	mymodel <- get(models[i], envir=environment())
	output <- c(output, list(mymodel))
	Xres <- Xres - mymodel$X
	rm(mymodel)
	gc()

}
names(output) <- models

#------------------------- PART III: Submodel abg -----------------------------------

Model.res<-ASCAfun.res(Xres,Facres)

#print("OK PartIII")

LIST<-list(Model.res)
names(LIST)<-c("Model.res")
output<-c(output,LIST)

#------------------------- Show varaibility associated to each model-----------------
if (showvar) {
	print("Variability associated to each model:")
	ev <- show.var(output, Xoff=Xoff)
	print(ev)
}
#------------------------- Show screeplot ---------------------------------------
if (showscree) {
	screeplot(output)
}
output
}

