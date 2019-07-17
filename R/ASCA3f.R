#' @include ASCA2f.R
NULL

#' ASCA for 3 experimental factors
#'
#' @param X Data matrix (p x n) that contains expression values of n genes (in columns) and p conditions (in rows).
#' @param Designa Matrix (p x I) that contains 0's and 1's for the time-points in the experiment.
#' @param Designb Matrix (p x J) experimental group factor 1.
#' @param Designc Matrix (p x K) another factor.
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
#'
#' @return A list containing as many slots as computed models, each one with the following values:
#' \describe{
#'    \item{data}{data matrix used for PCA on that factor, corresponds to the estimated effects for each level of the factor}
#'    \item{scores}{corresponding to the PCA of that factor}
#'    \item{loadings}{corresponding to the PCA of that factor}
#'    \item{var.exp}{exaplined varaibility for each component in the PCA o that factor}
#'    \item{X}{the data matrix associated to each model}
#'    \item{TP}{the reconstructed data matrix corresponding that factor}
#'    \item{E}{residuals for the PCA model of that factor}
#'    \item{leverage}{gene leverages for the computed PCA model}
#'    \item{SPE}{gene SPEs for the computed PCA model}
#'  }
#'
#' @examples
#' \dontrun{
#' my.asca <- ASCA.3f(X = t(data.example), Designa = Designa,
#'     Designb = Designb, Designc = Designc, Fac= c(1, 2, 2, 2, 2, 2, 2, 2),
#'     type = c(1, 1, 2))
#'  }
ASCA.3f<-function(X = X,Designa = Designa, Designb = Designb, Designc = Designc,Fac = c(1,2,2,2,2,2,2,2),type = c(1,1,2), showvar=TRUE, showscree=TRUE)

{
# Guidance ---------------------------------------------------------------------
#  Dimensions of the matrices:
#  X (p x n) contains expression values of n genes (in columns) and p conditions (in rows)
#  Designa (p x I) contains 0's and 1's for the TIME-POINTS in the experiment
#  Designb (p x J) EXPERIMENTAL GROUP FACTOR 1
#  Designc (p x K) ANOTHER FACTOR

#  type is a vector to indicate if the analyses of the models: (b.ab), (c.ac) and (bc) are studied jointly (1) or separately (2)
#  for model (bc) jointly means that these effects are added to the residuals.


n<-ncol(X)
p<-nrow(X)
I<-ncol(Designa)
J<-ncol(Designb)
K<-ncol(Designc)

Faca=Fac[1]# number components Model a (time)
Facb=Fac[2] # number components Model b  (second factor)
Facc=Fac[3] # number components Model c  (third factor)
Facab=Fac[4] # number components Model ab (interaction) to not consider interations set to 0
Facac=Fac[5] # number components Model ac (interaction) to not consider interations set to 0
Facbc=Fac[6] # number components Model bc (interaction) to not consider interations set to 0
Facabc=Fac[7]  # number components Model abc (interaction) to not consider interations set to 0
Facres=Fac[8]
#----------------------- Calculate Overall Mean -------------------------------------

offset<-apply(X,2,mean)
Xoff<-X-(cbind(matrix(1,nrow=p,ncol=1))%*%rbind(offset))

#-----------------------  PART I: Submodel a (TIME) ---------------------------------

Model.a<-ASCAfun1(Xoff,Designa,Faca)
Xres<-Xoff-Model.a$X

#print("OK PartI")
#-------------------------- PART II.1: Submodel b.ab-----------------------------------

if(type[1]==2)
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

if(type[1]==1)
{
	Model.bab<-ASCAfun12(Xoff,Designa,Designb,Facab)
	# output<-list(Model.a,Model.bab)
	# names(output)<-c("Model.a","Model.bab")
	# Xres<-Xres-Model.bab$X
}
#print("OK PartII.1")
#-------------------------- PART II.2: Submodel (c.ac) -------------------------------
if(type[2]==2)
{
	Model.c<-ASCAfun1(Xoff,Designc,Facc)
	if (Facac != 0 ) {
		Model.ac<-ASCAfun2(Xoff,Designa,Designc,Facac)
	}
	else {
		# Model.ac = 0
	}
	# LIST<-list(Model.c,Model.ac)
	# names(LIST)<-c("Model.c","Model.ac")
	# output<-c(output,LIST)

	# Xres<-Xres-Model.c$X-Model.ac$X
}

if(type[2]==1)
{
	Model.cac<-ASCAfun12(Xoff,Designa,Designc,Facac)

	# LIST<-list(Model.cac)
	# names(LIST)<-c("Model.cac")
	# output<-c(output,LIST)

	# Xres<-Xres-Model.cac$X
}
#print("OK PartII.2")
#-------------------------- PART II.3: Submodel (bc) --------------------------------
if(type[3]==2)
{
	if (Facbc != 0 ) {
		Model.bc<-ASCAfun2(Xoff,Designb,Designc,Facbc)
	}
	else {
		# Model.bc = 0
	}
	# LIST<-list(Model.bc)
	# names(LIST)<-c("Model.bc")
	# output<-c(output,LIST)

	# Xres<-Xres-Model.bc$X
}
#print("OK PartII.3")

#-------------------------- PART II.4: Submodel (abc) --------------------------------

	if (Facabc != 0 ) {
		Model.abc<-ASCAfun.triple(Xoff,Designa,Designb,Designc,Facabc)
	}
	else {
		# Model.abc = 0
	}
	# LIST<-list(Model.abc)
	# names(LIST)<-c("Model.abc")
	# output<-c(output,LIST)

	# Xres<-Xres-Model.abc$X
	#print("OK PartII.4")

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

#------------------------- PART III: Submodel res -----------------------------------

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

