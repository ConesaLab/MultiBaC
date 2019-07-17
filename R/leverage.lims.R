#' @include SPE.lims.R
NULL

#' Find leverage (importance of a variable in the PCA model) cut-off
#'
#' @param data Data matrix (p x n) that contains expression values of p genes (in rows) and n conditions (in columns).
#' @param R Number of repeats.
#' @param FUN Function to be called (ASCA.2f or ASCA.3f)
#' @param Designa First design matrix to pass to FUN.
#' @param Designb Second design matrix to pass to FUN.
#' @param Designc Third design matrix to pass to FUN.
#' @param Fac Vector with numbers of components to pass to FUN.
#' @param alpha Alpha value
#' @param type Vector indicating whether the analyses of the model of the factors that interact with time (b.ab and c.ac) are studied jointly(1) or separately(2).
#' @param showvar Show variability associated to each model or not.
#' @param showscree Show a screenplot.
#'
#' @return A list with elements:
#' \describe{
#'     \item{NullDistribution}{}
#'     \item{Cutoff}{}
#'     \item{Selection}{}
#' }
#'
#' @examples
#' \dontrun{
#' leverage.cutoff <- leverage.lims(data = data.example, R = 3, ASCA.3f, Designa = Designa,
#'     Designb = Designb, Designc = Designc,Fac = c(1,2,2,2,2,2,2,2), type = c(1,1,2),
#'     alpha = 0.95, showvar = FALSE, showscree = FALSE)
#'}
leverage.lims <- function(data = data, R = 100, FUN, Designa = Designa, Designb = Designb, Designc = NULL, Fac = c(1,2,2,2), type = 2, alpha = 0.01, showvar=FALSE, showscree=FALSE)
{
## Compute ASCA model for data
	Model <- FUN(X = t(data), Designa = Designa,Designb = Designb,Designc = Designc, Fac = Fac,type = type)
	n <- ncol(data)
	lim <- Selection <- vector(mode = "list", length = length(Model)-1)
	names(lim) <- names(Selection) <- names(Model)[1:length(lim)]

## Calculate the reference distribution of leverages
	for (i in 1:R) {
		place <- sample(1:n)
		permu <- data[,place]
		A <- FUN (X = t(permu), Designa = Designa, Designb = Designb, Designc = Designc, Fac = Fac, type = type, showvar=showvar, showscree=showscree)
		for (j in 1:length(lim)) {
			lim[[j]] <- cbind(lim[[j]],A[[j]]$leverage)
		}
	}
		rm(A)
		gc()
## Find alpha quantile value for reference distribution
	QC <- apply(sapply(lim,function(x) {apply(x,1,quantile,probs=1 - alpha)}),2,quantile,probs=1 - alpha)

# Gene selection
	for (h in 1:length(lim)) {
		Selection[[h]] <- rownames(data)[Model[[h]]$leverage > QC[[h]]]
	}
	output <- list(lim, QC, Selection)
	names(output) <- c("NullDistribution", "Cutoff", "Selection")
	output
}
