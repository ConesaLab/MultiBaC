#' @include show.var.R
NULL

#' Find squared prediction error cut-off
#'
#' @param my.asca Result of ASCA analysis.
#' @param alpha Alpha value.
#'
#'
#' @examples
#' \dontrun{
#' my.asca <- ASCA.3f(X = t(data.example), Designa = Designa,
#'     Designb = Designb, Designc = Designc, Fac= c(1, 2, 2, 2, 2, 2, 2, 2),
#'     type = c(1, 1, 2))
#'
#'  SPE.cutoff <- SPE.lims(my.asca, alpha = 0.01)
#'  }
SPE.lims <- function (my.asca, alpha)
{
limits <- vector("list", length(my.asca))
names(limits) <- names(my.asca)
	for (i in 1:length(my.asca)) {
		assign ("model", my.asca[[i]])
		if (!is.null(model$SPE)) {
			m <- mean(model$SPE)
			v <- var(model$SPE)
			g <- v/(2*m)
			h <- 2*m*m/v
			limits[[i]] <- g*qchisq(1-alpha, df=h)
		}
	}
limits
}

