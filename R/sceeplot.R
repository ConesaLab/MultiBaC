#' @include PCA-GENES.R
NULL

#' Show a scree plot of an ASCA analysis
#'
#' @param asca Result of ASCA analysis.
#'
#' @export
#'
#' @examples
#' my.asca <- ASCA.3f(X = t(data.example), Designa = Designa,
#'     Designb = Designb, Designc = Designc, Fac= c(1, 2, 2, 2, 2, 2, 2, 2),
#'     type = c(1, 1, 2))
#'
#' screeplot(my.asca)
screeplot <- function (asca) {

num <- length(asca)
r <- round(sqrt(num))
c <- num-r
n <- r*c

initpar <- par("mfrow")
on.exit(par(mfrow = initpar))
par(mfrow = c(1,num))

#layout(matrix(c(1:num), 1, num, byrow = TRUE))

   for ( i in 1: length (asca)) {

      plot(asca[[i]]$var.exp[,1], type="l", main=paste("Screeplot", names(asca)[[i]]), xlab="Component", ylab="Expalined variability")

   }

}
