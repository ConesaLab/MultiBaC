#' @include sceeplot.R
NULL

show.var <- function (asca, Xoff=Xoff) {

   ev <- NULL
   a <- length(asca)
   a <- a-1
   for ( i in 1: a)  {
      ev <- c(ev,sum(asca[[i]]$X^2)/sum(Xoff^2))
   }
   ev.t <- c(ev,1-sum(ev))
   names (ev.t) <- names(asca)
   ev.t
}
